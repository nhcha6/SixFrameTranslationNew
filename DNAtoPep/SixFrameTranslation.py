from DNA_CODON_TABLE import *
from Bio import SeqIO
from Bio.Seq import Seq
import psutil
from Bio.SeqRecord import SeqRecord
import multiprocessing
from multiprocessing import Queue
import tempfile
from removeSubsetSeq import *
from operator import itemgetter
import heapq
from time import time
from queue import Queue
import logging
import os


MEMORY_THRESHOLD = 0.2


def memory_usage_psutil():
    # return the memory usage in percentage like top
    process = psutil.Process(os.getpid())
    mem = process.memory_percent()
    return mem

# set of new function which don't require the storage of all forward and reverse frames to run
def buildForwProt(seq, minLen):

    proteins = []
    for j in range(0,3):
        proteinTemp = ""
        remainder = (len(seq[j:-1]) + 1) % 3
        for i in range(j, len(seq) - remainder, 3):
            compressed = seq[i:i+2]
            if compressed in COMPRESSED_TABLE.keys():
                amino = COMPRESSED_TABLE[compressed]
            else:
                codon = seq[i:i+3]
                if 'X' in codon:
                    amino = 'X'
                else:
                    amino = DNA_TABLE[codon]
            if amino == -1:
                if len(proteinTemp) >= minLen:
                    proteins.append(proteinTemp)
                proteinTemp = ""
            elif i == len(seq) - remainder - 3:
                proteinTemp += amino
                if len(proteinTemp) >= minLen:
                    proteins.append(proteinTemp)
                proteinTemp = ""
            else:
                proteinTemp += amino
    return proteins

def buildRevProt(seq, minLen):
    proteins = []
    for j in range(0,3):
        proteinTemp = ""
        remainder = (len(seq[j:-1]) + 1) % 3
        for i in range(j, len(seq) - remainder, 3):
            if i == 0:
                codon = createReverseSeq(seq[-3:])
                compressed = codon[0:2]
            else:
                codon = createReverseSeq(seq[-1*(i+3):-i])
                compressed = codon[0:2]

            if compressed in COMPRESSED_TABLE.keys():
                amino = COMPRESSED_TABLE[compressed]
            else:
                if 'X' in codon:
                    amino = 'X'
                else:
                    amino = DNA_TABLE[codon]
            if amino == -1:
                if len(proteinTemp) >= minLen:
                    proteins.append(proteinTemp)
                proteinTemp = ""
            elif i == len(seq) - remainder - 3:
                proteinTemp += amino
                if len(proteinTemp) >= minLen:
                    proteins.append(proteinTemp)
                proteinTemp = ""
            else:
                proteinTemp += amino
    return proteins

def seqToProteinNew(dnaSeq, minLen, name):

    # NEED TO COUNT HOW MANY N'S At start and end, and remove them.
    newSeq = str(dnaSeq).upper()
    newSeq = removeNsDNA(newSeq)

    newSeq = newSeq.upper().replace('N', 'X')

    start = time()

    proteins = buildForwProt(newSeq, minLen) + buildRevProt(newSeq, minLen)
    seqToProteinNew.toWriteQueue.put([name, proteins])
    end = time()
    print(str(dnaSeq[0:5]) + "took " + str(end-start))


def removeNsDNA(dnaSeq):

    fiveFrameCount = 0

    for base in dnaSeq:
        if base == 'N':
            fiveFrameCount +=1
        else:
            break

    reverseSeq = dnaSeq[::-1]

    threeFrameCount = 0
    for base in reverseSeq:
        if base == 'N':
            threeFrameCount +=1
        else:
            break

    return dnaSeq[fiveFrameCount: len(dnaSeq)-threeFrameCount]

def generateOutputNew(outputPath, minLen, input_path, removeSubFlag, writeSubFlag, originFlag):



    start = time()
    num_workers = multiprocessing.cpu_count()
    toWriteQueue = multiprocessing.Queue()
    writerProcess = multiprocessing.Process(target=writer, args=(toWriteQueue, outputPath, removeSubFlag, writeSubFlag, originFlag))
    writerProcess.start()



    pool = multiprocessing.Pool(processes=num_workers, initializer=poolInitialiser,
                                    initargs=(toWriteQueue,))

    with open(input_path, "rU") as handle:
        counter = 0
        for record in SeqIO.parse(handle, 'fasta'):
            counter += 1
            name = "rec" + str(counter)
            dnaSeq = record.seq

            print("Starting process for " + str(dnaSeq[0:5]))

            pool.apply_async(seqToProteinNew, args=(dnaSeq, minLen, name))
            #proteis = seqToProteinNew(dnaSeq, minLen)
            #toWriteQueue.put([name, proteins])
    pool.close()
    pool.join()

    toWriteQueue.put('stop')
    writerProcess.join()
    end = time()
    print("Altogether took " + str(end-start))



def writer(queue, outputPath, removeSubFlag, writeSubFlag, originFlag):
    seenProteins = {}
    saveHandle = outputPath + '-All.fasta'
    sortedTempFileNames = Queue()
    iterTempFileNames = []

    with open(saveHandle, "w") as output_handle:
        while True:
            tuple = queue.get()
            if tuple == 'stop':
                print("All proteins added to writer queue")
                break
            proteins = tuple[1]
            name = tuple[0]
            for protein in proteins:
                if protein not in seenProteins.keys():
                    seenProteins[protein] = [name]
                else:
                    if name not in seenProteins[protein]:
                        seenProteins[protein].append(name)

            if memory_usage_psutil() > MEMORY_THRESHOLD:

                sortedSeenProts = sorted([*seenProteins], key=len, reverse=True)
                sortedTempName = writeTempFasta(sortedSeenProts)
                sortedTempFileNames.put(sortedTempName)

                iterTempName = writeTempFasta(seenProteins)
                iterTempFileNames.append(iterTempName)

                seenProteins = {}
        if sortedTempFileNames.empty():
            # come back to this!!!!
            print('something')
        else:
            sortedPath = mergeSortedFiles(sortedTempFileNames)
            print("Sorted path is :" + sortedPath)
        # print("writing to fasta")
        # SeqIO.write(createSeqObj(seenProteins), output_handle, "fasta")

    if removeSubFlag:
        print('removing subset sequences')        # writeTempFasta(sortedSeenProts)
        refinedRemoveSubsetSeq(originFlag, writeSubFlag, sortedPath, iterTempFileNames, outputPath)
        os.remove(sortedPath)

def combineAllTempFasta(outputTempFiles, writeSubsets=False):

    while not outputTempFiles.empty():

        fileOne = outputTempFiles.get()
        fileTwo = outputTempFiles.get()

        if outputTempFiles.empty():
            break

        seenPeptides = combineTempFile(fileOne, fileTwo, writeSubsets)

        tempName = writeTempFasta(seenPeptides)
        outputTempFiles.put(tempName)

    finalSeenPeptides = combineTempFile(fileOne, fileTwo, writeSubsets)

    # Return the last combination of two files remaining
    return finalSeenPeptides

def combineTempFile(fileOne, fileTwo, writeSubsets=False):
    logging.info("Combining two files !")

    with open(fileOne, 'rU') as handle:
        if writeSubsets:
            seenPeptides = set()
            for line in fileOne:
                seenPeptides.add(line)
        else:
            seenPeptides = {}
            for record in SeqIO.parse(handle, 'fasta'):

                peptide = str(record.seq)
                protein = str(record.name)
                if peptide not in seenPeptides.keys():
                    seenPeptides[peptide] = [protein]
                else:
                    seenPeptides[peptide].append(protein)
    with open(fileTwo, 'rU') as handle:
        if writeSubsets:
            seenPeptides = set()
            for line in fileOne:
                seenPeptides.add(line)
        else:
            for record in SeqIO.parse(handle, 'fasta'):

                peptide = str(record.seq)
                protein = str(record.name)
                if peptide not in seenPeptides.keys():
                    seenPeptides[peptide] = [protein]
                else:
                    seenPeptides[peptide].append(protein)
    # Delete temp files as they are used up
    os.remove(fileOne)
    os.remove(fileTwo)
    return seenPeptides


def writeTempFasta(seenProteins):
    logging.info("Writing to temp")
    temp = tempfile.NamedTemporaryFile(mode='w+t', delete=False)
    try:
        for key, value in seenProteins.items():
            temp.writelines(">")
            for protein in value:
                temp.writelines(str(protein))
            temp.writelines("\n")
            temp.writelines(str(key))
            temp.writelines("\n")
    except AttributeError:

        for protein in seenProteins:
            temp.writelines(str(protein))
            temp.writelines("\n")

    print(temp.name)
    return temp.name

def mergeSortedFiles(tempSortedFiles):

    while not tempSortedFiles.empty():
        fileOne = tempSortedFiles.get()
        fileTwo = tempSortedFiles.get()
        if tempSortedFiles.empty():
            break
        tempName = combineTwoSorted(fileOne, fileTwo)
        tempSortedFiles.put(tempName)
        os.remove(fileOne)
        os.remove(fileTwo)
    finalTempName = combineTwoSorted(fileOne, fileTwo)
    os.remove(fileOne)
    os.remove(fileTwo)
    return finalTempName

def combineTwoSorted(fileOne, fileTwo):
    with open(fileOne) as f1, open(fileTwo) as f2:
        sources = [f1, f2]
        temp = tempfile.NamedTemporaryFile(mode='w+t', suffix='nodups.txt', delete=False)

        decorated = [
            ((len(line), line) for line in f)
            for f in sources]
        merged = heapq.merge(*decorated, reverse=True)

        undecorated = list(map(itemgetter(-1), merged))
        for i in range(len(undecorated) - 1, 0, -1):
            if undecorated[i] == undecorated[i - 1]:
                del undecorated[i]

        temp.writelines(undecorated)
    return temp.name

def poolInitialiser(toWriteQueue):
    seqToProteinNew.toWriteQueue = toWriteQueue


def createReverseSeq(dnaSeq):
    reverseDir = dnaSeq[::-1]
    _tab = str.maketrans(dict(zip('ATCG', 'TAGC')))
    reverseSeq = reverseDir.translate(_tab)
    return reverseSeq


#
# temp = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fasta", delete=False)
#                     allTempFiles.append(temp.name)
#                     counter = 0
#
#                 temp.writelines(">" + record.description)
#                 temp.writelines("\n")
#                 temp.writelines(record.seq)
#                 temp.writelines("\n")
#                 counter += 1