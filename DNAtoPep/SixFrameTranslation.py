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
import time
from queue import Queue
import logging
import os


MEMORY_THRESHOLD = 19
MEMORY_THRESHOLD_LOWER = 24


def memory_usage_psutil():
    # return the memory usage in percentage like top
    mem = psutil.virtual_memory()
    return mem.percent

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

    start = time.time()

    proteins = buildForwProt(newSeq, minLen) + buildRevProt(newSeq, minLen)
    seqToProteinNew.toWriteQueue.put([name, proteins])
    end = time.time()
    #print(str(dnaSeq[0:5]) + "took " + str(end-start))


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

def generateOutputNew(outputPath, minLen, input_path, writeSubFlag, originFlag):



    start = time.time()
    num_workers = multiprocessing.cpu_count()
    toWriteQueue = multiprocessing.Queue()
    writerProcess = multiprocessing.Process(target=writer, args=(toWriteQueue, outputPath, writeSubFlag, originFlag))
    writerProcess.start()



    pool = multiprocessing.Pool(processes=num_workers, initializer=poolInitialiser,
                                    initargs=(toWriteQueue,))

    with open(input_path, "rU") as handle:
        counter = 0
        for record in SeqIO.parse(handle, 'fasta'):
            counter += 1
            name = "rec" + str(counter) + ';'
            dnaSeq = record.seq

            print("Starting process for " + str(dnaSeq[0:5]))

            pool.apply_async(seqToProteinNew, args=(dnaSeq, minLen, name))
            if memory_usage_psutil() > MEMORY_THRESHOLD:
                print('Memory usage exceded. Waiting for processes to finish.')
                pool.close()
                pool.join()
                toWriteQueue.put("memFlag")
                break
                pool = multiprocessing.Pool(processes=num_workers, initializer=poolInitialiser,
                                            initargs=(toWriteQueue,))

                # print('in memory thresold process generation')
                # while memory_usage_psutil() > MEMORY_THRESHOLD_LOWER:
                #     time.sleep(5)
                # print('hit lower threshold process generation.')

            #proteis = seqToProteinNew(dnaSeq, minLen)
            #toWriteQueue.put([name, proteins])
    #pool.close()
    #pool.join()

    toWriteQueue.put('stop')
    writerProcess.join()
    end = time.time()
    print("Altogether took " + str(end-start))



def writer(queue, outputPath, writeSubFlag, originFlag):
    seenProteins = {}

    # two sets of tempfiles, the sorted temp files and the seenPeptides temp files
    sortedTempFileNames = Queue()
    iterTempFileNames = []


    while True:
        tuple = queue.get()
        if tuple == 'stop':
            print("All proteins added to writer queue")
            break
        if tuple == "memFlag":
            # sort the proteins by length, and write to a sorted tempFile
            print("memFlag recieved, all processes finished.")
            sortedSeenProts = sorted([*seenProteins], key=len, reverse=True)
            sortedTempName = writeTempFasta(sortedSeenProts)
            sortedTempFileNames.put(sortedTempName)
            # write the current seen proteins straight to temp files
            iterTempName = writeTempFasta(seenProteins)
            iterTempFileNames.append(iterTempName)
            seenProteins = {}
        else:
            proteins = tuple[1]
            name = tuple[0]
            for protein in proteins:
                if protein not in seenProteins.keys():
                    seenProteins[protein] = [name]
                else:
                    if name not in seenProteins[protein]:
                        seenProteins[protein].append(name)

        # # Ran over memory write to temp files
        # if memory_usage_psutil() > MEMORY_THRESHOLD:
        #
        #     # sort the proteins by length, and write to a sorted tempFile
        #     sortedSeenProts = sorted([*seenProteins], key=len, reverse=True)
        #     sortedTempName = writeTempFasta(sortedSeenProts)
        #     sortedTempFileNames.put(sortedTempName)
        #     # write the current seen proteins straight to temp files
        #     iterTempName = writeTempFasta(seenProteins)
        #     iterTempFileNames.append(iterTempName)
        #     seenProteins = {}
    # Make sure to write the sorted files (and iter temp files) if didn't run out of memory initially.
    if sortedTempFileNames.empty():
        sortedSeenProts = sorted([*seenProteins], key=len, reverse=True)
        sortedPath = writeTempFasta(sortedSeenProts)
    # Otherwise merge the sorted Temp files so that we have one file containing all seenPeptides in sorted order.
    else:
        sortedPath = mergeSortedFiles(sortedTempFileNames)

    # If there are proteins left over, finish it off, also deal with if there was enough memory initally
    if seenProteins:
        iterTempName = writeTempFasta(seenProteins)
        iterTempFileNames.append(iterTempName)

    # remove subsets if the user has input to do so. Takes flags for keeping origin data and writing subsets
    # to file. Also takes a list of the temp files which contain all seen peptides, and takes the
    # output path to the sorted protein file.
    #if removeSubFlag:
    refinedRemoveSubsetSeq(originFlag, writeSubFlag, sortedPath, iterTempFileNames, outputPath)
    os.remove(sortedPath)

def combineAllTempFasta(outputTempFiles, ignoreNames=False, writeSubsets=False):
    # file Two none just to deal with if only one file in the outputTempFIles queue
    fileTwo = None
    while not outputTempFiles.empty():

        fileOne = outputTempFiles.get()
        # Only one file left over, tying to get fileTwo from empty queue would fail
        if outputTempFiles.empty():
            break

        fileTwo = outputTempFiles.get()

        if outputTempFiles.empty():
            break

        seenPeptides = combineTempFile(fileOne, fileTwo, ignoreNames, writeSubsets)

        tempName = writeTempFasta(seenPeptides)
        outputTempFiles.put(tempName)

    finalSeenPeptides = combineTempFile(fileOne, fileTwo, ignoreNames, writeSubsets)

    # Return the last combination of two files remaining
    return finalSeenPeptides

def combineTempFile(fileOne, fileTwo, ignoreNames, writeSubsets):
    with open(fileOne, 'rU') as handle:

        if ignoreNames:
            # Messy :(
            if writeSubsets:
                seenPeptides = set()
                for line in handle:
                    line = line.strip()
                    seenPeptides.add(line)
            else:
                seenPeptides = set()
                for line in handle:
                    line = line.strip()
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
    # Only open file Two if not none
    if fileTwo:
        with open(fileTwo, 'rU') as handle:
            if ignoreNames:
                if writeSubsets:
                    for line in handle:
                        line = line.strip()
                        seenPeptides.add(line)
                else:
                    for line in handle:
                        line = line.strip()
                        seenPeptides.add(line)

            else:
                for record in SeqIO.parse(handle, 'fasta'):

                    peptide = str(record.seq)
                    protein = str(record.name)
                    if peptide not in seenPeptides.keys():
                        seenPeptides[peptide] = [protein]
                    else:
                        seenPeptides[peptide].append(protein)
        os.remove(fileTwo)

    # Delete temp files as they are used up
    os.remove(fileOne)
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
    # Adapted from stack overflow
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