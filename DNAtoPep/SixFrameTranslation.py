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
import math


MEMORY_THRESHOLD = 80
totalProcesses = 1000
numProcAlive = 40
logging.basicConfig(level=logging.DEBUG, format='%(message)s')

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

def seqToProteinNew(seqDict, minLen):
    start = time.time()
    nameProtTups = []
    for name, dnaSeq in seqDict.items():
        # NEED TO COUNT HOW MANY N'S At start and end, and remove them.
        newSeq = str(dnaSeq).upper()
        newSeq = removeNsDNA(newSeq)

        newSeq = newSeq.upper().replace('N', 'X')

        proteins = buildForwProt(newSeq, minLen) + buildRevProt(newSeq, minLen)
        nameProtTups.append((name, proteins))

    seqToProteinNew.toWriteQueue.put(nameProtTups)
    seqToProteinNew.protCompletedQueue.put(1)
    end = time.time()
    #print("process Completed")

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

    protCompletedQueue = multiprocessing.Queue()

    pool = multiprocessing.Pool(processes=num_workers, initializer=poolInitialiser,
                                    initargs=(toWriteQueue,protCompletedQueue))

    # count total size of input fasta
    totalSize = 0
    with open(input_path, "rU") as handle:
        for entry in SeqIO.parse(handle, 'fasta'):
            totalSize += 1

    # calculate number of proteins per process
    procSize = math.ceil(totalSize/totalProcesses)
    print("Process Size: " + str(procSize))

    with open(input_path, "rU") as handle:
        procGenerated = 0
        counter = 0
        completedProts = 0
        seqDict = {}
        for record in SeqIO.parse(handle, 'fasta'):
            name = "rec" + str(counter) + ';'
            dnaSeq = record.seq

            seqDict[name] = dnaSeq
            counter += 1

            if counter % procSize == 0:
                # slow process generation
                print(procGenerated)
                if procGenerated > numProcAlive:
                    while True:
                        if not protCompletedQueue.empty():
                            completedProts += protCompletedQueue.get()
                            break
                print(completedProts)

                # create process
                pool.apply_async(seqToProteinNew, args=(seqDict, minLen))
                procGenerated += 1
                seqDict = {}

    pool.close()
    pool.join()

    toWriteQueue.put('stop')
    writerProcess.join()
    end = time.time()
    print("Altogether took " + str(end-start))



def writer(queue, outputPath, writeSubFlag, originFlag):
    seenProteins = {}
    saveHandle = outputPath + '-All.fasta'
    with open(saveHandle, "w") as output_handle:
        while True:
            tuples = queue.get()
            if tuples == 'stop':
                print("All proteins added to writer queue")
                break

            for tuple in tuples:
                proteins = tuple[1]
                name = tuple[0]
                for protein in proteins:
                    if protein not in seenProteins.keys():
                        seenProteins[protein] = [name]
                    else:
                        if name not in seenProteins[protein]:
                            seenProteins[protein].append(name)
        print("writing to fasta")
        SeqIO.write(createSeqObjMain(seenProteins), output_handle, "fasta")

    print('removing subset sequences')
    removeSubsetSeq(originFlag, writeSubFlag, outputPath)

# create sequence object adapted from the Mers code to account for the input of either a dict or a a set
def createSeqObjMain(seenPeptides):
    """
    Given the set of matchedPeptides, converts all of them into SeqRecord objects and passes back a generator
    """
    count = 1
    for sequence, names in seenPeptides.items():
        finalId = "ipd|pep" + str(count) + ';'
        for name in names:
            finalId += name
        yield SeqRecord(Seq(sequence), id=finalId, description="")
        count += 1

def poolInitialiser(toWriteQueue, protCompletedQueue):
    seqToProteinNew.toWriteQueue = toWriteQueue
    seqToProteinNew.protCompletedQueue = protCompletedQueue


def createReverseSeq(dnaSeq):
    reverseDir = dnaSeq[::-1]
    _tab = str.maketrans(dict(zip('ATCG', 'TAGC')))
    reverseSeq = reverseDir.translate(_tab)
    return reverseSeq

