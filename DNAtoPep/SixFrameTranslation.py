from DNA_CODON_TABLE import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
import multiprocessing
from multiprocessing import Queue
import logging


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

def generateOutputNew(outputPath, minLen, input_path):

    start = time.time()
    num_workers = multiprocessing.cpu_count()
    toWriteQueue = multiprocessing.Queue()
    writerProcess = multiprocessing.Process(target=writer, args=(toWriteQueue, outputPath))
    writerProcess.start()
    pool = multiprocessing.Pool(processes=num_workers, initializer=poolInitialiser,
                                initargs=(toWriteQueue,))
    with open(input_path, "rU") as handle:
        counter = 0
        for record in SeqIO.parse(handle, 'fasta'):
            counter += 1
            name = record.name
            dnaSeq = record.seq

            print("Starting process for " + str(dnaSeq))

            pool.apply_async(seqToProteinNew, args=(dnaSeq, minLen, name))
            #proteis = seqToProteinNew(dnaSeq, minLen)
            #toWriteQueue.put([name, proteins])
            break

    pool.close()
    pool.join()
    toWriteQueue.put('stop')
    writerProcess.join()
    end = time.time()
    print("Altogether took " + str(end-start))

def createSeqObj(finalPeptides):
    """
    Given the set of matchedPeptides, converts all of them into SeqRecord objects and passes back a generator
    """
    count = 1
    seqRecords = []

    for key, value in finalPeptides.items():

        finalId = "dna|pro"+str(count)+';'

        # used to convey where the protein was derived from. We may need to do something similar
        for protein in value:
             finalId+=protein+';'

        yield SeqRecord(Seq(key), id=finalId, description="")

        count += 1

    return seqRecords


def writer(queue, outputPath):
    seenProteins = {}
    saveHandle = outputPath + '/DNAFastaProteins.fasta'
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
                    seenProteins[protein].append(name)

        print("writing to fasta")
        SeqIO.write(createSeqObj(seenProteins), output_handle, "fasta")

def poolInitialiser(toWriteQueue):
    seqToProteinNew.toWriteQueue = toWriteQueue


def createReverseSeq(dnaSeq):
    reverseDir = dnaSeq[::-1]
    _tab = str.maketrans(dict(zip('ATCG', 'TAGC')))
    reverseSeq = reverseDir.translate(_tab)
    return reverseSeq

