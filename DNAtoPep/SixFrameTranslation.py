from DNA_CODON_TABLE import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
import multiprocessing
from multiprocessing import Queue
import math
from removeSubsetSeq import *
import logging
import traceback
import io

numProc = 3

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

def seqToProteinNew(seqDict, minLen, procNum):
    try:
        nameProtTups = []
        for name, dnaSeq in seqDict.items():
            # NEED TO COUNT HOW MANY N'S At start and end, and remove them.
            newSeq = str(dnaSeq).upper()
            newSeq = removeNsDNA(newSeq)

            newSeq = newSeq.upper().replace('N', 'X')

            proteins = buildForwProt(newSeq, minLen) + buildRevProt(newSeq, minLen)
            nameProtTups.append((name, proteins))

        seqToProteinNew.toWriteQueue.put(nameProtTups)
        print("Process number " + str(procNum) + " completed!")

    except Exception as e:

        exc_buffer = io.StringIO()

        traceback.print_exc(file=exc_buffer)

        errorString = 'Uncaught exception in worker process: ' + str(procNum) + '\n%s'

        logging.error(

            errorString,

            exc_buffer.getvalue())

        raise e


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

    start = time.time()
    num_workers = multiprocessing.cpu_count()
    toWriteQueue = multiprocessing.Queue()
    writerProcess = multiprocessing.Process(target=writer, args=(toWriteQueue, outputPath, removeSubFlag, writeSubFlag, originFlag))
    writerProcess.start()
    pool = multiprocessing.Pool(processes=num_workers, initializer=poolInitialiser,
                                initargs=(toWriteQueue,))

    # calculate total size of input fasta
    with open(input_path, "rU") as handle:
        totalProt = 0
        for entry in SeqIO.parse(handle, 'fasta'):
            totalProt += 1


    pepPerProc = math.ceil(totalProt/numProc)
    print("Process Size: " + str(pepPerProc))

    with open(input_path, "rU") as handle:
        counter = 0
        seqDict = {}
        procNum = 0
        for record in SeqIO.parse(handle, 'fasta'):
            name = "rec" + str(counter) + ';'
            dnaSeq = record.seq

            seqDict[name] = dnaSeq
            counter += 1
            if counter % pepPerProc == 0:
                procNum += 1
                # create process
                print("Starting process number: " + str(procNum))
                pool.apply_async(seqToProteinNew, args=(seqDict, minLen, procNum))
                seqDict = {}

        procNum += 1
        # create process
        print("Starting process number: " + str(procNum))
        pool.apply_async(seqToProteinNew, args=(seqDict, minLen, procNum))
        seqDict = {}

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
             finalId+=protein

        yield SeqRecord(Seq(key), id=finalId, description="")

        count += 1

    return seqRecords


def writer(queue, outputPath, removeSubFlag, writeSubFlag, originFlag):
    seenProteins = {}
    saveHandle = outputPath + '-All.fasta'
    with open(saveHandle, "w") as output_handle:
        counter = 0
        while True:
            counter += 1
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
                        seenProteins[protein].append(name)
            print("Got from Queue: " + str(counter))
        print("writing to fasta")
        SeqIO.write(createSeqObj(seenProteins), output_handle, "fasta")

    # if removal of subsets is selected, run code.
    if removeSubFlag:
        removeSubsetSeq(originFlag, writeSubFlag, outputPath)

def poolInitialiser(toWriteQueue):
    seqToProteinNew.toWriteQueue = toWriteQueue


def createReverseSeq(dnaSeq):
    reverseDir = dnaSeq[::-1]
    _tab = str.maketrans(dict(zip('ATCG', 'TAGC')))
    reverseSeq = reverseDir.translate(_tab)
    return reverseSeq

