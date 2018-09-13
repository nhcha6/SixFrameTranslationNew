from DNA_CODON_TABLE import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
import multiprocessing
from multiprocessing import Queue
import logging

#MIN_PEPTIDE_LEN = 7

# set of new function which don't require the storage of all forward and reverse frames to run
def buildForwProt(seq, minLen):
    proteins = []
    for j in range(0,3):
        proteinTemp = ""
        remainder = (len(seq[j:-1]) + 1) % 3
        for i in range(j, len(seq) - remainder, 3):
            codon = seq[i:i+3]
            amino = DNA_TABLE[codon]
            if amino == -1:
                if len(proteinTemp) > minLen:
                    proteins.append(proteinTemp)
                    proteinTemp = ""
            elif i == len(seq) - remainder - 3:
                proteinTemp += amino
                if len(proteinTemp) > minLen:
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
                codon  = createReverseSeq(seq[-3:])
            else:
                codon = createReverseSeq(seq[-1*(i+3):-i])
            amino = DNA_TABLE[codon]
            if amino == -1:
                if len(proteinTemp) > minLen:
                    proteins.append(proteinTemp)
                    proteinTemp = ""
            elif i == len(seq) - remainder - 3:
                proteinTemp += amino
                if len(proteinTemp) > minLen:
                    proteins.append(proteinTemp)
                proteinTemp = ""
            else:
                proteinTemp += amino
    return proteins

def seqToProteinNew(dnaSeq, minLen, name):
    newSeq = dnaSeq.upper().replace('N', '')
    start = time.time()
    proteins = buildForwProt(newSeq, minLen) + buildRevProt(newSeq, minLen)
    seqToProteinNew.toWriteQueue.put([name, proteins])
    end = time.time()
    #print(end - start)
    return proteins

def generateOutputNew(outputPath, minLen, input_path):
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
            seq = record.seq
            dnaSeq = str(seq).upper()
            pool.apply_async(seqToProteinNew, args=(dnaSeq, minLen, name))
            #proteins = seqToProteinNew(dnaSeq, minLen)
            #toWriteQueue.put([name, proteins])
            if counter%10000==0:
                print(counter)
                if counter == 100000:
                    break
    pool.close()
    pool.join()
    toWriteQueue.put('stop')

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

# old functions (slightly more time efficient) which store forward and reverse frames in memory and don't use
# any multiprocessing.
def seqToProtein(dnaSeq, minLen):

    newSeq = dnaSeq.upper().replace('N', '')
    start = time.time()
    forwFrames, revFrames = seqToFrames(newSeq)
    peptides = []

    for frame in forwFrames:
        peptide = tripletToAmino(frame, minLen)
        peptides += peptide

    for frame in revFrames:
        peptide = tripletToAmino(frame, minLen)
        peptides += peptide

    end = time.time()

    print(end-start)

    return peptides


def seqToFrames(dnaSeq):
    forward = dnaSeq
    reverse = createReverseSeq(dnaSeq)

    forwardFrames = createFrames(forward)
    reverseFrames = createFrames(reverse)
    return forwardFrames, reverseFrames


def createFrames(dnaSeq):
    frames = [[],[],[]]
    for i in range(0,3):
        frame = frames[i]
        for j in range(0, len(dnaSeq), 3):
            if i+j < len(dnaSeq)-2:
                triplet = dnaSeq[i+j:i+j+3]
                frame.append(triplet)

    return frames


def createReverseSeq(dnaSeq):
    reverseDir = dnaSeq[::-1]
    _tab = str.maketrans(dict(zip('ATCG', 'TAGC')))
    reverseSeq = reverseDir.translate(_tab)
    return reverseSeq

# incorporate start triplet
def tripletToAmino(frame, minLen):
    aminoList = []
    peptideList = []

    for triplet in frame:
        amino = DNA_TABLE[triplet]
        if amino == -1:
            if len(aminoList) > minLen:
                peptideList.append(''.join(aminoList))
                aminoList.clear()
        else:
            aminoList.append(amino)
    if len(aminoList) > minLen:
        peptideList.append(''.join(aminoList))
        aminoList.clear()
    return peptideList

def parseFastaDna(input_path):
    # fasta_sequences = SeqIO.parse(open(input_path), 'fasta')
    sequenceDictionary = {}
    # for fasta in fasta_sequences:
    #     name, sequence = fasta.id, str(fasta.seq)
    #     sequence = sequence.upper().replace("N", "")
    #     #sequenceDictionary[name] = sequence.upper()
    #     print(seqToProtein(sequence))
    #     break;

    with open(input_path, "rU") as handle:
        counter = 0

        # convert to tuple and look to start multiprocessing from here
        for record in SeqIO.parse(handle, 'fasta'):
            counter += 1
            print(counter)
            sequenceDictionary[record.name] = record.seq

    return sequenceDictionary

def generateOutput(outputPath, minLen, inputFile):
    finalPeptides = {}
    seqDict = parseFastaDna(inputFile)
    for key, value in seqDict.items():
        dnaSeq = str(value).upper()
        peptides = seqToProteinNew(dnaSeq, minLen)
        for peptide in peptides:
            if peptide not in finalPeptides.keys():
                finalPeptides[peptide] = [key]
            else:
                finalPeptides[peptide].append(key)
    print(finalPeptides)
    saveHandle = outputPath + '/DNAFastaProteins.fasta'
    with open(saveHandle, "w") as output_handle:
        SeqIO.write(createSeqObj(finalPeptides), output_handle, "fasta")

# parseFastaDna('C:/Users/Arpit/Desktop/DNAtoPep/InputData/hg38.fa')

#seqDict = parseFastaDna('/Users/nicolaschapman/Documents/UROP/6FrameTranslation/DNAtoPep/DNAsmall.fasta')
#generateProteins('/Users/nicolaschapman/Documents/UROP/6FrameTranslation/DNAtoPep/DNAsmall.fasta')
#seqRecords = createSeqObj(finPep)
#print(seqDict)
#
# aminoFrames = seqToProtein('NNNNNNNNNNNNNNNNNNNNNNNNNACTGACTGATCTGACTANNNNNNNN')
# # print(aminoFrames)
# print(aminoFrames)