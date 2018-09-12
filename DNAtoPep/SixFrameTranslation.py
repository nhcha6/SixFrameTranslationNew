from DNA_CODON_TABLE import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time

MIN_PEPTIDE_LEN = 7

# hello
def seqToProtein(dnaSeq):

    newSeq = dnaSeq.upper().replace('N', '')
    start = time.time()
    forwFrames, revFrames = seqToFrames(newSeq)
    peptides = []

    for frame in forwFrames:
        peptide = tripletToAmino(frame)
        peptides += peptide

    for frame in revFrames:
        peptide = tripletToAmino(frame)
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
def tripletToAmino(frame):
    aminoList = []
    peptideList = []

    for triplet in frame:
        amino = DNA_TABLE[triplet]
        if amino == -1:
            if len(aminoList) > MIN_PEPTIDE_LEN:
                peptideList.append(''.join(aminoList))
                aminoList.clear()
        else:
            aminoList.append(amino)
    if len(aminoList) > MIN_PEPTIDE_LEN:
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
        for record in SeqIO.parse(handle, 'fasta'):
            sequenceDictionary[record.name] = record.seq

    return sequenceDictionary

def generateProteins(input_path):
    seqDict = parseFastaDna(input_path)
    finalPeptides = {}
    for key, value in seqDict.items():
        dnaSeq = str(value).upper()
        peptides = seqToProtein(dnaSeq)
        for peptide in peptides:
            if peptide not in finalPeptides.keys():
                finalPeptides[peptide] = [key]
            else:
                finalPeptides[peptide].append(key)
    print(finalPeptides)
    saveHandle = '/Users/nicolaschapman/Documents/UROP/6FrameTranslation/DNAtoPep/dnaFasta.fasta'
    with open(saveHandle, "w") as output_handle:
        SeqIO.write(createSeqObj(finalPeptides), output_handle, "fasta")


def createSeqObj(finalPeptides):
    """
    Given the set of matchedPeptides, converts all of them into SeqRecord objects and passes back a generator
    """
    count = 1
    seqRecords = []

    for key, value in finalPeptides.items():

        finalId = "dna|pep"+str(count)+';'

        # used to convey where the protein was derived from. We may need to do something similar
        for protein in value:
             finalId+=protein+';'

        yield SeqRecord(Seq(key), id=finalId, description="")

        count += 1

    return seqRecords




# parseFastaDna('C:/Users/Arpit/Desktop/DNAtoPep/InputData/hg38.fa')

seqDict = parseFastaDna('/Users/nicolaschapman/Documents/UROP/6FrameTranslation/DNAtoPep/DNAsmall.fasta')
generateProteins('/Users/nicolaschapman/Documents/UROP/6FrameTranslation/DNAtoPep/DNAsmall.fasta')
#seqRecords = createSeqObj(finPep)
#print(seqDict)
#
# aminoFrames = seqToProtein('NNNNNNNNNNNNNNNNNNNNNNNNNACTGACTGATCTGACTANNNNNNNN')
# # print(aminoFrames)
# print(aminoFrames)