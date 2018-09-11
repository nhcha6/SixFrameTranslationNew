from DNA_CODON_TABLE import *
from Bio import SeqIO
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

        peptides.append(peptide)
    for frame in revFrames:
        peptide = tripletToAmino(frame)
        peptides.append(peptide)
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
            seqToProtein(str(record.seq))

    return sequenceDictionary

# parseFastaDna('C:/Users/Arpit/Desktop/DNAtoPep/InputData/hg38.fa')

#parseFastaDna('C:/Users/Arpit/Desktop/DNAtoPep/InputData/smallDNA.fa')

#
# aminoFrames = seqToProtein('NNNNNNNNNNNNNNNNNNNNNNNNNACTGACTGATCTGACTANNNNNNNN')
# print(aminoFrames)