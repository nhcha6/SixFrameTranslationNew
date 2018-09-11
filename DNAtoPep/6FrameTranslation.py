from DNA_CODON_TABLE import *
from Bio import SeqIO
import time


# hello
def seqToProtein(dnaSeq):
    start = time.time()
    forwFrames, revFrames = seqToFrames(dnaSeq)
    end = time.time()
    #return 0
    print(forwFrames)
    print(revFrames)
    aminoFrames = []
    for frame in forwFrames:
        amino = tripletToAmino(frame)
        aminoFrames.append(amino)
    for frame in revFrames:
        amino = tripletToAmino(frame)
        aminoFrames.append(amino)
    print(end-start)

    return aminoFrames

def seqToFrames(dnaSeq):
    forward = dnaSeq
    reverse = createReverseSeq(dnaSeq)
    print("Forward is: " + forward)
    print("Reverse is: " + reverse)
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
    aminoList = ""
    for triplet in frame:
        amino = DNA_TABLE[triplet]
        if amino == -1:
            aminoList+= " STOP "
        else:
            aminoList += str(amino)

    return aminoList

def parseFastaDna(input_path):
    fasta_sequences = SeqIO.parse(open(input_path), 'fasta')
    sequenceDictionary = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequence = sequence.upper().replace("N", "")
        #sequenceDictionary[name] = sequence.upper()
        print(seqToProtein(sequence))

        break
    return sequenceDictionary

#parseFastaDna('C:/Users/Arpit/Desktop/DNAtoPep/InputData/hg38.fa')

aminoFrames = seqToProtein('ACTGACTGATCTGACTA')
print(aminoFrames)