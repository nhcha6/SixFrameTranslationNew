from DNA_CODON_TABLE import *

# hello
def seqToProtein(dnaSeq):
    forwFrames, revFrames = seqToFrames(dnaSeq)
    print(forwFrames)
    print(revFrames)
    aminoFrames = []
    for frame in forwFrames:
        amino = tripletToAmino(frame)
        aminoFrames.append(amino)
    for frame in revFrames:
        amino = tripletToAmino(frame)
        aminoFrames.append(amino)
    return aminoFrames

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
    reverseSeq = ""
    for i in range(0,len(reverseDir)):
        reverseSeq += nucleotidePair[reverseDir[i]]
    return reverseSeq

def tripletToAmino(frame):
    aminoList = []
    for triplet in frame:
        amino = DNA_TABLE[triplet]
        aminoList.append(amino)
    return aminoList

aminoFrames = seqToProtein('ACTGACTG')
print(aminoFrames)