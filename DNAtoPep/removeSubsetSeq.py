from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from time import time
from queue import Queue
import os
import tempfile
import SixFrameTranslation as sf

# boolean which sets if the origin data in the fasta record name is ignored or not. True means it is ignored.
ignoreNames = False
# input file name and sortedFile name
outputPath = 'RemoveSubset/100,000-Record'
writeSubseqs = True

# create sequence object adapted from the Mers code to account for the input of either a dict or a a set
def createSeqObj(seenPeptides):
    """
    Given the set of matchedPeptides, converts all of them into SeqRecord objects and passes back a generator
    """
    count = 1
    seqRecords = []
    try:
        for sequence, name in seenPeptides.items():
            finalId = "dna|pro" + str(count) + ';'
            count+=1

            # used to convey where the protein was derived from. We may need to do something similar
            for protein in name:
                finalId += protein + ';'
            yield SeqRecord(Seq(str(sequence)), id=finalId, description="")

    except AttributeError:
        for sequence in seenPeptides:
            finalId = "ipd|pep" + str(count) + ';'
            yield SeqRecord(Seq(str(sequence)), id=str(finalId), description="")
            count += 1
    return seqRecords

def seenPepList(filePath):
    seenPeptides = []
    # iterate through the 6FrameTranslation output file and add all sequences to a list
    with open(filePath, "rU") as handle:
        origNo = 0
        for record in SeqIO.parse(handle, 'fasta'):
            origNo += 1
            seenPeptides.append(str(record.seq))
    return seenPeptides, origNo

def sortList(savePath, seenPeptides):
    # sort the list of sequences so that it is ordered from longest to shortest
    seenPeptides.sort(key=len, reverse=True)
    seenPeptides.reverse()
    # write sorted list to fasta file so it need not be stored in memory
    with open(savePath, "w") as output_handle:
        SeqIO.write(createSeqObj(seenPeptides, False), output_handle, "fasta")

def seenPepSetDict(filePath, ignoreOrigin):
    # if we have chose to ignore the sequence names (which stores information about the origin peptides) we need only
    # convert the sorted list to a set before preceding.
    # if ignoreOrigin:
    #     seenPeptides = set(seenPeptides)
    # # however, if we wish to keep the origin data in the final output, we need to make seenPeptides and empty set,
    # # re-iterate through the 6FrameTranslation output file and build up a dictionary which stores the sequences as
    # # a key and the origins as a value.
    # else:
    #     seenPeptides = {}
    #     with open(filePath, "rU") as handle:
    #         for record in SeqIO.parse(handle, 'fasta'):
    #             seenPeptides[str(record.seq)] = record.name
    # return seenPeptides

    with open(filePath, "rU") as handle:
        if ignoreOrigin:
            seenPeptides = set()
        else:
            seenPeptides = {}
        for record in SeqIO.parse(handle, 'fasta'):
            if ignoreOrigin:
                seenPeptides.add(str(record.seq))
            else:
                seenPeptides[str(record.seq)] = record.name
    return seenPeptides


def pepRemoveNoOrigin(handle, seenPeptides, writeSubsets):
    # iterate through each record in the sorted.fasta
    print('ignoring origins')
    counter = 1
    seenSubsets = set()
    for pep in handle:
        # print the counter ever 10000 records to track progress of code without slowing it by printing too much.
        if counter % 10000 == 0:
            print(counter)
        counter += 1
        pep = pep.strip()
        # if the current peptide has already been deleted from seenPeptides because it is a subset of a longer
        # peptide, we simply move on as any subsets of this peptide will have been deleted also.
        # if pep not in seenPeptides:
        #     continue
        # for peptides which are not subsets of any others, we create all the possible subset sequences (same as
        # linear splits) in order to check if they match any other sequence in seenPeptides
        for i in range(len(pep)):
            for j in range(i + 1, len(pep) + 1):
                # if a given split/subsequences is in seenPeptides it is deleted from seenPeptides. Additional check
                # is to ensure that we do not delete the original peptide which has been split up from seenPeptides
                if pep[i:j] in seenPeptides and pep[i:j] != pep:
                    seenPeptides.remove(pep[i:j])
                    # Add the deleted subsequence to seenSubsets if the user has set writeSubsets to True
                    if writeSubsets:
                        seenSubsets.add(pep[i:j])
    # if writeSubsets is True, two sets need to be returned
    if writeSubsets:
        return seenPeptides, seenSubsets
    else:
        return seenPeptides

def pepRemoveOrigin(handle, seenPeptides, writeSubsets):
    # this code is the same as pepRemoveNoOrigin, except it is written expecting seenPeptides to be a dictionary so that
    # the origin/name data can be included.
    counter = 1
    seenSubsets = {}
    for pep in handle:
        # Get rid of new line from tempfile
        pep = pep.strip()
        if counter % 10000 == 0:
            print(counter)
        counter += 1
        # if pep not in seenPeptides.keys():
        #     continue
        for i in range(len(pep)):
            for j in range(i + 1, len(pep)+1):
                currSplit = pep[i:j]

                if currSplit in seenPeptides.keys() and currSplit != pep:
                    if writeSubsets:
                        seenSubsets[pep[i:j]] = seenPeptides[pep[i:j]]

                    del seenPeptides[pep[i:j]]
    # if writeSubsets is True, two sets need to be returned
    if writeSubsets:
        return seenPeptides, seenSubsets
    else:
        return seenPeptides


def refinedRemoveSubsetSeq(ignoreNames, writeSubsets, sortedPath, iterTempFiles, outputPath):

    noSubseqPath = outputPath + "-NoSubsets.fasta"
    onlySubseqPath = outputPath + "-OnlySubsets.fasta"

    subSeqTempFiles = Queue()
    seenPepTempFiles = Queue()
    for currentFile in iterTempFiles:

        with open(sortedPath, 'rU') as handle:

            seenPeptides = seenPepSetDict(currentFile, ignoreNames)
            if ignoreNames:
                # if we want to writeSubsets to a separate file, we will receive to outputs.
                if writeSubsets:
                    seenPeptides, seenSubseqs = pepRemoveNoOrigin(handle, seenPeptides, writeSubsets)
                    print("seen peptides ignore is: " + str(seenPeptides))
                    print("seen subseqs ignore is: " + str(seenSubseqs))
                    if seenSubseqs:
                        subseqTemp = sf.writeTempFasta(seenSubseqs)
                        subSeqTempFiles.put(subseqTemp)
                else:
                    seenPeptides = pepRemoveNoOrigin(handle, seenPeptides, writeSubsets)
                    print("seen peptides: " + str(seenPeptides))

            else:
                if writeSubsets:
                    seenPeptides, seenSubseqs = pepRemoveOrigin(handle, seenPeptides, writeSubsets)
                    print("seen peptides false ignore is: " + str(seenPeptides))
                    print("seen subseqs false ignore is: " + str(seenSubseqs))
                    subseqTemp = sf.writeTempFasta(seenSubseqs)
                    subSeqTempFiles.put(subseqTemp)
                else:
                    seenPeptides = pepRemoveOrigin(handle, seenPeptides, writeSubsets)
            if seenPeptides:
                print("Writing seen peps to temp")
                seenPepTemp = sf.writeTempFasta(seenPeptides)
                seenPepTempFiles.put(seenPepTemp)
        # Remove the initial files as they are used
        os.remove(currentFile)

    finalSeenPeptides = sf.combineAllTempFasta(seenPepTempFiles, ignoreNames)

    # write the new, smaller seenPeptides to file
    with open(noSubseqPath, "w") as output_handle:
        SeqIO.write(createSeqObj(finalSeenPeptides), output_handle, "fasta")
    if writeSubsets:
        print("AABOUT TO COMBINE WRITE SUBSETS")
        subSeqList = list(subSeqTempFiles.queue)
        print(subSeqList)
        finalSubSeqs = sf.combineAllTempFasta(subSeqTempFiles, ignoreNames, writeSubsets)
        print("Final subseqs: " + str(finalSubSeqs))

        with open(onlySubseqPath, "w") as output_handle:
            SeqIO.write(createSeqObj(finalSubSeqs), output_handle, "fasta")


def removeSubsetSeq(ignoreNames, writeSubsets, outputPath):
    # establish inputPath, sortedPath and noSubseq/onlySubseq outputPath
    inputPath = outputPath + "-All.fasta"
    sortedPath = outputPath + "-Sorted.fasta"
    noSubseqPath = outputPath + "-NoSubsets.fasta"
    onlySubseqPath = outputPath + "-OnlySubsets.fasta"

    # multiple timing variables used to time separate parts of the code
    time1 = time()

    # create list of all peptides recorded in the file
    seenPeptides, origNo = seenPepList(inputPath)

    time2 = time()

    print('Number of original entries: ' + str(origNo))
    print('Time to read all entries to list: ' + str(time2-time1))

    # sort list of peptides and write to new file
    sortList(sortedPath, seenPeptides)

    time3 = time()
    print('All sequences sorted and written to new dictionary in: ' + str(time3-time2))

    # create seenPeptides dict or set (depending on if we care about the origin data)
    seenPeptides = seenPepSetDict(seenPeptides, inputPath, ignoreNames)

    time4 = time()
    print('Second read of input file took: ' + str(time4-time3))

    # open sorted fasta to iterate through the sequences in order from largest to smallest. We do so because
    # a large peptide will not be a subset of a smaller one, and thus we can delete peptides sooner and reduce
    # the runtime of the algorithm by starting with the largest.
    with open(sortedPath, "rU") as handle:
        counter = 1
        # only difference between using ignoring names and not ignoring them is that seenPeptides will be a different data
        # structure (set if ignoreNames is True, dict if False) and thus requires different syntax.
        if ignoreNames:
            # if we want to writeSubsets to a separate file, we will receive to outputs.
            if writeSubsets:
                seenPeptides, seenSubseqs = pepRemoveNoOrigin(handle, seenPeptides, writeSubsets)
            else:
                seenPeptides = pepRemoveNoOrigin(handle, seenPeptides, writeSubsets)
        else:
            if writeSubsets:
                seenPeptides, seenSubseqs = pepRemoveOrigin(handle, seenPeptides, writeSubsets)
            else:
                seenPeptides = pepRemoveOrigin(handle, seenPeptides, writeSubsets)

    # remove sorted.fasta from where it is saved
    os.remove(sortedPath)

    print('Time to delete subset sequences: ' + str(time()-time4))
    print('No. of sequences reduced from ' + str(origNo) + ' to ' + str(len(seenPeptides)))

    # write the new, smaller seenPeptides to file
    with open(noSubseqPath, "w") as output_handle:
        SeqIO.write(createSeqObj(seenPeptides, False), output_handle, "fasta")

    # if writeSubsets is True, write seenSubsets to file
    if writeSubsets:
        with open(onlySubseqPath, "w") as output_handle:
            SeqIO.write(createSeqObj(seenSubseqs, False), output_handle, "fasta")

#removeSubsetSeq(ignoreNames, writeSubseqs, outputPath)