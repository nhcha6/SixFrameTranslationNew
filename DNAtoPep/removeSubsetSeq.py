from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
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
                finalId += protein
            yield SeqRecord(Seq(str(sequence)), id=finalId, description="")

    except AttributeError:
        for sequence in seenPeptides:
            finalId = "ipd|pep" + str(count) + ';'
            yield SeqRecord(Seq(str(sequence)), id=str(finalId), description="")
            count += 1
    return seqRecords

def seenPepSetDict(filePath, ignoreOrigin):
    # if we have chose to ignore the sequence names (which stores information about the origin peptides) we need only
    # convert the sorted list to a set before preceding.

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

    # create a queue to store subseq temp files and seenPeptide temp files.
    subSeqTempFiles = Queue()
    seenPepTempFiles = Queue()
    for currentFile in iterTempFiles:

        with open(sortedPath, 'rU') as handle:
            # read the current temp file into seenPeptides as either a set or dict depending on
            # ignoreNames (also called originsFlag).
            seenPeptides = seenPepSetDict(currentFile, ignoreNames)

            # code for ignoreNames = True and ignoreNames = False is very similar, except different
            # syntax is needed to to deal with a set/dict.
            if ignoreNames:
                # if we want to writeSubsets to a separate file, we will receive two outputs.
                if writeSubsets:
                    # takes the seenPeptides in the current temp file, and removes the ones which are subsets
                    # by iterating through the sortedlist file
                    seenPeptides, seenSubseqs = pepRemoveNoOrigin(handle, seenPeptides, writeSubsets)
                    # write the seenSubseqs to a temporary file if any exist.
                    if seenSubseqs:
                        subseqTemp = sf.writeTempFasta(seenSubseqs)
                        subSeqTempFiles.put(subseqTemp)
                # does the same as above except pepRemoveNoOrigin() only has the one output as we are
                # not concerned with seenSubseqs
                else:
                    seenPeptides = pepRemoveNoOrigin(handle, seenPeptides, writeSubsets)
            # same as above but written for a dictionary.
            else:
                if writeSubsets:
                    seenPeptides, seenSubseqs = pepRemoveOrigin(handle, seenPeptides, writeSubsets)

                    subseqTemp = sf.writeTempFasta(seenSubseqs)
                    subSeqTempFiles.put(subseqTemp)
                else:
                    seenPeptides = pepRemoveOrigin(handle, seenPeptides, writeSubsets)

            # write seenPeptides after subseq deletion to a temp file, and add temp file to the queue.
            if seenPeptides:
                seenPepTemp = sf.writeTempFasta(seenPeptides)
                seenPepTempFiles.put(seenPepTemp)
        # Remove the initial files as they are used
        os.remove(currentFile)

    # read all the seenPeptides from the created temp files into finalSeenPeptides.
    finalSeenPeptides = sf.combineAllTempFasta(seenPepTempFiles, ignoreNames)

    # write the new, smaller seenPeptides to file
    with open(noSubseqPath, "w") as output_handle:
        SeqIO.write(createSeqObj(finalSeenPeptides), output_handle, "fasta")

    # if writeSubsets is True, we also need to combine the subseq temp files and write them to file.
    if writeSubsets:
        finalSubSeqs = sf.combineAllTempFasta(subSeqTempFiles, ignoreNames, writeSubsets)
        with open(onlySubseqPath, "w") as output_handle:
            SeqIO.write(createSeqObj(finalSubSeqs), output_handle, "fasta")


