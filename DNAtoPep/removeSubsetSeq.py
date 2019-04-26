from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import time
import os

# boolean which sets if the origin data in the fasta record name is ignored or not. True means it is ignored.
IGNORE_NAMES = False
# input file name and sortedFile name
OUTPUT_PATH = 'RemoveSubset/100,000-Record'
WRITE_SUBSEQ = True

def removeSubsetSeq(ignoreNames, writeSubsets, outputPath):
    """
    This function manages the flow of the script, calling the functions required to remove from the input file
    all proteins which exist as subsets in longer proteins.

    :param outputPath: the input fasta file from which subset proteins are to be deleted.
    :return:
    """
    # establish inputPath, sortedPath and noSubseq/onlySubseq outputPath
    inputPath = outputPath + "-All.fasta"
    sortedPath = outputPath + "-Sorted.fasta"
    noSubseqPath = outputPath + "-NoSubsets.fasta"
    onlySubseqPath = outputPath + "-OnlySubsets.fasta"

    # multiple timing variables used to time separate parts of the code
    time1 = time.time()

    # create list of all proteins recorded in the file
    seenProteins, origNo = seenPepList(inputPath)

    time2 = time.time()

    print('Number of original entries: ' + str(origNo))
    print('Time to read all entries to list: ' + str(time2-time1))

    # sort list of proteins and write to new file
    sortList(sortedPath, seenProteins)

    time3 = time.time()
    print('All sequences sorted and written to new dictionary in: ' + str(time3-time2))

    # create seenProteins dict or set (depending on if we care about the origin data)
    seenProteins = seenPepSetDict(seenProteins, inputPath, ignoreNames)

    time4 = time.time()
    print('Second read of input file took: ' + str(time4-time3))

    # open sorted fasta to iterate through the sequences in order from largest to smallest. We do so because
    # a large protein will not be a subset of a smaller one, and thus we can delete proteins sooner and reduce
    # the runtime of the algorithm by starting with the largest.
    with open(sortedPath, "rU") as handle:
        counter = 1
        # only difference between using ignoring names and not ignoring them is that seenProteins will be a different data
        # structure (set if ignoreNames is True, dict if False) and thus requires different syntax.
        if ignoreNames:
            # if we want to writeSubsets to a separate file, we will receive to outputs.
            if writeSubsets:
                seenProteins, seenSubseqs = pepRemoveNoOrigin(handle, seenProteins, writeSubsets)
            else:
                seenProteins = pepRemoveNoOrigin(handle, seenProteins, writeSubsets)
        else:
            if writeSubsets:
                seenProteins, seenSubseqs = pepRemoveOrigin(handle, seenProteins, writeSubsets)
            else:
                seenProteins = pepRemoveOrigin(handle, seenProteins, writeSubsets)

    # remove sorted.fasta from where it is saved
    os.remove(sortedPath)

    print('Time to delete subset sequences: ' + str(time.time()-time4))
    print('No. of sequences reduced from ' + str(origNo) + ' to ' + str(len(seenProteins)))

    # write the new, smaller seenProteins to file
    with open(noSubseqPath, "w") as output_handle:
        SeqIO.write(createSeqObj(seenProteins), output_handle, "fasta")

    # if writeSubsets is True, write seenSubsets to file
    if writeSubsets:
        with open(onlySubseqPath, "w") as output_handle:
            SeqIO.write(createSeqObj(seenSubseqs), output_handle, "fasta")

def sortList(savePath, seenProteins):
    """
    Called by removeSubsetSeq(), this function takes an output file location and a list of proteins, and write the
    proteins to file in order of longest to shortest.

    :param savePath: the location that the sorted proteins are to be written to.
    :param seenProteins: the unsorted list of proteins to be written to sorted file.
    :return:
    """
    # sort the list of sequences so that it is ordered from longest to shortest
    seenProteins.sort(key=len)
    seenProteins.reverse()
    # write sorted list to fasta file so it need not be stored in memory
    with open(savePath, "w") as output_handle:
        SeqIO.write(createSeqObj(seenProteins), output_handle, "fasta")

def seenPepList(filePath):
    """
    Called by Called by removeSubsetSeq(). From the file path input, this function adds all proteins to a list and
    returns the list along with its size.

    :param filePath: the path of the input protein file.

    :return seenProteins: a list of all the protein sequences in the input file.
    :return origNo: the number of proteins in the list.
    """
    seenProteins = []
    # iterate through the 6FrameTranslation output file and add all sequences to a list
    with open(filePath, "rU") as handle:
        origNo = 0
        for record in SeqIO.parse(handle, 'fasta'):
            origNo += 1
            seenProteins.append(str(record.seq))
    return seenProteins, origNo

def seenPepSetDict(seenProteins, filePath, ignoreOrigin):
    """
    Called by removeSubsetSeq(), this function returns either a set containing all the protein sequences in the
    input file, or a dictionary storing the protein sequences as keys and the sequence name as values. It will store
    them as a set if ingoreOrigins is True, and a dictionary if it is False.

    :param seenProteins: a list of the proteins contained within seenProteins.
    :param filePath: the path to the input protein file. Needed if ignoreOrigins is False to extract the name of
    each protein sequence.
    :param ignoreOrigin: True if the user does not require the name of each protein in the output file which has
    subsequences removed.

    :return seenProteins: the set or dictionary containing the input protein sequence data.
    """
    # if we have chose to ignore the sequence names (which stores information about the origin proteins) we need only
    # convert the sorted list to a set before preceding.
    if ignoreOrigin:
        seenProteins = set(seenProteins)
    # however, if we wish to keep the origin data in the final output, we need to make seenProteins and empty set,
    # re-iterate through the 6FrameTranslation output file and build up a dictionary which stores the sequences as
    # a key and the origins as a value.
    else:
        seenProteins = {}
        with open(filePath, "rU") as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                seenProteins[str(record.seq)] = record.name
    return seenProteins

def pepRemoveNoOrigin(handle, seenProteins, writeSubsets):
    """
    Called by removeSubsetSeq() if the user has chosen to ignore protein names, this function iterates through each
    protein in the sorted proteinList (from longest to shortest) and deletes the proteins from seenProteins
    which are subsets of the current protein.

    :param handle: the filePath of the fasta file containing the sorted proteins.
    :param seenProteins: the set of all proteins from the initial input file.
    :param writeSubsets: True if the user wishes to output a file containing all the subsets which were removed from
    the input file.

    :return seenProteins: the edited set of proteins which has all subsets removed.
    """
    # iterate through each record in the sorted.fasta
    counter = 1
    seenSubsets = set()
    for record in SeqIO.parse(handle, 'fasta'):
        # print the counter ever 10000 records to track progress of code without slowing it by printing too much.
        if counter % 10000 == 0:
            print(counter)
        counter += 1
        pep = str(record.seq)
        # if the current protein has already been deleted from seenProteins because it is a subset of a longer
        # protein, we simply move on as any subsets of this protein will have been deleted also.
        if pep not in seenProteins:
            continue
        # for proteins which are not subsets of any others, we create all the possible subset sequences (same as
        # linear splits) in order to check if they match any other sequence in seenProteins
        for i in range(len(pep)):
            for j in range(i + 1, len(pep) + 1):
                # if a given split/subsequences is in seenProteins it is deleted from seenProteins. Additional check
                # is to ensure that we do not delete the original protein which has been split up from seenProteins
                if pep[i:j] in seenProteins and pep[i:j] is not pep:
                    seenProteins.remove(pep[i:j])
                    # Add the deleted subsequence to seenSubsets if the user has set writeSubsets to True
                    if writeSubsets:
                        seenSubsets.add(pep[i:j])
    # if writeSubsets is True, two sets need to be returned
    if writeSubsets:
        return seenProteins, seenSubsets
    else:
        return seenProteins

def pepRemoveOrigin(handle, seenProteins, writeSubsets):
    """
    Called by removeSubsetSeq() if the user has chosen not to ignore protein names, this function iterates through each
    protein in the sorted proteinList (from longest to shortest) and deletes the proteins from seenProteins
    which are subsets of the current protein.

    :param handle: the filePath of the fasta file containing the sorted proteins.
    :param seenProteins: the set of all proteins from the initial input file.
    :param writeSubsets: True if the user wishes to output a file containing all the subsets which were removed from
    the input file.

    :return seenProteins: the edited set of proteins which has all subsets removed.
    """
    # this code is the same as pepRemoveNoOrigin, except it is written expecting seenProteins to be a dictionary so that
    # the origin/name data can be included.
    counter = 1
    seenSubsets = {}
    for record in SeqIO.parse(handle, 'fasta'):
        if counter % 10000 == 0:
            print(counter)
        counter += 1
        pep = str(record.seq)
        if pep not in seenProteins.keys():
            continue
        for i in range(len(pep)):
            for j in range(i + 1, len(pep) + 1):
                if pep[i:j] in seenProteins.keys() and pep[i:j] is not pep:
                    if writeSubsets:
                        seenSubsets[pep[i:j]] = seenProteins[pep[i:j]]
                    del seenProteins[pep[i:j]]
    # if writeSubsets is True, two sets need to be returned
    if writeSubsets:
        return seenProteins, seenSubsets
    else:
        return seenProteins

def createSeqObj(seenProteins):
    """
    Called by the writer() function, this function takes either a dictionary or a set of final proteins, converts all
    of them into SeqRecord objects and adds them to an iterable via the yield statement. This iterable is passed back
    to the SeqIO.write() method which writes the created SeqRecord to fasta file.

    :param seenProteins: either a dictionary with protein sequences as keys and a list of origin DNA sequences as
    values or a set of protein sequences. Which one is passed through depends on if the user has selected to ignore
    origin data or not.
    :return:
    """
    count = 1
    seqRecords = []
    try:
        for sequence, name in seenProteins.items():
            yield SeqRecord(Seq(sequence), id=name, description="")

    except AttributeError:
        for sequence in seenProteins:
            finalId = "ipd|pep" + str(count) + ';'
            yield SeqRecord(Seq(sequence), id=finalId, description="")
            count += 1

#removeSubsetSeq(IGNORE_NAMES, WRITE_SUBSEQ, OUTPUT_PATH)