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
import psutil

NUM_PROC = 1000
NUM_PROC_ALIVE = 40
MEMORY_THRESHOLD = 30

def generateOutputNew(outputPath, minLen, input_path, removeSubFlag, writeSubFlag, originFlag):
    """
    Called from Example.createOutput() in SixFrameGUI.py to split the six frame translation of an entire DNA fasta
    file into smaller processes which each deal with a small number of DNA sequences.

    :param outputPath: the full path of the output file.
    :param minLen: tbe minimum length an output protein must be to be included in the output.
    :param input_path: the path of the input fasta file containing DNA sequences.
    :param removeSubFlag: True if the user wants a second output file, which has sequences which are a subset of another
    longer sequence removed.
    :param writeSubFlag: True if removeSubFlag is True and the user wishes to create a fasta file containing all
    of the sequences which were subsets of others and were subsequently deleted from the second output.
    :param originFlag: True if removeSubFlag is True and the user wishes to ignore the origin data (data on which
    DNA sequence a protein was derived from) when creating the second output fasta file with subsets removed.
    :return:
    """
    # initialise some of the queues and processes required for multiprocessing.
    start = time.time()
    num_workers = multiprocessing.cpu_count()
    toWriteQueue = multiprocessing.Queue()
    writerProcess = multiprocessing.Process(target=writer, args=(toWriteQueue, outputPath, removeSubFlag, writeSubFlag, originFlag))
    writerProcess.start()

    # initialise the pool.
    protCompletedQueue = multiprocessing.Queue()
    pool = multiprocessing.Pool(processes=num_workers, initializer=poolInitialiser,
                            initargs=(toWriteQueue,protCompletedQueue))

    # calculate total size of input fasta
    with open(input_path, "rU") as handle:
        totalProt = 0
        for entry in SeqIO.parse(handle, 'fasta'):
            totalProt += 1

    # calculate the number of DNA sequences which will be translated in each process.
    pepPerProc = math.ceil(totalProt/NUM_PROC)
    print("Process Size: " + str(pepPerProc))

    # open the input file.
    with open(input_path, "rU") as handle:
        # initialise required counters and the dictionary to store the DNA sequences for each process.
        counter = 0
        seqDict = {}
        procNum = 0
        completedProts = 0

        # iterate through each record in the input file.
        for record in SeqIO.parse(handle, 'fasta'):
            name = "rec" + str(counter) + ';'
            dnaSeq = record.seq
            # add the DNA sequence name and sequence to the dictionary.
            seqDict[name] = dnaSeq
            counter += 1
            # if the length of seqDict equals the calculated pepPerProc value, we want to start another process
            # and reset seqDict to be an empty dictionary, ready to store the sequences for the next process.
            if counter % pepPerProc == 0:
                procNum += 1
                # once the numProcAlive value has been exceeded, only create process once an
                # alive process has been finished.
                if procNum > NUM_PROC_ALIVE:
                    while True:
                        if not protCompletedQueue.empty():
                            completedProts += protCompletedQueue.get()
                            break
                # create process once while loop is broken, and reset seqDict.
                print("Starting process number: " + str(procNum))
                pool.apply_async(seqToProteinNew, args=(seqDict, minLen, procNum))
                seqDict = {}
                # Check the memory usage. If it exceeds a certain level close the pool as this will clear
                # a lot of old memory.
                if memory_usage_psutil() > MEMORY_THRESHOLD:
                    print('Memory usage exceded. Waiting for processes to finish.')
                    pool.close()
                    pool.join()
                    pool = multiprocessing.Pool(processes=num_workers, initializer=poolInitialiser,
                                                initargs=(toWriteQueue, protCompletedQueue))

        # once the end of the input file is reached, we need to start a final process if seqDict contains the last
        # few processes to be computed.
        if seqDict:
            procNum += 1
            # create process
            print("Starting process number: " + str(procNum))
            pool.apply_async(seqToProteinNew, args=(seqDict, minLen, procNum))
            seqDict = {}

    # close the pool and end the function.
    pool.close()
    pool.join()
    toWriteQueue.put('stop')
    writerProcess.join()
    end = time.time()
    print("Altogether took " + str(end-start))

def memory_usage_psutil():
    """
    Called by generateOutputNew() after a process to check if the program is using close the maximum RAM of
    the computer it is being run on.
    :return mem.percent: the percentage of the computer's RAM being used.
    """
    # return the memory usage in percentage like top
    mem = psutil.virtual_memory()
    return mem.percent

def seqToProteinNew(seqDict, minLen, procNum):
    """
    Called from generateOutputNew() as the worker function to each process in the multiprocessing.Pool(). This
    function takes a dictionary containing DNA seqeunces as values and the name of these sequences as keys. It applies
    the forward three frame and reverse three frame translation (amounting to a six frame translation) to each
    DNA sequence to create a list of proteins that can be created from reading the DNA sequence. It then adds this
    list of proteins to a queue accessible by the writer() function.

    :param seqDict: a dictionary containing DNA sequences as values and the names of these sequences as keys.
    :param minLen: the minimum length a protein must be to be included in the final output.
    :param procNum: the sequential number of the process for which this function has been called as the worker.
    :return:
    """
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
        seqToProteinNew.protCompletedQueue.put(1)
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
    """
    Called from seqToProteinNew() to remove the N character from the DNA sequence. These characters are ignored in
    this program.
    :param dnaSeq: the initial DNA sequence, which may have Ns at the start and/or end of the sequence.
    :return dnaSeq[fiveFrameCount: len(dnaSeq)-threeFrameCount]: the DNA sequence with Ns removed from start and end.
    """

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

def buildForwProt(seq, minLen):
    """
    Called by seqToProteinNew(), this function takes a DNA sequence and returns all the proteins which could be
    formed by the forward three frames of reading.
    :param seq: the DNA sequence to be converted to proteins using the forward three frame translation.
    :param minLen: the minimum length a protein sequence must be to be included in the final protein list.
    :return proteins: a list of proteins generated from the forward three frame translation.
    """

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
    """
    Called by seqToProteinNew(), this function takes a DNA sequence and returns all the proteins which could be
    formed by the reverse three frames of reading.
    :param seq: the DNA sequence to be converted to proteins using the reverse three frame translation.
    :param minLen: the minimum length a protein sequence must be to be included in the final protein list.
    :return proteins: a list of proteins generated from the reverse three frame translation.
    """
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

def createReverseSeq(dnaSeq):
    """
    Called by createReverseProt() to create the reversed sequence of each codon and ultimately conduct three frame
    translation in the reverse direction.
    :param dnaSeq: the small sequence of DNA which is to be reversed.
    :return reverseSeq: the reversed DNA sequence.
    """
    reverseDir = dnaSeq[::-1]
    _tab = str.maketrans(dict(zip('ATCG', 'TAGC')))
    reverseSeq = reverseDir.translate(_tab)
    return reverseSeq

def writer(queue, outputPath, removeSubFlag, writeSubFlag, originFlag):
    """
    Called as the worker function for the writerProcess from generateOutputNew(). This function takes the translated
    proteins that are put to queue and compiles them into a dictionary for writing to file. It ensures no duplication
    of output proteins and also initiates the generation of a file with subsequences removed if requested to do so by
    the user.

    :param queue: the queue (called toWriteQueue in each process) which each process pushes their final output to once
    six frame translation has been comleted.
    :param outputPath: the final path of the output fasta file.
    :param removeSubFlag: True if the user wants a second output file, which has sequences which are a subset of another
    longer sequence removed.
    :param writeSubFlag: True if removeSubFlag is True and the user wishes to create a fasta file containing all
    of the sequences which were subsets of others and were subsequently deleted from the second output.
    :param originFlag: True if removeSubFlag is True and the user wishes to ignore the origin data (data on which
    DNA sequence a protein was derived from) when creating the second output fasta file with subsets removed.
    :return:
    """
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

def createSeqObj(finalProteins):
    """
    Called by the writer() function, this function takes a dictionary of final proteins, converts all of them into
    SeqRecord objects and adds them to an iterable via the yield statement. This iterable is passed back to the
    SeqIO.write() method which writes the created SeqRecord to fasta file.

    :param finalProteins: a dictionary with protein sequences as keys and a list of origin DNA sequences as values.
    :return:
    """
    count = 1
    seqRecords = []

    for key, value in finalProteins.items():

        finalId = "dna|pro"+str(count)+';'

        # used to convey where the protein was derived from. We may need to do something similar
        for protein in value:
             finalId+=protein

        yield SeqRecord(Seq(key), id=finalId, description="")

        count += 1

def poolInitialiser(toWriteQueue, protCompletedQueue):
    """
    Called by generateOutputNew() before the multiprocessing pool is created to set up a global lock for a child
    processes and to give all processes access to important queues and shared variables.

    :param toWriteQueue: a multiprocessing.Queue() which the writer() function has access to. The output proteins of
    each process are pushed to this queue at the end of a process.
    :param protCompletedQueue: a queue which is used to track when processes are finished, ensuring that processes
    are only created as another process is finished.
    :return:
    """
    seqToProteinNew.toWriteQueue = toWriteQueue
    seqToProteinNew.protCompletedQueue = protCompletedQueue



