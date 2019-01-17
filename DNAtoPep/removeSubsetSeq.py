from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from time import time

# create sequence object taken directly from the Mers code.
def createSeqObj(seenPeptides, transFlag):
    """
    Given the set of matchedPeptides, converts all of them into SeqRecord objects and passes back a generator
    """
    count = 1
    seqRecords = []
    try:
        for sequence, name in seenPeptides.items():
            yield SeqRecord(Seq(sequence), id=name, description="")

    except AttributeError:
        for sequence in seenPeptides:
            finalId = "ipd|pep" + str(count) + ';'
            yield SeqRecord(Seq(sequence), id=finalId, description="")
            count += 1
    return seqRecords


# seenPeptides = []
# with open('100,000-Record.fasta', "rU") as handle:
#     counter = 0
#     for record in SeqIO.parse(handle, 'fasta'):
#         counter += 1
#         seenPeptides.append(str(record.seq))

time1 = time()

seenPeptides = []
with open('1mil.fasta', "rU") as handle:
    origNo = 0
    for record in SeqIO.parse(handle, 'fasta'):
        origNo += 1
        seenPeptides.append(str(record.seq))

time2 = time()

print('Number of original entries: ' + str(origNo))
print('Time to read all entries to list: ' + str(time2-time1))

seenPeptides.sort(key=len)
seenPeptides.reverse()

with open('sorted.fasta', "w") as output_handle:
    SeqIO.write(createSeqObj(seenPeptides, False), output_handle, "fasta")

time3 = time()
print('All sequences sorted and written to new dictionary in: ' + str(time3-time2))

seenPeptides = {}
with open('1mil.fasta', "rU") as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        seenPeptides[str(record.seq)] = record.name

time4 = time()
print('Second read of input file took: ' + str(time4-time3))

with open('sorted.fasta', "rU") as handle:
    counter = 1
    for record in SeqIO.parse(handle, 'fasta'):
        if counter % 10000 == 0:
            print(counter)
        counter+=1
        pep = str(record.seq)
        if pep not in seenPeptides.keys():
            continue
        for i in range(len(pep)):
            for j in range(i + 1, len(pep) + 1):
                if pep[i:j] in seenPeptides.keys() and pep[i:j] is not pep:
                    del seenPeptides[pep[i:j]]

print('Time to delete subset sequences: ' + str(time()-time4))
print('No. of sequences reduced from ' + str(origNo) + ' to ' + str(len(seenPeptides)))

with open('Output.fasta', "w") as output_handle:
    SeqIO.write(createSeqObj(seenPeptides, False), output_handle, "fasta")