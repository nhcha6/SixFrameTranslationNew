from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from time import time

def createSeqObj(matchedPeptides, transFlag):
    """
    Given the set of matchedPeptides, converts all of them into SeqRecord objects and passes back a generator
    """
    count = 1
    seqRecords = []
    for sequence in matchedPeptides:

        finalId = "ipd|pep"+str(count)+';'

        yield SeqRecord(Seq(sequence), id=finalId, description="")

        count += 1

    return seqRecords

seenPeptides = []
with open('blahh2.fasta', "rU") as handle:
    counter = 0
    for record in SeqIO.parse(handle, 'fasta'):
        counter += 1
        seenPeptides.append(str(record.seq))

start = time()
seenPeptides.sort(key=len)
seenPeptides.reverse()

with open('sorted.fasta', "w") as output_handle:
    SeqIO.write(createSeqObj(seenPeptides, False), output_handle, "fasta")

seenPeptides = set(seenPeptides)

with open('sorted.fasta', "rU") as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        print(len(seenPeptides))
        pep = str(record.seq)
        if pep not in seenPeptides:
            continue
        for i in range(len(pep)):
            for j in range(i + 1, len(pep) + 1):
                if pep[i:j] in seenPeptides and pep[i:j] is not pep:
                    seenPeptides.remove(pep[i:j])

print(time()-start)
print(len(seenPeptides))

with open('Output.fasta', "w") as output_handle:
    SeqIO.write(createSeqObj(seenPeptides, False), output_handle, "fasta")