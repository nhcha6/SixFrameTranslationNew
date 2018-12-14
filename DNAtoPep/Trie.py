from typing import Tuple
import random
import string
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class TrieNode(object):
    """
    Our trie node implementation. Very basic. but does the job
    """

    def __init__(self, char: str):
        self.char = char
        self.children = []
        # Is it the last character of the word.`
        self.word_finished = False
        # How many times this character appeared in the addition process
        self.counter = 1


def add(root, word: str):
    """
    Adding a word in the trie structure
    """
    node = root
    for char in word:
        found_in_child = False
        # Search for the character in the children of the present `node`
        for child in node.children:
            if child.char == char:
                # We found it, increase the counter by 1 to keep track that another
                # word has it as well
                child.counter += 1
                # And point the node to the child that contains this char
                node = child
                found_in_child = True
                break
        # We did not find it so add a new chlid
        if not found_in_child:
            new_node = TrieNode(char)
            node.children.append(new_node)
            # And then point node to the new child
            node = new_node
    # Everything finished. Mark it as the end of a word.
    node.word_finished = True


def find_prefix(root, prefix: str) -> Tuple[bool, int]:
    """
    Check and return 
      1. If the prefix exsists in any of the words we added so far
      2. If yes then how may words actually have the prefix
    """
    node = root
    # If the root node has no children, then return False.
    # Because it means we are trying to search in an empty trie
    if not root.children:
        return False, 0
    for char in prefix:
        char_not_found = True
        # Search through all the children of the present `node`
        for child in node.children:
            if child.char == char:
                # We found the char existing in the child.
                char_not_found = False
                # Assign node as the child containing the char and break
                node = child
                break
        # Return False anyway when we did not find a char.
        if char_not_found:
            return False, 0
    # Well, we are here means we have found the prefix. Return true to indicate that
    # And also the counter of the last node. This indicates how many words have this
    # prefix
    return True, node.counter

def list_words(node, oldString = "", strings = []):
    for child in node.children:
        newString = oldString + child.char
        if child.children == []:
            strings.append(newString)
            continue
        list_words(child, newString, strings)
    return strings

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

root = TrieNode('*')

with open('Blah.fasta', "rU") as handle:
    # counter = 0
    for record in SeqIO.parse(handle, 'fasta'):
        # counter += 1
        seq = str(record.seq)
        add(root, seq)

seenPeptides = set(list_words(root))

with open('Output.fasta', "w") as output_handle:
    SeqIO.write(createSeqObj(seenPeptides, False), output_handle, "fasta")

# strings = []
# for i in range(0,100000):
#     st = ""
#     for j in range(0, random.randint(1,10)):
#         st += random.choice(string.ascii_uppercase)
#     strings.append(st)
# print(strings)
#
# for string in strings:
#     add(root, string)
#
# print(len(list_words(root)))