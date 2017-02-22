import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
import bandwagon as bw
import re

PYTHON3 = (sys.version_info > (3, 0))
if PYTHON3:
    from io import StringIO
else:
    from StringIO import StringIO

def string_to_record(string):
    """Convert a string of a fasta, genbank... into a simple ATGC string.

    Can also be used to detect a format.

    """
    matches = re.match("([ATGC][ATGC]*)", string)
    if (matches is not None) and (matches.groups()[0] == string):
        return SeqRecord(Seq(string, DNAAlphabet()), "ATGC")

    for fmt in ("fasta", "genbank"):
        try:
            stringio = StringIO(string)
            return (SeqIO.read(stringio, fmt), fmt)
        except:
            pass
    raise ValueError("Invalid sequence format")

LADDERS = {
   "100-4k": bw.ladders.LADDER_100_to_4k
}
