import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
import bandwagon as bw
import re
from base64 import b64encode

PYTHON3 = (sys.version_info > (3, 0))
if PYTHON3:
    from io import StringIO, BytesIO
    StringByteIO = BytesIO
else:
    from StringIO import StringIO
    StringByteIO = StringIO

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

def matplotlib_figure_to_svg_base64_data(fig, **kwargs):
    """Return a string of the form 'data:image/svg+xml;base64,XXX' where XXX
    is the base64-encoded svg version of the figure."""
    output = StringByteIO()
    fig.savefig(output, format='svg', **kwargs)
    svg_txt = output.getvalue().decode("ascii")
    svg_txt = "\n".join(svg_txt.split("\n")[4:])
    svg_txt = "".join(svg_txt.split("\n"))

    if PYTHON3:
        content = b64encode(svg_txt.encode("ascii"))
        result = (b"data:image/svg+xml;base64," + content).decode("utf-8")

        return result
    else:
        return "data:image/svg+xml;base64," + b64encode(svg_txt)
