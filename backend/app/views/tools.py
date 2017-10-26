import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
import bandwagon as bw
import re
from base64 import b64encode, b64decode
from matplotlib.backends.backend_pdf import PdfPages
import flametree


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
            records = list(SeqIO.parse(stringio, fmt))
            if len(records) > 0:
                return (records, fmt)
        except:
            pass
    raise ValueError("Invalid sequence format")

def file_to_filelike_object(file_):
    content = file_.content.split("base64,")[1]
    return StringByteIO(b64decode(content))

def records_from_zip_file(zip_file):
    zip_file = flametree.file_tree(file_to_filelike_object(zip_file))
    records = []
    for f in zip_file._all_files:
        if f._extension.lower() in ['gb', 'fa']:
            new_records, fmt = string_to_record(f.read().decode('utf-8'))
            single_record = len(new_records) == 1
            for i, record in enumerate(new_records):
                name = record.id
                if name in [None, '', "<unknown id>"]:
                    number = ('' if single_record else ("%04d" % i))
                    name = f._name_no_extension + number
                record.id = name
            records += new_records
    return records


def records_from_data_file(data_file):
    content = data_file.content.split("base64,")[1]
    content = b64decode(content).decode("utf-8")
    records, fmt = string_to_record(content)
    return records, fmt

def records_from_data_files(data_files):
    records = []
    for file_ in data_files:
        if file_.name.lower().endswith('zip'):
            records += records_from_zip_file(file_)
            continue
        recs, fmt = records_from_data_file(file_)
        single_record = len(recs) == 1
        for i, record in enumerate(recs):
            name = record.id
            if name in [None, '', "<unknown id>"]:
                name = file_.name + ('' if single_record else ("%04d" % i))
            record.id = name
        records += recs
    return records

def data_to_html_data(data, datatype):
    datatype = {
        'zip': 'application/zip',
        'genbank': 'application/genbank'

    }.get(datatype, datatype)
    return 'data:%s;base64,%s' % (datatype, b64encode(data).decode("utf-8"))

def zip_data_to_html_data(data):
    return data_to_html_data(data, 'application/zip')

LADDERS = {
   "100_to_4k": bw.ladders.LADDER_100_to_4k
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

def figures_to_pdf_report_data(figures, filename='report.pdf'):
    pdf_io = BytesIO()
    with PdfPages(pdf_io) as pdf:
        for fig in figures:
            pdf.savefig(fig, bbox_inches="tight")
    return {
        'data': ('data:application/pdf;base64,' +
                 b64encode(pdf_io.getvalue()).decode("utf-8")),
        'name': filename,
        'mimetype': 'application/pdf'
    }
