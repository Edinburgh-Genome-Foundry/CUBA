import sys
from io import StringIO, BytesIO
import re
from base64 import b64encode, b64decode
from copy import deepcopy

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
from matplotlib.backends.backend_pdf import PdfPages

import flametree
from snapgene_reader import snapgene_file_to_seqrecord
import bandwagon as bw
import crazydoc
from dnachisel.biotools import sequence_to_biopython_record

from fuzzywuzzy import process
import pandas


def did_you_mean(name, other_names, limit=5, min_score=50):
    results = process.extract(name, list(other_names), limit=limit)
    return [e for (e, score) in results if score >= min_score]


def fix_ice_genbank(genbank_txt):
    lines = genbank_txt.splitlines()
    lines[0] += max(0, 80 - len(lines[0])) * " "
    return "\n".join(lines)


def write_record(record, target, fmt="genbank"):
    """Write a record as genbank, fasta, etc. via Biopython, with fixes"""
    record = deepcopy(record)
    if fmt == "genbank":
        if isinstance(record, (list, tuple)):
            for r in record:
                r.name = r.name[:20]
        else:
            record.name = record.name[:20]
    if hasattr(target, "open"):
        target = target.open("w")
    SeqIO.write(record, target, fmt)


def autoname_genbank_file(record):
    return record.id.replace(".", "_") + ".gb"


def string_to_record(string):
    """Convert a string of a fasta, genbank... into a simple ATGC string.

    Can also be used to detect a format.
    """
    matches = re.match("([ATGC][ATGC]*)", string)
    # print("============", len(matches.groups()[0]), len(string))
    # print (matches.groups()[0] == string)
    if (matches is not None) and (matches.groups()[0] == string):
        return SeqRecord(Seq(string, DNAAlphabet())), "ATGC"

    for fmt in ("fasta", "genbank"):
        if fmt == "genbank":
            string = fix_ice_genbank(string)
        try:
            stringio = StringIO(string)
            records = list(SeqIO.parse(stringio, fmt))
            if len(records) > 0:
                return (records, fmt)
        except:
            pass
    try:
        record = snapgene_file_to_seqrecord(filecontent=StringIO(string))
        return record
    except:
        pass
    raise ValueError("Invalid sequence format")


def file_to_filelike_object(file_, type="byte"):
    content = file_.content.split("base64,")[1]
    filelike = BytesIO if (type == "byte") else StringIO
    return filelike(b64decode(content))


def spreadsheet_file_to_dataframe(filedict, header="infer"):
    filelike = file_to_filelike_object(filedict)
    if filedict.name.endswith(".csv"):
        return pandas.read_csv(filelike, header=header)
    else:
        return pandas.read_excel(filelike, header=header)


def records_from_zip_file(zip_file):
    zip_file = flametree.file_tree(file_to_filelike_object(zip_file))
    records = []
    for f in zip_file._all_files:
        ext = f._extension.lower()
        if ext in ["gb", "gbk", "fa", "dna"]:
            try:
                new_records, fmt = string_to_record(f.read())
            except:
                content_stream = BytesIO(f.read("rb"))
                try:
                    record = snapgene_file_to_seqrecord(
                        fileobject=content_stream
                    )
                    new_records, fmt = [record], "snapgene"
                except:
                    try:
                        parser = crazydoc.CrazydocParser(
                            ["highlight_color", "bold", "underline"]
                        )
                        new_records = parser.parse_doc_file(content_stream)
                        fmt = "doc"
                    except:
                        raise ValueError(
                            "Format not recognized for file " + f._path
                        )

            single_record = len(new_records) == 1
            for i, record in enumerate(new_records):
                name = record.id
                if name in [
                    None,
                    "",
                    "<unknown id>",
                    ".",
                    " ",
                    "<unknown name>",
                ]:
                    number = "" if single_record else ("%04d" % i)
                    name = f._name_no_extension.replace(" ", "_") + number
                record.id = name
                record.name = name
                record.file_name = f._name_no_extension
            records += new_records
    return records


def records_from_data_file(data_file):
    content = b64decode(data_file.content.split("base64,")[1])
    try:
        records, fmt = string_to_record(content.decode("utf-8"))
    except:
        try:
            record = snapgene_file_to_seqrecord(fileobject=BytesIO(content))
            records, fmt = [record], "snapgene"
        except:
            try:
                parser = crazydoc.CrazydocParser(
                    ["highlight_color", "bold", "underline"]
                )
                records = parser.parse_doc_file(BytesIO(content))
                fmt = "doc"
            except:
                try:
                    df = spreadsheet_file_to_dataframe(data_file, header=None)
                    records = [
                        sequence_to_biopython_record(
                            sequence=seq, id=name, name=name
                        )
                        for name, seq in df.values
                    ]
                    fmt = "spreadsheet"
                except:
                    raise ValueError(
                        "Format not recognized for file " + data_file.name
                    )
    if not isinstance(records, list):
        records = [records]
    return records, fmt


def record_to_formated_string(record, fmt="genbank", remove_descr=False):
    if remove_descr:
        record = deepcopy(record)
        if isinstance(record, (list, tuple)):
            for r in record:
                r.description = ""
        else:
            record.description = ""
    fileobject = StringIO()
    write_record(record, fileobject, fmt)
    return fileobject.getvalue().encode("utf-8")


def records_from_data_files(data_files):
    records = []
    for file_ in data_files:
        circular = ("circular" not in file_) or file_.circular
        if file_.name.lower().endswith("zip"):
            records += records_from_zip_file(file_)
            continue
        recs, fmt = records_from_data_file(file_)
        single_record = len(recs) == 1
        for i, record in enumerate(recs):
            record.circular = circular
            record.linear = not circular
            name_no_extension = "".join(file_.name.split(".")[:-1])
            name = name_no_extension + ("" if single_record else ("%04d" % i))
            name = name.replace(" ", "_")
            UNKNOWN_IDS = [
                "None",
                "",
                "<unknown id>",
                ".",
                "EXPORTED",
                "<unknown name>",
                "Exported",
            ]
            record.seq.alphabet = DNAAlphabet()
            # Sorry for this parts, it took a lot of "whatever works".
            # keep your part names under 20c and pointless, and everything
            # will be good
            if str(record.id).strip() in UNKNOWN_IDS:
                record.id = name
            if str(record.name).strip() in UNKNOWN_IDS:
                record.name = name
            record.file_name = name_no_extension
        records += recs
    return records


def data_to_html_data(data, datatype, filename=None):
    """Data types: zip, genbank, fasta, pdf"""
    datatype = {
        "zip": "application/zip",
        "genbank": "application/genbank",
        "fasta": "application/fasta",
        "pdf": "application/pdf",
        "xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    }.get(datatype, datatype)
    datatype = "data:%s;" % datatype
    data64 = "base64,%s" % b64encode(data).decode("utf-8")
    headers = ""
    if filename is not None:
        headers += "headers=filename%3D" + filename + ";"
    return datatype + headers + data64


def zip_data_to_html_data(data):
    return data_to_html_data(data, "application/zip")


LADDERS = {"100_to_4k": bw.ladders.LADDER_100_to_4k}


def matplotlib_figure_to_svg_base64_data(fig, **kwargs):
    """Return a string of the form 'data:image/svg+xml;base64,XXX' where XXX
    is the base64-encoded svg version of the figure."""
    output = BytesIO()
    fig.savefig(output, format="svg", **kwargs)
    svg_txt = output.getvalue().decode("utf-8")
    svg_txt = "\n".join(svg_txt.split("\n")[4:])
    svg_txt = "".join(svg_txt.split("\n"))

    content = b64encode(svg_txt.encode("utf-8"))
    result = (b"data:image/svg+xml;base64," + content).decode("utf-8")

    return result


def matplotlib_figure_to_bitmap_base64_data(fig, fmt="png", **kwargs):
    """Return a string of the form 'data:image/png;base64,XXX' where XXX
    is the base64-encoded svg version of the figure."""
    output = BytesIO()
    fig.savefig(output, format=fmt, **kwargs)
    bitmap = output.getvalue()
    content = b64encode(bitmap)
    result = (
        b"data:image/%s;base64,%s" % (fmt.encode("utf-8"), content)
    ).decode("utf-8")
    return result


def figures_to_pdf_report_data(figures, filename="report.pdf"):
    pdf_io = BytesIO()
    with PdfPages(pdf_io) as pdf:
        for fig in figures:
            pdf.savefig(fig, bbox_inches="tight")
    return {
        "data": (
            "data:application/pdf;base64,"
            + b64encode(pdf_io.getvalue()).decode("utf-8")
        ),
        "name": filename,
        "mimetype": "application/pdf",
    }


def csv_to_list(csv_string, sep=","):
    return [
        element.strip()
        for line in csv_string.split("\n")
        for element in line.split(sep)
        if len(element.strip())
    ]


def set_record_topology(record, topology):
    """Set the Biopython record's topology, possibly passing if already set.

    This actually sets the ``record.annotations['topology']``.The ``topology``
    parameter can be "circular", "linear", "default_to_circular" (will default
    to circular if ``annotations['topology']`` is not already set) or
    "default_to_linear".
    """
    valid_topologies = [
        "circular",
        "linear",
        "default_to_circular",
        "default_to_linear",
    ]
    if topology not in valid_topologies:
        raise ValueError(
            "topology (%s) should be one of %s."
            % (topology, ", ".join(valid_topologies))
        )
    annotations = record.annotations
    default_prefix = "default_to_"
    if topology.startswith(default_prefix):
        if "topology" not in annotations:
            annotations["topology"] = topology[len(default_prefix) :]
    else:
        annotations["topology"] = topology

