from copy import deepcopy
import sys
from collections import defaultdict
from base64 import b64encode

from Bio import Restriction, SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

import dna_features_viewer as dfw

import bandwagon as bw

import flametree

PYTHON3 = (sys.version_info > (3, 0))
if PYTHON3:
    from io import BytesIO
    StringByteIO = BytesIO
else:
    from StringIO import StringIO
    StringByteIO = StringIO

def annotate_record(seqrecord, location="full", feature_type="misc_feature",
                margin=0, **qualifiers):
    """Add a feature to a Biopython SeqRecord. (also returns that same record)
    """
    if location == "full":
        location = (margin, len(seqrecord)-margin)

    strand = location[2] if (len(location) == 3) else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type
        )
    )
    return seqrecord

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


class DigestGraphicTranslator(dfw.BiopythonTranslator):

    def compute_feature_color(self, feature):
        if "band_label" in feature.qualifiers:
            return "#e5d740"
        elif "source" in feature.qualifiers:
            return "#ff8282"
        else:
            return "#eaedff"


def annotate_digestion_bands(record, enzymes, ladder):

    linear = record.linear if hasattr(record, 'linear') else False
    batch = Restriction.RestrictionBatch(enzymes)
    cuts_dict = batch.search(record.seq)
    all_cuts = sorted(
        set([0, len(record)] + [c for cc in cuts_dict.values() for c in cc]))
    bands = list(zip(all_cuts, all_cuts[1:]))
    if (not linear) and len(bands) > 1:
        start, end = bands.pop()
        band0 = [-(end - start), bands[0][1]]
        if bands == []:
            bands = [band0]
        else:
            bands[0] = band0
    sorted_bands = sorted(bands, key=lambda b: b[0] - b[1])
    new_record = deepcopy(record)
    for (band, label) in zip(sorted_bands, "abcdefghijkl"):
        band_size = abs(band[1] - band[0])
        formatted_size = bw.Band._format_dna_size(band_size)
        annotate_record(
            new_record,
            location=band,
            label="%s - %s" % (label, formatted_size),
            feature_type="misc_feature", band_label=label,
            band_size=band_size
        )
    return new_record


def plot_record_digestion(record_digestion, ladder,
                          record_label, digestion_label):
    gs = gridspec.GridSpec(4, 10)

    ax_bands, ax_record, ax_digest = [plt.subplot(
        z) for z in (gs[:, 0], gs[:3, 3:], gs[3, 3:])]
    ax_digest.figure.set_size_inches(20, 7)
    ax_digest.set_title("BANDS")
    ax_record.set_title(record_label, fontsize=22)

    bands = sorted([
        (feature.qualifiers["band_label"], feature.qualifiers[
         "band_size"], feature.location.start)
        for feature in record_digestion.features
        if feature.qualifiers.get("band_label", False)
    ])

    for _, _, start in bands:
        for ax in (ax_digest, ax_record):
            ax.axvline(start, ls=":", color="k", lw=0.5)

    gr_record, gr_digestion = [
        DigestGraphicTranslator([fl]).translate_record(record_digestion)
        for fl in [lambda f: ("band_label" not in f.qualifiers)
                   and ("homology" not in f.qualifiers.get("label", [])),
                   lambda f: ("band_label" in f.qualifiers)]
    ]
    gr_digestion.split_overflowing_features_circularly()
    for f in gr_digestion.features:
        ax_record.axvline(f.start, ls=":", color="k", lw=0.5)
    gr_digestion.plot(ax=ax_digest)
    gr_record.plot(ax=ax_record, with_ruler=False)

    pattern = bw.BandsPattern([bw.Band(dnasize, ladder=ladder, label=label)
                               for label, dnasize, _ in bands],
                              ladder=ladder)
    patternset = bw.BandsPatternsSet([pattern], ladder=ladder, ladder_ticks=3,
                                     ticks_fontdict=dict(size=12),
                                     label=digestion_label,
                                     label_fontdict=dict(size=18))
    patternset.plot(ax_bands)
    return (ax_bands, ax_record, ax_digest)


def generate_report(records, digestions, ladder, group_by="digestions",
                    full_report=True):
    """Yeah !"""

    files_contents = []
    all_patterns = defaultdict(lambda *a: {})
    zip_root = flametree.file_tree("@memory")
    # with PdfPages(details_pdf_sio) as details_pdf:
    with PdfPages(zip_root._file("Details.pdf").open("wb")) as details_pdf:
        for record in records:
            record_label = record.id
            for enzymes in digestions:

                enzymes_label = " + ".join(sorted(enzymes))
                basename = "%s--%s" % (record_label.replace(" ", "_"),
                                       "+".join(enzymes))
                record_digestion = annotate_digestion_bands(
                    record, enzymes, ladder)
                files_contents.append([
                    ("genbanks", basename + ".gb"),
                    lambda fh: SeqIO.write(record_digestion, fh, "genbank")
                ])
                if full_report:
                    (ax, _, _) = plot_record_digestion(
                        record_digestion, ladder, record_label, enzymes_label)
                    details_pdf.savefig(ax.figure, bbox_inches="tight")
                    plt.close(ax.figure)
                bands = sorted([
                    f.qualifiers["band_size"]
                    for f in record_digestion.features
                    if f.qualifiers.get("band_size", False)
                ])
                if group_by == "digestions":
                    all_patterns[enzymes_label][record_label] = bands
                else:
                    all_patterns[record_label][enzymes_label] = bands

    Y = len(all_patterns)
    X = len(list(all_patterns.values())[0])
    fig, axes = plt.subplots(Y, 1, figsize=(0.9 * X, 3 * Y))
    if Y == 1:
        axes = [axes]
    for ax, (cat1, cat2s) in zip(axes, sorted(all_patterns.items())):
        pattern_set = bw.BandsPatternsSet(
            patterns=[
                bw.BandsPattern(
                    _bands,
                    ladder=ladder,
                    label=cat2 if (ax == axes[0]) else None,
                    label_fontdict=dict(rotation=70),
                    global_bands_props={"band_thickness": 2.5})
                for cat2, _bands in cat2s.items()
            ],
            ladder=ladder,
            ladder_ticks=4,
            ticks_fontdict=dict(size=9),
            label=cat1
        )
        pattern_set.plot(ax)

    preview = matplotlib_figure_to_svg_base64_data(fig, bbox_inches="tight")
    if full_report:
        # zip_root._file("details.pdf").write(details_pdf_sio.getvalue())
        fig.savefig(zip_root._file("summary.pdf").open("wb"), format="pdf",
                    bbox_inches="tight")
        report = zip_root._close()
    else:
        report = None
    return preview, report
