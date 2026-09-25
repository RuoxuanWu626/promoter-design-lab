"""Cell-type-specific regulatory elements — a second, separate library.

Why this exists
---------------
Puffin's ten motifs are **core promoter** elements: TATA, the initiator, the
GC-box, CCAAT, NRF1, ETS, CREB, YY1, ZNF143, the U1 site. They largely set how
much initiation a promoter produces and where it starts, and they do that in
much the same way across cell types. They are not what makes a promoter
*cell-type specific*.

What drives cell-type specificity is other sequence in the promoter region:
binding sites for lineage-determining transcription factors. Those are a
different set of motifs and they live alongside, not inside, the core promoter.

So the app carries two libraries:

* ``motifs.get_library()`` — Puffin's ten core promoter motifs, learned from
  the real checkpoint. These drive the Puffin profile models.
* this module — lineage TF sites, which drive the cell-type model.

Keeping them separate is the point. Mode 1's attribution panel compares how
much of the *activity* attribution versus the *specificity* attribution falls
on each, and the honest answer is that they come apart.

PROVENANCE: the consensus strings below are the standard published consensus
for each factor, and the cell-type assignments follow the usual lineage
associations. They are a caricature at the level of "GATA works in erythroid
cells", not a quantitative model, and they drive the MOCK cell-type adapter
only. The real AlphaGenome adapter ignores all of this and reads its own CAGE
tracks.
"""

from __future__ import annotations

from motifs import Motif

# id, name, consensus, colour, factors, cell types it activates, cell types it represses
_SPEC = [
    ("gata",   "GATA", "AGATAAGA", "#d64550",
     "GATA1/GATA2", ["K562"], []),
    ("hnf4",   "HNF4A / DR1", "AGGTCAAAGGTCA", "#e09f3e",
     "HNF4A, NR2F2", ["HepG2", "Hepatocyte"], []),
    ("cebp",   "C/EBP", "TTGCGCAAT", "#c1666b",
     "CEBPA/CEBPB", ["HepG2", "A549", "Hepatocyte"], []),
    ("spi1",   "PU.1 / SPI1", "AAAGAGGAAGTG", "#4f9d69",
     "SPI1", ["GM12878"], []),
    ("ebf1",   "EBF1", "TCCCAAGGGA", "#2a9d8f",
     "EBF1", ["GM12878"], []),
    ("sox_oct", "SOX2–OCT4", "TTTGCATAACAA", "#8367c7",
     "POU5F1/SOX2", ["H1"], []),
    # The bare CAGCTG core is 6 bp and turns up roughly once per kb of random
    # sequence, which would flood the background with chance specificity
    # sites. Written with flanking context, as the bound site actually is.
    ("ebox_neuro", "E-box (bHLH)", "GCCAGCTGGC", "#5465ff",
     "NEUROD1, ASCL1", ["SK-N-SH"], []),
    ("ere",    "ERE", "AGGTCACAGTGACCT", "#e36414",
     "ESR1", ["MCF-7"], []),
    ("tead",   "TEAD / MCAT", "GGCATTCCAG", "#3d9dd6",
     "TEAD1/YAP", ["HUVEC", "HeLa-S3"], []),
    ("runx",   "RUNX", "TGTGGTTT", "#9a6fb0",
     "RUNX1/RUNX3", ["HCT116", "GM12878"], []),
    ("tp53",   "p53", "AGACATGCCTAGACATGCCT", "#b5179e",
     "TP53", ["HCT116", "IMR-90"], []),
    ("ap1",    "AP-1 / TRE", "ATGACTCAT", "#f28f3b",
     "FOS/JUN", ["A549", "IMR-90", "HeLa-S3"], []),
    ("nfkb",   "NF-κB", "GGGACTTTCC", "#d90429",
     "RELA/NFKB1", ["GM12878", "HUVEC"], []),
    # REST/NRSF silences neuronal genes in NON-neuronal cells. Neuronal cells
    # express little REST, so an NRSE site represses almost everywhere except
    # there. That makes it the textbook way to win specificity by silencing
    # off-targets rather than by raising the target, which is the distinction
    # Mode 1 is built to expose.
    ("rest",   "REST / NRSE", "TTCAGCACCACGGACAGCGCC", "#6c757d",
     "REST/NRSF", [],
     ["K562", "HepG2", "GM12878", "HeLa-S3", "A549", "MCF-7",
      "IMR-90", "HCT116", "HUVEC", "Hepatocyte"]),
]

CELLTYPE_ELEMENT_IDS = [s[0] for s in _SPEC]


def build_elements() -> dict[str, Motif]:
    """Return {id: Motif} for the cell-type element library."""
    out: dict[str, Motif] = {}
    for eid, name, consensus, color, factors, up, down in _SPEC:
        targets = ", ".join(up) if up else "none"
        represses = f"; represses {', '.join(down)}" if down else ""
        out[eid] = Motif(
            id=eid,
            name=name,
            aliases=[factors],
            consensus=consensus,
            # These are enhancer-like sites: they act over a broad window
            # rather than from a fixed distance, so there is no single
            # canonical offset. -120 just puts them upstream of the core.
            typical_offset=-120,
            strand_specific=False,
            color=color,
            notes=(f"Binding site for {factors}. Activates {targets}{represses}. "
                   "A lineage TF site, not a core promoter motif — this is the "
                   "kind of element that makes a promoter cell-type specific."),
            # Mock parameters. Deliberately small amplitude: these sites are
            # not what sets bulk initiation strength.
            mock_amplitude=0.22,
            mock_peak_shift=110,
            mock_width=40.0,
            mock_broad_fraction=0.75,
            mock_strand_asymmetry=0.9,
        )
        out[eid].kind = "celltype_element"
        out[eid].factors = factors
        out[eid].activates = list(up)
        out[eid].represses = list(down)
        out[eid].source = "literature_consensus"
    return out


def activation_map() -> dict[str, dict[str, float]]:
    """{element_id: {cell_type: weight}} for the mock cell-type adapter.

    Strongly positive in the lineages the factor belongs to, mildly negative
    where it is repressive, near zero elsewhere. This is what makes tau move.
    """
    out: dict[str, dict[str, float]] = {}
    for eid, _n, _c, _col, _f, up, down in _SPEC:
        w: dict[str, float] = {}
        for ct in up:
            w[ct] = 1.0
        for ct in down:
            w[ct] = -0.85
        out[eid] = w
    return out
