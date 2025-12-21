MICRO = [
    "Apoe",
    "C1qa",
    "C1qb",
    "C1qc",
    "Cst3",
    "Csf1r",
    "Cx3cr1",
    "Trem2",
    "Tyrobp",
    "Lyz2",
    "Gpnmb",
    "P2ry12",
    "P2ry13",
    "Hexb",
    "Ctss",
    "Cd68",
    "Cd74",
]

ASTRO = [
    "Gfap",
    "Aqp4",
    "Slc1a3",
    "Aldh1l1",
    "Aldh1a2",
    "Clu",
    "Megf10",
    "Gja1",
    "Igf1",
    "Vim",
    "Serpina3n",
    "C4b",
    "C3",
]

OLIGO = ["Mbp", "Plp1", "Mobp", "Cnp", "Olig2", "Opalin", "Mog", "Spag16", "Igf2"]

ENDOTHELIAL = ["Kdr", "Pecam1", "Emcn", "Rgs5", "Acta2", "Pdgfra"]

NEURON = [
    "Slc17a7",
    "Slc17a6",
    "Gad1",
    "Gad2",
    "Pvalb",
    "Sst",
    "Ntsr2",
    "Neurod6",
    "Rorb",
    "Tle4",
    "Foxp2",
    "Calb1",
    "Calb2",
    "Chat",
    "Sncg",
]


def build_marker_sets(gene_cols):
    """
    Filter marker gene lists to keep only genes that exist in gene_cols.
    """
    return {
        "micro": [g for g in MICRO if g in gene_cols],
        "astro": [g for g in ASTRO if g in gene_cols],
        "oligo": [g for g in OLIGO if g in gene_cols],
        "endo": [g for g in ENDOTHELIAL if g in gene_cols],
        "neuron": [g for g in NEURON if g in gene_cols],
    }
