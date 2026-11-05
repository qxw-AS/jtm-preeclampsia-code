#!/usr/bin/env python3
"""Build MAGMA --set-annot file for placental cell-type marker sets.

Maps marker symbols -> Entrez (via hg19 NCBI37.3.gene.loc), writes a set-annot
file (one gene set per line: SETID <entrez...>), and reports unmapped markers.
"""
GENE_LOC = "/workspace/data/annotation/NCBI37.3.gene.loc"
OUT = "/workspace/work/magma_celltype_sets.txt"

MARKERS = {
    "Syncytiotrophoblast": ["CGA", "CGB", "CGB3", "PSG1", "PSG3", "PSG4", "CSH1", "CSH2",
                            "ERVFRD-1", "GCM1", "MUC15", "CYP19A1", "TFAP2A"],
    "Cytotrophoblast": ["KRT7", "GATA3", "TEAD4", "EGFR", "MSX2", "TP63", "PEG10", "HAND1"],
    "Extravillous_trophoblast": ["HLA-G", "ITGA1", "MMP2", "PAPPA2", "ASCL2", "NOTUM",
                                 "MMP11", "ADAM12", "PLAC8"],
    "Stromal_mesenchymal": ["COL1A1", "COL3A1", "COL1A2", "COL5A1", "DCN", "LUM",
                            "PDGFRA", "ACTA2", "FN1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "ENG", "PLVAP", "FLT1", "CD34"],
    "Macrophage_Hofbauer": ["CD68", "CD14", "CSF1R", "LYZ", "CD163", "MRC1", "MSR1", "AIF1"],
    "T_NK": ["CD3D", "CD3E", "NKG7", "KLRD1", "GNLY", "CCL5", "CD2", "TRAC"],
    "B_Plasma": ["CD79A", "MS4A1", "IGHG1", "MZB1", "JCHAIN", "CD79B"],
    "Erythrocyte_fetal": ["HBB", "HBA1", "HBA2", "ALAS2"],
}

# symbol -> entrez (hg19)
sym2entrez = {}
with open(GENE_LOC) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 6:
            sym2entrez[f[-1]] = f[0]

lines, missing = [], {}
for ct, genes in MARKERS.items():
    ids = []
    for g in genes:
        e = sym2entrez.get(g)
        if e:
            ids.append(e)
        else:
            missing.setdefault(ct, []).append(g)
    lines.append(f"{ct} " + " ".join(ids))
    print(f"  {ct:32s} mapped {len(ids)}/{len(genes)}")

with open(OUT, "w") as out:
    out.write("\n".join(lines) + "\n")
print(f"[write] {OUT}")
if missing:
    print("[unmapped]", missing)
