# Published library

The SmartSepsis-Oph multimodal AMR variants dataset is publicly available on
Hugging Face:

> **<https://huggingface.co/datasets/JVLegend/smartsepsis-oph>**

License: **CC BY 4.0**.

A Zenodo DOI is planned alongside the framework preprint:

> Zenodo DOI: *to be assigned*.

## Two configs

| Config | Rows | Size | Contents |
|---|---|---|---|
| `panel` | 45 (43 with full pipeline) | ~3.2 MB | Multimodal: DNA + protein + ESM-2 + ProtT5 + PDB + structural features + drug class. |
| `extended` | **9 034** | ~34 MB | AMRFinderPlus catalog (8 991) + panel (43): variant_id + source + drug_classes + ESM-2 embedding. |

## Coverage (panel)

| Family | Variants | Drug class |
|---|---|---|
| **mecA** | mecA1, mecA2 | penam, cephalosporin, methicillin (MRSA) |
| **blaKPC** | KPC-2, 3, 4, 5, 11, 30, 31 | carbapenem, cephalosporin, penam |
| **blaNDM** | NDM-1, 2, 5, 7 | carbapenem, cephalosporin, penam |
| **blaOXA-48** | OXA-48, 181, 232 | carbapenem, cephalosporin, penam |
| **blaVIM** | VIM-1, 2, 4 | carbapenem, cephalosporin, penam |
| **blaIMP** | IMP-1, 6 | carbapenem, cephalosporin, penam |
| **blaGES** | GES-1, 5 | carbapenem, cephalosporin |
| **blaCTX-M** | CTX-M-2, 8, 9, 14, 27 | cephalosporin, penam (ESBL) |
| **vanA** | vanA | glycopeptide (vancomycin) |
| **mcr** | mcr-1, mcr-1.1, mcr-5 | polymyxin, peptide |
| **qnrS** | qnrS1, qnrS2 | fluoroquinolone |
| **armA** | armA | aminoglycoside |

## Schema (panel)

Each row contains:

- `variant_id`, `gene_family`, `dna_accession`.
- `dna_sequence`, `protein_sequence`, `protein_length`.
- `drug_classes`, `resistance_mechanism`.
- `esm2_embedding` (640d), `esm2_model`.
- `prott5_embedding` (1024d), `prott5_model`.
- `structure_pdb` (full PDB text), `structure_source`, `struct_length`.
- Structural descriptors: `struct_rg`, `struct_compactness`,
  `struct_contact_density`, `struct_mean_plddt`.

For the extended config, only `variant_id`, `source`, `drug_classes`, and
`esm2_embedding` are guaranteed.

## How to load

```python
from datasets import load_dataset
import numpy as np

# Multimodal panel (curated)
ds = load_dataset("JVLegend/smartsepsis-oph", "panel", split="train")
row = ds[0]
emb = np.concatenate([
    np.array(row["esm2_embedding"]),
    np.array(row["prott5_embedding"]),
])  # 1664d ensemble

# Extended (AMRFinderPlus + panel, ESM-2 only)
ds_ext = load_dataset("JVLegend/smartsepsis-oph", "extended", split="train")
```

## Construction (summary)

1. Variants pulled from NCBI RefSeq via Entrez, curated against
   AMRFinderPlus / CARD.
2. Translation: longest ORF, bacterial table 11.
3. ESM-2 and ProtT5 mean-pooled embeddings.
4. 3D structure: ESMFold / ColabFold AF2 / AlphaFold Server AF3 (tiered).
5. Structure descriptors from Cα coordinates.
6. Drug class / mechanism from CARD ontology.

See the [dataset README on Hugging Face](https://huggingface.co/datasets/JVLegend/smartsepsis-oph)
for the full construction notes.

## License notes

- Data licensed CC BY 4.0.
- **AlphaFold Server (AF3) terms** apply to the mecA1 / mecA2 PDBs — research
  use only, no commercial redistribution. Regenerate via ColabFold AF2 or
  AlphaFold 2 / OpenFold for commercial use.

## See also

- [Reproducibility](reproducibility.md).
- [Structural analysis](../methods/structural.md).
- [Citing](../citing.md).
