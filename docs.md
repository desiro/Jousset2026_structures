# RNA Interaction Detection Pipeline (SHAPE-guided)

This repository contains two Python scripts for detecting and comparing short inter-molecular RNA–RNA interactions using sequence information, SHAPE reactivity data, and thermodynamic modeling via the ViennaRNA package.

---

## Overview

The pipeline identifies short, perfectly base-paired RNA duplexes between different RNA sequences (e.g. viral genome segments), filters them based on structural and thermodynamic criteria, and optionally compares interaction sets across different datasets.

Two main scripts are provided:

- **`getStructuresSHAPE.py`** — detects and annotates RNA–RNA interactions  
- **`compareStructuresSHAPE.py`** — compares interaction sets between two datasets  

---

## `getStructuresSHAPE.py`

### Purpose

This script detects candidate inter-molecular RNA interactions between sequences by combining k-mer scanning, SHAPE-based structural filtering, constrained RNA cofolding, thermodynamic evaluation, positional filtering, and redundancy removal.

---

### Workflow

#### 1. K-mer extraction

RNA sequences are scanned using a sliding window to extract k-mers (default: 4–20 nt). For each k-mer, the corresponding SHAPE reactivity values are retrieved.

Only k-mers where all nucleotides have SHAPE reactivity ≤ 0.8 are retained, selecting structurally constrained regions.

---

#### 2. Duplex candidate identification

K-mers of equal length from different sequences are compared pairwise. Each pair is evaluated using constrained RNA cofolding with the ViennaRNA package.

Folding is performed with:

- temperature: **37.0 °C (default)**  
- dangling end model: **2**  
- no lonely base pairs  

Hard constraints enforce exclusive intermolecular pairing, preventing intra-strand interactions.

Only perfectly base-paired duplexes (no mismatches or loops) are retained.

---

#### 3. Thermodynamic filtering

Candidate duplexes are re-evaluated using partition function–based dimer ensemble calculations.

Interactions are retained if:

- ensemble free energy **ΔG ≤ −10.0 kcal/mol (default)**  

This ensures thermodynamic stability beyond the minimum free energy structure.

---

#### 4. Positional filtering

To account for structural organization (e.g. proximity of sequence ends), the relative position of interacting regions is evaluated.

Interactions are retained only if the difference in distance from the sequence ends between both partners does not exceed a defined threshold (default: 200 nt).

---

#### 5. Redundancy removal

Shorter interactions fully contained within longer interactions (on both sequences) are removed, resulting in a non-redundant set of maximal interaction regions.

---

#### 6. Annotation

Each interaction is annotated with:

- sequence identifiers  
- genomic positions (start and end)  
- interaction length (k-mer size)  
- free energy (ΔG)  
- SHAPE reactivity values and mean SHAPE  
- distance-to-end metrics  
- interacting RNA sequences  
- simplified duplex structure  

---

## `compareStructuresSHAPE.py`

### Purpose

This script compares two sets of RNA interaction junctions (e.g. from different genomes or conditions) and identifies shared and dataset-specific interactions.

---

### Workflow

#### 1. Grouping by interaction length

Junctions from both input files are grouped by k-mer size to ensure comparisons are made between interactions of equal length.

---

#### 2. Exact matching

Interactions are compared based on:

- sequence pair (aSeq, bSeq)  
- positions (ai, bi)  
- thermodynamic stability (ΔG)  

An interaction is considered identical only if all attributes match exactly.

---

#### 3. Classification

The script outputs three sets:

- **Unique to dataset 1** — interactions only present in the first input  
- **Unique to dataset 2** — interactions only present in the second input  
- **Overlap** — interactions shared between both datasets  

---

## References

```
A.-C. Jousset, A. Hache, C. Jakob, B. Chane-Woon-Ming, D. Ferhadian, D. Desirò, R.P. Stansilaus, M. Marz, M. Schwemmle, H. Bolte, and R. Marquet.
"The wild-type packaging signal network cooperatively induces a metastable conformation of the influenza A virus genome during packaging."
In review, 2026.

R. Lorenz, S.H. Bernhart, C. Höner Zu Siederdissen, H. Tafer, C. Flamm C., P.F. Stadler, and I.L. Hofacker.
ViennaRNA Package 2.0.
Algorithms Mol. Biol. 2011.
```
