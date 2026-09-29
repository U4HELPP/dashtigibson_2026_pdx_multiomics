# Pan-cancer comparisons of multi-omics data for therapeutic selection — Analysis Code Repository

This repository consolidates the analysis code referenced in the manuscript below. Code is either included directly in this repo or linked to its original source repository, depending on where it was originally developed and published.

## Citation

Dashti-Gibson N, Olex AL, Tyc KM, Barber AM, Huynh KLA, Smith OM, Hu B, Bhoopathi P, Bala MMI, Walker CJ, Yajaman DR, Myrick RK, Sallee AJ, Murugavel N, Wang L, Puchalapalli M, Gibson AE, Thompson GR, Vudatha V, Bandy VL, Punjala A, Bear HD, Robila V, McGuire KP, Shah RD, Fernandez LJ, Lee SD, Reyna J, Liu J, Koblinski JE, Dozmorov MG, Woods NT, Trevino JG, Harrell JC. **Pan-cancer comparisons of multi-omics data for therapeutic selection.** *Genome Medicine*. [In submission — volume/issue/DOI to be added upon publication.]

## License Note

Code included directly in this repository is covered by this repository's license (see `LICENSE`). Repositories linked from this README (marked **Linked** or **Private** in the table below) are separate projects maintained elsewhere and are governed by their own licenses. This repository's license does **not** extend to or supersede the licensing terms of any linked external repository.

## Status Legend
- Included — Included in this repo
- Linked — External repo (linked, not copied)
- Private — External repo, currently private, needs to be made public
- Needs organizing — Needs organizing/cleanup before inclusion
- Pending — To be collected from contributor

---

## 1. WGS Data Analysis and Ancestry
- **Status:** Pending
- **Location:** TBD
- **Contributor(s):** Katarzyna Tyc
- **Notes:** Will share repo once it is polished.

## 2. Bulk RNA-Seq
- **Status:** Linked (JAX pipeline) / Private (SOP) / Pending (PCA/PLS-DA scripts)
- **Location:**
  - Pipeline: https://github.com/TheJacksonLaboratory/jds-nf-workflows
  - Processing SOP: https://github.com/U4HELPP/SOP_bulk-rnaseq-processing *(currently private — needs to be made public)*
  - PCA plots / PLS-DA analysis script: TBD
- **Contributor(s):** Mikail Bala, Amy Olex and Nina Dashti-Gibson
- **Notes:** Nina will provide the PCA/PLS-DA scripts. SOP repo needs visibility changed to public before publication.

## 3. Proteomics Analysis
- **Status:** Pending
- **Location:** TBD
- **Contributor(s):** TBD — contributor to confirm
- **Notes:**

## 4. Gene Set Variation Analysis (GSVA) and Figures
- **Status:** Pending
- **Location:** TBD
- **Contributor(s):** Nina Dashti-Gibson
- **Notes:** Nina will provide.

## 5. Molecular Subtyping
- **Status:** Private (Breast) / Pending (PDAC)
- **Location:**
  - Breast (PAM50): https://github.com/U4HELPP/SOP_breast-cancer-PAM50-classification
  - PDAC (PurIST): https://github.com/naimurashid/PurIST
  - CRC (CMS): https://github.com/peterawe/CMScaller 

- **Contributor(s):** Amy Olex, Nina Dashti-Gibson
- **Notes:** Need PDAC and CRC script that uses these packages. Breast PAM50 repo needs visibility changed to public before publication.

## 6. High-Throughput Drug Screening Analysis
- **Status:** Pending
- **Location:** TBD
- **Contributor(s):** Nina Dashti-Gibson
- **Notes:** Nina will provide.

## 7. scRNA-Seq Analysis
- **Status:** Needs organizing (alignment/clustering) / Included-pending-copy (cell filtering, merging, pseudobulk)
- **Location:**
  - Alignment, cell-level filtering, sample merging, ambient RNA adjustment and Pseudobulk PAM50 subtyping are described here: https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/scRNASeq
- **Contributor(s):** Amy Olex
- **Notes:** Organization completed

## 8. Drug Response Prediction Modeling
- **Status:** Pending
- **Location:** TBD
- **Contributor(s):** TBD — contributor to confirm
- **Notes:**
