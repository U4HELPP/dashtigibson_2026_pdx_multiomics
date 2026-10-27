# Cross-cancer comparisons of multi-omics data for therapeutic selection — Analysis Code Repository

This repository consolidates the analysis code referenced in the manuscript below. Code is either included directly in this repo or linked to its original source repository, depending on where it was originally developed and published.

## Citation

Dashti-Gibson N, Olex AL, Tyc KM, Barber AM, Huynh KLA, Smith OM, Hu B, Bhoopathi P, Bala MMI, Walker CJ, Yajaman DR, Myrick RK, Sallee AJ, Murugavel N, Wang L, Puchalapalli M, Gibson AE, Thompson GR, Vudatha V, Bandy VL, Punjala A, Bear HD, Robila V, McGuire KP, Shah RD, Fernandez LJ, Lee SD, Reyna J, Liu J, Koblinski JE, Dozmorov MG, Woods NT, Trevino JG, Harrell JC. **Cross-cancer comparisons of multi-omics data for therapeutic selection.** *Genome Medicine*. [In submission — volume/issue/DOI to be added upon publication.]

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
- **Status:** Private
- **Location:** https://github.com/TYCGenomicsLab/publications/tree/main/2026_DashtiGibson_pancancer
- **Contributor(s):** Katarzyna M. Tyc
- **Notes:** The linked repository is currently private and will be made public upon publication.

## 2. Bulk RNA-Seq Processing
- **Status:** Included
- **Location:** https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/bulk-rnaseq
- **Contributor(s):** Mikail Bala and Amy Olex
- **Notes:**

## 3. Proteomics Processing
- **Status:** Linked
- **Location:** https://github.com/Alexa-Barber/U54_Proteomic_Workflows
- **Version/Release:** [v1.0](https://github.com/Alexa-Barber/U54_Proteomic_Workflows/releases/tag/v1.0)
- **Contributor(s):** Alexa Barber
- **Notes:** Pipeline and code documentation is located in the above repo's README.

## 4. Comparative Omics Analysis
- **Status:** Included
- **Location:** https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/comparative_omics
- **Contributor(s):** Nina Dashti-Gibson
- **Notes:** 

## 5. Molecular Subtyping
- **Status:** Private (Breast) / Included (PDAC and CRC)
- **Location:**
  - Breast (PAM50): Refer to the [Bulk-RNASeq README](https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/bulk-rnaseq). Subtyping script is located at https://github.com/AmyOlex/WCCTR_RNASeq_Pipeline/blob/master/WorkstationScripts/06_PAM50_subtyping.R
  - PDAC (PurIST) and CRC (CMS): https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/bulk-rnaseq/analysis
- **Contributor(s):** Amy Olex and Nina Dashti-Gibson
- **Notes:** Breast PAM50 repo needs visibility changed to public before publication.

## 6. High-Throughput Drug Screening Analysis
- **Status:** Pending
- **Location:** https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/drug_screening
- **Contributor(s):** Nina Dashti-Gibson
- **Notes:** 

## 7. scRNA-Seq Analysis
- **Status:** Included / Linked
- **Location:**
  - Alignment, cell-level filtering, sample merging, ambient RNA adjustment and Pseudobulk PAM50 subtyping are described here: https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/scRNASeq
  - The primary code repo is here: https://github.com/AmyOlex/WCCTR_RNASeq_Pipeline/tree/v1.0.0/SingleCell
- **Contributor(s):** Amy Olex
- **Notes:** Organization completed.

## 8. Drug Response Prediction Modeling
- **Status:** Linked
- **Location:**
    - Code, model run, and results are described here: https://github.com/huynhkl953/2026_pdx_multiomics_drug_model
- **Contributor(s):** Khoa Huynh
- **Notes:**
