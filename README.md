# SCRATCH-TumorHeterogeneity (Meta-Programs)

## Introduction
SCRATCH-TumorHeterogeneity identifies tumor meta-programs (MPs) from scRNA-seq datasets by:

* inputting annotated Seurat object 
* extracting and preprocessing per-sample matrices
* running NMF per sample across a rank grid
* aggregating programs across samples to derive reproducible meta-programs with consistent gene signatures

This module provides a clean three-stage workflow using QMD notebooks, fully orchestrated through Nextflow for scalability and portability using Docker or Singularity. It serves as both a standalone workflow and a core component of the broader SCRATCH single-cell analysis ecosystem.

---

## Prerequisites

* Nextflow ≥ 21.04.0  
* Java ≥ 8  
* Docker or Singularity/Apptainer  
* Git  
* R packages (automatically handled via container):  
  Seurat, Matrix, NMF, ggplot2, reshape2, viridis, pheatmap, data.table  

---

## Installation

```bash
git clone https://github.com/WangLab-ComputationalBiology/SCRATCH-TumorHeterogeneity.git
cd SCRATCH-TumorHeterogeneity
```

---

## Architecture

### Pipeline Stages (QMD notebooks)

| Stage | Notebook | Description |
|------:|----------|-------------|
| 1 | prep.qmd | Extract per-sample counts, apply log-CPM/10, filter, center, clip → `<sample>_preprocessed.rds` |
| 2 | nmf.qmd | Per-sample NMF on HVGs across rank grid → `<sample>_nmf_fit.rds` or `.SKIP.txt` |
| 3 | aggregate.qmd | Aggregate NMF fits into meta-programs and figures/tables |

### Orchestration (Nextflow components)

* main.nf — pipeline entrypoint  
* subworkflows/local/SCRATCH_MetaProg.nf — scatter/gather logic  
* modules/local/main.nf — QMD execution modules  
* nextflow.config — default runtime and container settings  

Parallelization is handled at the sample level to maximize HPC/cloud utilization while ensuring reproducibility.

---

## Quick Start

### Minimal example (Docker profile)

```bash
nextflow run main.nf -profile docker \
  --input_seurat_object /path/to/project_Azimuth_annotation_object.RDS \
  --project_name MyProject \
  --subset_col azimuth_labels \
  --subset_value Epithelial \
  -resume
```

---

## Typical workflow execution

1. prep: Runs once on full Seurat object  
2. nmf: Scattered execution per sample  
3. aggregate: Gathers all NMF fits into unified MPs  

---

## Key Parameters

### Shared Parameters

| Parameter | Description |
|----------|-------------|
| `--project_name` | Label for outputs |
| `--work_directory` | Output root (default: `./output`) |
| `--seed` | Reproducibility |

### Subsetting (prep stage)

| Parameter | Description |
|----------|-------------|
| `--subset_col`, `--subset_value` | Metadata-based selection (e.g. epithelial cells only) |

### Per-sample NMF (nmf stage)

| Parameter | Default | Purpose |
|----------|---------|---------|
| `--hvg_keep` | 5000 | Max HVGs retained |
| `--rank_lb`, `--rank_ub` | 3–7 | Rank search range |
| `--nrun` | 10 | NMF restarts |
| `--min_cells` | 100 | Skip small samples |

### Aggregation stage

| Parameter | Purpose |
|----------|---------|
| `--intra_min`, `--intra_max` | Within-sample filtering |
| `--inter_filter`, `--inter_min` | Cross-sample retention |
| `--min_intersect_initial`, `--min_intersect_cluster`, `--min_group_size` | MP clustering thresholds |

---

## Expected Input

A Seurat `.RDS` containing:
* multiple samples
* a metadata column allowing clean subsetting (exact string match required)

---

## Outputs

All outputs are stored in the `work_directory`:

```
data/per_sample_mat/
  <sample>_raw.rds
  <sample>_preprocessed.rds
  nmf_fit/<sample>_nmf_fit.rds
  nmf_fit/<sample>_nmf_fit.SKIP.txt

nmf_intersect.rds
nmf_programs_sig_filtered.rds
MP_table_50genes_per_MP.{rds,csv}
Cluster_list_final.rds
MP_list_final.rds

figures/metaprog/
  jaccard_heatmap_dendrogram.pdf
  NMF_cluster_pheatmap.pdf
```

These include:

* Meta-program signatures  
* Similarity heatmaps  
* Filtering artifacts and diagnostics  

---

## Example Full Run

```bash
nextflow run main.nf -profile singularity \
  --input_seurat_object project_Azimuth_annotation_object.RDS \
  --project_name Lung_MP \
  --subset_col azimuth_labels \
  --subset_value Epithelial \
  --hvg_keep 5000 --rank_lb 3 --rank_ub 7 --nrun 15 \
  --min_cells 150 \
  --intra_min 35 --intra_max 10 --inter_filter true --inter_min 10 \
  -resume
```

---

## Documentation
For more detailed documentation and advanced usage, refer to the Nextflow documentation and the comments within the subworkflow script (main.nf).

## Contributing
Contributions are welcome! Please submit a pull request or open an issue to discuss any changes.

## License
This project is available under the GNU General Public License v3.0. See the LICENSE file for more details.

## Contact
For questions or issues, please contact:

sazaidi@mdanderson.org

lwang22@mdanderson.org
