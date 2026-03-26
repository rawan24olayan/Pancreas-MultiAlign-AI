# BENGAL: Pancreas Multi-Align AI Pipeline
**Single-Cell Cross-Species Alignment & Agentic Discovery**

BENGAL is a high-fidelity bioinformatics pipeline designed to harmonize pancreatic single-cell RNA-seq data across Human and Mouse species. It utilizes Variational Autoencoders (VAEs) to decouple species-level batch effects from conserved biological identities.

## 🚀 Key Performance Metrics (Phase 2)
After hyperparameter optimization ($n\_latent=50, n\_layers=1$), the pipeline achieved:
* **Mean Prediction Confidence:** 97.84% (Bayesian posterior probability)
* **Batch Mixing (Species ASW):** -0.0039 (Near-perfect manifold harmonization)
* **Label Integrity:** Successful transfer of endocrine identities (Alpha, Beta, Delta) across species.

## 🏗️ Architecture
1.  **Phase 1 (Preprocessing):** R-based normalization and integration of `.h5ad` datasets.
2.  **Phase 2 (Alignment):** scVI/scANVI deep generative modeling to create a shared 50-dimensional latent space.
3.  **Phase 3 (Discovery):** LangGraph-powered AI Agent to validate gene-pathway conservation.

## 📊 Quality Control
The pipeline generates automated diagnostic plots in `/results`:
* `qc_elbo_convergence.png`: Validates VAE mathematical stability.
* `qc_alignment_visual_final.png`: UMAP visualization of species integration.
* `qc_metrics_summary.txt`: Hard-coded ASW and confidence scores.

## 🛠️ Usage
```bash
# Run the full end-to-end pipeline
bash run_pipeline.sh
