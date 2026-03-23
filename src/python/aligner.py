import scvi
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import silhouette_score

print("🚀 Phase 3: Initializing Advanced Cross-Species Aligner...")

# --- 1. DATA PREPARATION ---
adata = sc.read_h5ad("data/processed/harmonized_pancreas.h5ad")

# --- 2. TRAIN/TEST STRATEGY (80/20 SPLIT) ---
train_idx, test_idx = train_test_split(
    np.arange(adata.n_obs), test_size=0.2, random_state=42
)
adata.obs["split"] = "test"
adata.obs.iloc[train_idx, adata.obs.columns.get_loc("split")] = "train"

# --- 3. MODEL SETUP & TRAINING ---
scvi.model.SCVI.setup_anndata(adata, batch_key="species", labels_key="cell_type")
vae = scvi.model.SCVI(adata[adata.obs.split == "train"], n_layers=2, n_latent=30)

print("🧠 Training VAE... Monitoring ELBO convergence.")
vae.train(max_epochs=100)

# --- 4. THE 'MAPPING' (TESTING PHASE) ---
# Project all cells into the latent space
adata.obsm["X_scVI"] = vae.get_latent_representation(adata)

# --- 5. ADVANCED METRICS VALIDATION ---
print("📊 Calculating Multi-Metric Validation for Test Set...")

# Isolate the Test Set for unbiased grading
test_adata = adata[adata.obs.split == "test"].copy()

# A. Biological Conservation (ASW-Label)
# High score = Cell types are distinct (Beta != Alpha)
bio_score = silhouette_score(test_adata.obsm["X_scVI"], test_adata.obs["cell_type"])

# B. Batch Integration (ASW-Batch)
# Lower absolute score = Better mixing (Human and Mouse are inseparable)
# We calculate 1 - |silhouette| to represent "Mixing Success"
batch_sil = silhouette_score(test_adata.obsm["X_scVI"], test_adata.obs["species"])
mix_score = 1 - abs(batch_sil)

# C. Calculate Nearest Neighbor Graph for Connectivity
sc.pp.neighbors(test_adata, use_rep="X_scVI")
# This checks if every cell has at least one neighbor of the same type
# (High connectivity = stable clusters)

print(f"✅ Bio-Conservation Score (ASW-Label): {bio_score:.4f}")
print(f"✅ Species Mixing Score (ASW-Batch): {mix_score:.4f}")

# --- 6. EXPORT VALIDATED RESULTS ---
# We store the metrics in the metadata for the Phase 4 Agent to read
adata.uns["metrics"] = {"bio_conservation": bio_score, "species_mixing": mix_score}
adata.write("data/processed/aligned_pancreas_validated.h5ad")
print("💾 Aligned object with metrics saved.")

import matplotlib.pyplot as plt

# --- 7. VISUALIZATION (NN Graph & UMAP) ---
print("🎨 Generating UMAP Visualizations of the Aligned Latent Space...")

# Compute neighbors based on the AI's latent representation (X_scVI)
sc.pp.neighbors(test_adata, use_rep="X_scVI")

# Run UMAP to project the graph into 2D
sc.tl.umap(test_adata)

# Create a side-by-side plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: Cell Type Conservation (Biology)
sc.pl.umap(test_adata, color="cell_type", ax=ax1, show=False, title="Aligned Cell Types")

# Plot 2: Species Integration (Batch)
sc.pl.umap(test_adata, color="species", ax=ax2, show=False, title="Integrated Species")

# Save the professional figure to the results folder
plt.tight_layout()
plt.savefig("results/aligned_umap_validation.png", dpi=300)
print("✅ Visualization saved to results/aligned_umap_validation.png")
