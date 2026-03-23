import scvi
import scanpy as sc
import numpy as np
import torch
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score

print("🚀 Phase 3: Initializing Best-Practice scANVI Aligner...")

# --- 1. DATA PREPARATION ---
adata = sc.read_h5ad("data/processed/harmonized_pancreas.h5ad")

# BEST PRACTICE: We treat Human labels as "Known" and Mouse labels as "Unknown"
# This allows the model to 'Transfer' the human cell-type identity to the mouse.
adata.obs["cell_type_scanvi"] = adata.obs["cell_type"].astype(str)
adata.obs.loc[adata.obs["species"] == "mouse", "cell_type_scanvi"] = "Unknown"

# --- 2. STEP 1: PRE-TRAIN WITH scVI (Unsupervised) ---
# This learns the basic background noise and species differences.
scvi.model.SCVI.setup_anndata(adata, batch_key="species", labels_key="cell_type_scanvi")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30)

print("🧠 Step 1: Pre-training scVI Background Model...")
vae.train(max_epochs=100)

# --- 3. STEP 2: REFINE WITH scANVI (Semi-Supervised) ---
# We initialize scANVI using the weights from the scVI model.
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    vae,
    labels_key="cell_type_scanvi",
    unlabeled_category="Unknown"
)

print("🧠 Step 2: Training scANVI Classifier (Label Transfer)...")
scanvi_model.train(max_epochs=20)

# --- 4. PREDICTION & CONFIDENCE ---
# Predict the Mouse labels based on the Human reference
adata.obs["predicted_label"] = scanvi_model.predict(adata)
# Get the probability (confidence) for each prediction
adata.obs["prediction_confidence"] = scanvi_model.predict(adata, soft=True).max(axis=1)

# Extract the latent space for visualization
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)

# --- 5. VALIDATION METRICS ---
print("📊 Calculating Final Validation...")
# Calculate Bio-Conservation (Higher is better)
bio_score = silhouette_score(adata.obsm["X_scANVI"], adata.obs["predicted_label"])
print(f"✅ Bio-Conservation (ASW-Label): {bio_score:.4f}")

# --- 6. VISUALIZATION ---
sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
sc.pl.umap(adata, color="predicted_label", ax=ax1, show=False, title="scANVI Predicted Labels")
sc.pl.umap(adata, color="prediction_confidence", ax=ax2, show=False, title="Prediction Confidence")
plt.savefig("results/scanvi_alignment_final.png")

# --- 7. EXPORT ---
adata.write("data/processed/aligned_pancreas_final.h5ad")
print("💾 Final scANVI-aligned object saved.")
