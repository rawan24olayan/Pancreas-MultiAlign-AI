import os
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import silhouette_score
import warnings

warnings.filterwarnings("ignore")

def run_bengal_alignment(input_path, output_path, model_path="models/scvi_model"):
    """
    Phase 2: Perform scANVI-based cross-species alignment and label transfer.
    Includes automated QC diagnostics for ELBO convergence and manifold mixing.
    """
    print(f"🚀 Loading harmonized data from {input_path}...")
    adata = sc.read_h5ad(input_path)

    # --- Step 1: scVI Training (Unsupervised Manifold Learning) ---
    # We define 'species' as the batch key to regress out technical noise
    scvi.model.SCVI.setup_anndata(adata, batch_key="species")
    
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="zinb")
    
    print("🧠 Training scVI model (Step A: Manifold Alignment)...")
    vae.train(max_epochs=250, early_stopping=True)

    # --- Step 2: Quality Control - ELBO Convergence ---
    print("📊 Generating Training Diagnostics...")
    plt.figure(figsize=(8, 5))
    plt.plot(vae.history["elbo_train"], label="Train")
    plt.plot(vae.history["elbo_validation"], label="Validation")
    plt.title("Model Convergence (ELBO)")
    plt.ylabel("-ELBO")
    plt.xlabel("Epoch")
    plt.legend()
    plt.savefig("results/qc_elbo_convergence.png")
    plt.close()

    # --- Step 3: scANVI Fine-tuning (Semi-Supervised Label Transfer) ---
    # We use the 'cell_type' column where Mouse is 'Unknown'
    print("🧪 Training scANVI model (Step B: Label Transfer)...")
    lfm = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=adata,
        labels_key="cell_type",
        unlabeled_category="Unknown"
    )
    lfm.train(max_epochs=100)

    # --- Step 4: Extract Predictions and Latent Space ---
    adata.obsm["X_scANVI"] = lfm.get_latent_representation()
    adata.obs["predictions"] = lfm.predict()
    # Extract Bayesian posterior probabilities (Confidence Scores)
    adata.obs["prediction_confidence"] = lfm.predict(soft=True).max(axis=1)

    # --- Step 5: Statistical QC - Mixing & Separation ---
    print("📈 Calculating Manifold Metrics...")
    sc.pp.neighbors(adata, use_rep="X_scANVI")
    sc.tl.umap(adata)

    # Calculate Silhouette Scores (ASW)
    # Batch ASW near 0 = Good mixing. Label ASW near 1 = Good biology.
    latent_space = adata.obsm["X_scANVI"]
    batch_asw = silhouette_score(latent_space, adata.obs["species"])
    
    # Label ASW check (on Human cells only to verify teacher quality)
    human_mask = adata.obs["species"] == "human"
    label_asw = silhouette_score(latent_space[human_mask], adata.obs["cell_type"][human_mask])

    print(f"\n--- QC REPORT ---")
    print(f"Batch Mixing (ASW): {batch_asw:.4f} (Target < 0.1)")
    print(f"Label Separation (ASW): {label_asw:.4f} (Target > 0.3)")
    
    # --- Step 6: Visual QC - UMAP Generation ---
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    sc.pl.umap(adata, color="species", title="Species Mixing (Alignment)", show=False, ax=ax[0])
    sc.pl.umap(adata, color="predictions", title="Predicted Cell Types", show=False, ax=ax[1])
    plt.tight_layout()
    plt.savefig("results/qc_alignment_visual.png")
    plt.close()

    # --- Step 7: Save Results ---
    print(f"💾 Saving aligned data to {output_path}...")
    adata.write_h5ad(output_path)
    
    # Optional: Save model for future use
    lfm.save(model_path, overwrite=True)
    print("✅ Phase 2 Complete.")

if __name__ == "__main__":
    # Ensure results directory exists
    os.makedirs("results", exist_ok=True)
    os.makedirs("models", exist_ok=True)
    
    run_bengal_alignment(
        input_path="data/processed/harmonized_pancreas.h5ad",
        output_path="data/processed/aligned_pancreas_final.h5ad"
    )
