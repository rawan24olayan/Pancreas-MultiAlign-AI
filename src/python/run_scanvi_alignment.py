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
    Phase 2 Final: High-Resolution Cross-Species Alignment.
    Optimized for: High Label Separation (Biology) and Zero Batch Effect (Species).
    """
    print(f"🚀 Loading harmonized data from {input_path}...")
    adata = sc.read_h5ad(input_path)

    # --- Step 1: Optimized scVI Training (High Resolution) ---
    scvi.model.SCVI.setup_anndata(adata, batch_key="species")
    
    # SETTINGS: n_layers=1 to prevent over-smoothing; n_latent=50 for detail
    vae = scvi.model.SCVI(
        adata, 
        n_layers=1, 
        n_latent=50, 
        gene_likelihood="zinb",
        dispersion="gene-batch"
    )
    
    print("🧠 Training scVI (Optimizing Manifold Geometry)...")
    vae.train(max_epochs=400, early_stopping=True, batch_size=256)

    # --- Step 2: scANVI Fine-tuning (Focused Label Transfer) ---
    print("🧪 Training scANVI (Label Transfer Pass)...")
    lfm = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=adata,
        labels_key="cell_type",
        unlabeled_category="Unknown"
    )
    lfm.train(max_epochs=150)

    # --- Step 3: Extract Latent Space & Probabilistic Predictions ---
    adata.obsm["X_scANVI"] = lfm.get_latent_representation()
    adata.obs["predictions"] = lfm.predict()
    adata.obs["prediction_confidence"] = lfm.predict(soft=True).max(axis=1)

    # --- Step 4: Quality Control & Diagnostics ---
    print("📈 Generating Final QC Reports...")
    sc.pp.neighbors(adata, use_rep="X_scANVI")
    sc.tl.umap(adata)

    latent_space = adata.obsm["X_scANVI"]
    batch_asw = silhouette_score(latent_space, adata.obs["species"])
    
    # Only calculate Label ASW on labeled (Human) cells for benchmark
    human_mask = adata.obs["species"] == "human"
    label_asw = silhouette_score(latent_space[human_mask], adata.obs["cell_type"][human_mask])

    # Save Metrics to File for GitHub/CV Evidence
    with open("results/qc_metrics_summary.txt", "w") as f:
        f.write(f"BENGAL Phase 2 QC Report\n")
        f.write(f"-----------------------\n")
        f.write(f"Batch Mixing (Species ASW): {batch_asw:.4f}\n")
        f.write(f"Label Separation (Biology ASW): {label_asw:.4f}\n")
        f.write(f"Mean Prediction Confidence: {adata.obs['prediction_confidence'].mean():.4f}\n")

    # --- Step 5: Visualizations ---
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    sc.pl.umap(adata, color="species", title="Species Alignment (Mixing)", show=False, ax=ax[0])
    sc.pl.umap(adata, color="predictions", title="Cell Type Identity (Clustering)", show=False, ax=ax[1])
    plt.tight_layout()
    plt.savefig("results/qc_alignment_visual_final.png")
    plt.close()

    # --- Step 6: Save and Exit ---
    print(f"💾 Saving final aligned object to {output_path}...")
    adata.write_h5ad(output_path)
    lfm.save(model_path, overwrite=True)
    print("✅ Phase 2 Finalized.")

if __name__ == "__main__":
    os.makedirs("results", exist_ok=True)
    os.makedirs("models", exist_ok=True)
    
    run_bengal_alignment(
        input_path="data/processed/harmonized_pancreas.h5ad",
        output_path="data/processed/aligned_pancreas_final.h5ad"
    )
