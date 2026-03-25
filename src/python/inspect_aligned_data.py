import scanpy as sc
import pandas as pd
import numpy as np

def audit_aligned_data(file_path):
    print(f"🔬 Auditing AI-Aligned Data: {file_path}")
    print("-" * 50)
    
    try:
        # 1. Load the AI-processed object
        adata = sc.read_h5ad(file_path)
        
        # 2. Check for the "AI Confidence" Score
        print("--- 🤖 AGENTIC METADATA (obs) ---")
        if "prediction_confidence" in adata.obs.columns:
            avg_conf = adata.obs["prediction_confidence"].mean()
            print(f"✅ Prediction Confidence Found. Average Score: {avg_conf:.4f}")
            print(adata.obs[['cell_type', 'prediction_confidence']].head())
        else:
            print("⚠️ Warning: 'prediction_confidence' column missing!")

        # 3. Check for the Latent Space (The "Hidden Features")
        print("\n--- 🧠 LATENT EMBEDDINGS (obsm) ---")
        if "X_scANVI" in adata.obsm.keys():
            latent_shape = adata.obsm["X_scANVI"].shape
            print(f"✅ X_scANVI Embedding Found. Shape: {latent_shape}")
            print(f"Sample Latent Vector (First 5 dims): {adata.obsm['X_scANVI'][0, :5]}")
        else:
            print("⚠️ Warning: No scANVI latent space found in .obsm!")

        # 4. Check for Differential Expression (Ranked Genes)
        print("\n--- 🧬 RANKED GENE GROUPS (uns) ---")
        if "rank_genes_groups" in adata.uns:
            print("✅ Gene Rankings Found.")
            # Show the top 5 genes for the 'Beta' cell type
            top_beta = adata.uns['rank_genes_groups']['names']['Beta'][:5]
            print(f"Top 5 Beta Cell Markers: {list(top_beta)}")
        else:
            print("ℹ️ Note: No gene rankings found yet. (The Research Agent will calculate these).")

    except FileNotFoundError:
        print(f"❌ Error: File '{file_path}' not found. Did you run the simulation or aligner?")

if __name__ == "__main__":
    PATH = "data/processed/aligned_pancreas_final.h5ad"
    audit_aligned_data(PATH)
