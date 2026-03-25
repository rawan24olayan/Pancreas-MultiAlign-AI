import pandas as pd
import numpy as np
import scanpy as sc
from anndata import AnnData
import os

# 1. Setup Directory Structure
os.makedirs("data/processed", exist_ok=True)
os.makedirs("results", exist_ok=True)

def generate_bengal_mock_data():
    print("🧪 Generating Simulated Cross-Species Data...")
    
    # Constants for a lightweight laptop test
    n_cells = 200
    n_genes = 500
    
    # Define Gene Names (Include key markers for the Agent to find)
    genes = [f"GENE_{i}" for i in range(n_genes)]
    genes[0:4] = ["INS", "SLC30A8", "GCK", "KCNJ11"] # Pancreatic Beta Cell Markers

    # --- PHASE 1: Harmonized Data (Simulates R Preprocessing Output) ---
    # Create random Poisson counts (mimics scRNA-seq sparsity)
    X1 = np.random.poisson(lam=0.5, size=(n_cells, n_genes))
    
    obs1 = pd.DataFrame({
        "species": np.random.choice(["human", "mouse"], n_cells),
        "cell_type": np.random.choice(["Alpha", "Beta", "Delta"], n_cells),
        "batch": np.random.choice(["Batch1", "Batch2", "Batch3"], n_cells)
    }, index=[f"cell_{i}" for i in range(n_cells)])
    
    adata_p1 = AnnData(X=X1.astype(np.float32), obs=obs1, var=pd.DataFrame(index=genes))
    
    # Save Phase 1 File
    adata_p1.write("data/processed/harmonized_pancreas.h5ad")
    print(f"✅ Phase 1 Created: 200 cells x 500 genes (harmonized_pancreas.h5ad)")

    # --- PHASE 2: Aligned Data (Simulates scANVI/BENGAL Output) ---
    # We copy Phase 1 and add the "AI-generated" columns
    adata_p2 = adata_p1.copy()
    
    # Simulate 'Prediction Confidence' (The Agent uses this to filter genes)
    # We make Beta cells have high confidence so the Agent picks them up
    conf_scores = np.random.uniform(0.5, 0.99, n_cells)
    adata_p2.obs["prediction_confidence"] = conf_scores
    
    # Simulate the Latent Space (X_scANVI) - 30 dimensions
    adata_p2.obsm["X_scANVI"] = np.random.normal(size=(n_cells, 30))
    
    # Save Phase 2 File
    adata_p2.write("data/processed/aligned_pancreas_final.h5ad")
    print(f"✅ Phase 2 Created: Includes latent space & confidence (aligned_pancreas_final.h5ad)")

if __name__ == "__main__":
    # Set seed for reproducibility
    np.random.seed(42)
    generate_bengal_mock_data()
    print("✨ Simulation Complete. Ready for run_pipeline.sh")
