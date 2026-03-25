import scanpy as sc
import pandas as pd
import numpy as np

def audit_harmonized_data(file_path):
    print(f"🔍 Auditing: {file_path}")
    print("-" * 40)
    
    try:
        # 1. Load the AnnData object
        adata = sc.read_h5ad(file_path)
        
        # 2. Structural Summary
        print(f"✅ Dimensions: {adata.n_obs} cells x {adata.n_vars} genes")
        print(f"✅ Metadata Columns (obs): {list(adata.obs.columns)}")
        
        # 3. Data Peek (Top 5 Results)
        print("\n--- 🧫 CELL METADATA PEEK ---")
        print(adata.obs[['species', 'cell_type']].head())
        
        print("\n--- 🧬 GENE LIST PEEK ---")
        print(f"First 10 genes: {list(adata.var_names[:10])}")
        
        # 4. Expression Matrix Audit
        # Handling Sparse vs Dense matrices (Standard Bioinformatics Safety)
        sample_matrix = adata.X[:5, :5]
        if hasattr(sample_matrix, "toarray"):
            dense_sample = sample_matrix.toarray()
        else:
            dense_sample = sample_matrix
            
        print("\n--- 📉 COUNTS MATRIX (5x5 Sample) ---")
        print(dense_sample)
        
        # 5. Biological Integrity Check
        print("\n--- 🧪 BIOLOGICAL SUMMARY ---")
        print("Species Distribution:")
        print(adata.obs['species'].value_counts())
        
        print("\nCell Type Distribution:")
        print(adata.obs['cell_type'].value_counts())

        # Check for our "Gold Standard" research gene
        target = "SLC30A8"
        if target in adata.var_names:
            print(f"\n✅ Target Gene '{target}' found in dataset.")
        else:
            print(f"\n⚠️ Warning: Target Gene '{target}' MISSING from dataset.")

    except FileNotFoundError:
        print(f"❌ Error: File '{file_path}' not found. Run the simulation script first!")
    except Exception as e:
        print(f"❌ An unexpected error occurred: {e}")

if __name__ == "__main__":
    PATH = "data/processed/harmonized_pancreas.h5ad"
    audit_harmonized_data(PATH)
