import scanpy as sc
import scvi
import os
import torch

def run_scanvi_alignment(input_path, output_path):
    """
    Performs scANVI-based label transfer from Human (Teacher) to Mouse (Student).
    Outputs Bayesian confidence scores for downstream agentic research.
    """
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"❌ Input data missing at {input_path}")

    print(f"🧬 Loading harmonized data from {input_path}...")
    adata = sc.read_h5ad(input_path)
    
    # 1. Setup scVI Model (Base Generative Model)
    # batch_key='species' aligns the Human and Mouse technical variation
    scvi.model.SCVI.setup_anndata(adata, batch_key="species", labels_key="cell_type")
    
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    
    print("🤖 Training scVI Base Model (Mapping Latent Space)...")
    vae.train(max_epochs=100, check_val_every_n_epoch=10)
    
    # 2. Setup scANVI (The Alignment & Label Transfer Model)
    # We use the Human labels to 'teach' the model what the 'Unknown' Mouse cells are
    print("🧠 Initializing scANVI for Cross-Species Label Transfer...")
    lvae = scvi.model.SCANVI.from_scvi_model(vae, unlabeled_category="Unknown")
    lvae.train(max_epochs=20)
    
    # 3. Extract Bayesian Results (The BENGAL Critical Step)
    print("📊 Calculating Bayesian Posterior Probabilities...")
    
    # FIX: Use predict(soft=True) for scvi-tools 1.x compatibility
    predictions_soft = lvae.predict(soft=True) 
    
    adata.obs["predictions"] = lvae.predict()
    adata.obs["prediction_confidence"] = predictions_soft.max(axis=1)
    adata.obsm["X_scANVI"] = lvae.get_latent_representation()
    
    # 4. Save the Final Aligned Object
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path)
    
    print(f"✅ Alignment Complete.")
    print(f"💾 File Saved: {output_path}")
    print(f"📈 Mean Confidence: {adata.obs['prediction_confidence'].mean():.4f}")

if __name__ == "__main__":
    # Ensure paths match your run_pipeline.sh
    IN_FILE = "data/processed/harmonized_pancreas.h5ad"
    OUT_FILE = "data/processed/aligned_pancreas_final.h5ad"
    run_scanvi_alignment(IN_FILE, OUT
