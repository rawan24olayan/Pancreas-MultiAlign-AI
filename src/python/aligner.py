import scvi
import scanpy as sc

# 1. Load the data created by your R script
# We expect a 'Harmonized' file where genes are mapped to Human symbols
adata = sc.read_h5ad("data/processed/harmonized_pancreas.h5ad")

# 2. Setup the Model
# We tell the AI: "species" is a batch effect to ignore, 
# and "cell_timport scvi
import scanpy as sc

# 1. Load the data created by your R script
# We expect a 'Harmonized' file where genes are mapped to Human symbols
adata = sc.read_h5ad("data/processed/harmonized_pancreas.h5ad")

# 2. Setup the Model
# We tell the AI: "species" is a batch effect to ignore, 
# and "cell_type" is the biology to preserve.
scvi.model.SCVI.setup_anndata(adata, batch_key="species", labels_key="cell_type")

# 3. Initialize the VAE (Variational Autoencoder)
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)

# 4. Train the AI to "Align" the species
print("🧠 Training the Cross-Species Aligner...")
model.train(max_epochs=100)

# 5. Save the 'Aligned' Coordinates
adata.obsm["X_scVI"] = model.get_latent_representation()
adata.write("data/processed/aligned_pancreas.h5ad")
print("✅ Done! Human and Mouse cells are now aligned.")
