#!/bin/bash
set -e # Exit on any error

echo "🚀 BENGAL: Starting End-to-End Pipeline..."

# Step 1: Data Generation
echo "🧪 Phase 1: Generating Simulated Pancreas Data..."
python src/python/generate_simulated_data.py

# Step 2: AI Alignment
echo "🤖 Phase 2: Running scANVI Species Alignment..."
python src/python/run_scanvi_alignment.py

# Step 3: Agentic Discovery
echo "🔍 Phase 3: Launching LangGraph Research Agent..."
python src/python/research_agent.py

# Step 4: Statistical Reporting
echo "📊 Phase 4: Generating R Publication Report..."
Rscript src/R/NormConserved_Common_vs_Reference_Report.R

echo "✨ PIPELINE COMPLETE. Final findings in results/BENGAL_Discovery_Report.pdf"
