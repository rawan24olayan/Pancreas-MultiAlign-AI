import requests
import os
from langgraph.graph import StateGraph, END
from typing import TypedDict, List

# --- CONFIGURATION ---
OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
GENE_FILE = "data/genes_to_test.txt"

class ValidatorState(TypedDict):
    genes: List[str]
    valid_count: int

# --- 1. THE VALIDATION LOGIC ---
def validate_genes_from_file():
    if not os.path.exists(GENE_FILE):
        print(f"❌ ERROR: {GENE_FILE} not found. Please create it first.")
        return

    with open(GENE_FILE, "r") as f:
        genes = [line.strip() for line in f if line.strip()]

    print(f"🧬 Validating {len(genes)} genes against Open Targets v4...")
    valid_genes = 0

    for symbol in genes:
        # Step A: Resolve Symbol to ID
        search_query = """
        query search($symbol: String!) {
          search(queryString: $symbol, entityNames: ["target"]) {
            hits { id }
          }
        }
        """
        try:
            r = requests.post(OPENTARGETS_URL, json={"query": search_query, "variables": {"symbol": symbol}}, timeout=5)
            hits = r.json().get("data", {}).get("search", {}).get("hits", [])
            
            if hits:
                print(f"  ✅ {symbol}: Valid (ID: {hits[0]['id']})")
                valid_genes += 1
            else:
                print(f"  ⚠️ {symbol}: Not found in Human database. (Check spelling/orthologs)")
        except Exception as e:
            print(f"  ❌ {symbol}: Connection error ({e})")

    print(f"\n📊 SUMMARY: {valid_genes}/{len(genes)} genes are ready for the Agent.")

# --- 2. THE LANGGRAPH SMOKE TEST (System Check) ---
def test_system_logic():
    print("\n🧠 Checking LangGraph Orchestration...")
    workflow = StateGraph(ValidatorState)
    workflow.add_node("check", lambda x: {"status": "System OK"})
    workflow.set_entry_point("check")
    workflow.add_edge("check", END)
    
    try:
        app = workflow.compile()
        print("✅ System Logic: PASS")
    except:
        print("❌ System Logic: FAIL")

if __name__ == "__main__":
    print("--- 🚀 STARTING CUSTOM DATA VALIDATION ---")
    validate_genes_from_file()
    test_system_logic()
