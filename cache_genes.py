import json
import os
from tauso.genome.read_human_genome import get_locus_to_data_dict

def main():
    print("Starting gene caching process...")
    db_dir = os.environ.get("TAUSO_DATA_DIR", "/home/mambauser/.tauso_data")

    try:
        # Fetch the genes using the correct TAUSO import
        gene_to_data = get_locus_to_data_dict(include_introns=False)
        genes = sorted(list(gene_to_data.keys()))

        # Save to JSON
        output_path = os.path.join(db_dir, "available_genes.json")
        with open(output_path, "w") as f:
            json.dump(genes, f)

        print(f"Successfully cached {len(genes)} genes to {output_path}")

    except Exception as e:
        print(f"CRITICAL ERROR caching genes: {e}")

if __name__ == "__main__":
    main()