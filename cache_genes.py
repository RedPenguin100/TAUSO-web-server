import json
import os
import re

from tauso.genome.read_human_genome import get_locus_to_data_dict

# What an ASO can sensibly be designed against here. A mature miRNA is ~22 nt -- shorter than a
# gapmer with its flanks -- and snRNA, snoRNA, misc_RNA and the pseudogenes are not knockdown
# targets either, so offering them only invites a job that cannot mean anything. lncRNA stays:
# MALAT1 and NEAT1 are lncRNA and among the most commonly targeted transcripts there are.
TARGETABLE_BIOTYPES = {"protein_coding", "lncRNA"}


def gene_biotypes(gtf_path):
    """gene_name -> gene_type, read from the GTF's gene rows."""
    biotypes = {}
    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            name = re.search(r'gene_name "([^"]+)"', fields[8])
            biotype = re.search(r'gene_type "([^"]+)"', fields[8])
            if name and biotype:
                biotypes[name.group(1)] = biotype.group(1)
    return biotypes


def main():
    print("Starting gene caching process...")
    db_dir = os.environ.get("TAUSO_DATA_DIR", "/home/mambauser/.tauso_data")

    try:
        gene_to_data = get_locus_to_data_dict(include_introns=False)
        genes = sorted(gene_to_data.keys())

        gtf_path = os.path.join(db_dir, "GRCh38.gtf")
        if os.path.exists(gtf_path):
            biotypes = gene_biotypes(gtf_path)
            kept = [g for g in genes if biotypes.get(g) in TARGETABLE_BIOTYPES]
            # A shorter list is the intent; an empty one is an outage, so an unreadable or
            # unexpected GTF leaves the full list rather than emptying the picker.
            if kept:
                print(f"Filtered {len(genes)} genes to {len(kept)} targetable ones.")
                genes = kept
            else:
                print("WARNING: biotype filter matched nothing; keeping the full gene list.")
        else:
            print("WARNING: no GTF to read biotypes from; keeping the full gene list.")

        output_path = os.path.join(db_dir, "available_genes.json")
        with open(output_path, "w") as f:
            json.dump(genes, f)

        print(f"Successfully cached {len(genes)} genes to {output_path}")

    except Exception as e:
        print(f"CRITICAL ERROR caching genes: {e}")


if __name__ == "__main__":
    main()
