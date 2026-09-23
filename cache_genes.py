import json
import os
import re

from tauso.genome.read_human_genome import CANONICAL_CHROMS

# What an ASO can sensibly be designed against here. A mature miRNA is ~22 nt -- shorter than a
# gapmer with its flanks -- and snRNA, snoRNA, misc_RNA and the pseudogenes are not knockdown
# targets either, so offering them only invites a job that cannot mean anything. lncRNA stays:
# MALAT1 and NEAT1 are lncRNA and among the most commonly targeted transcripts there are.
TARGETABLE_BIOTYPES = {"protein_coding", "lncRNA"}


def gene_biotypes(gtf_path):
    """gene_name -> gene_type on the canonical chromosomes, read from the GTF's gene rows.

    Two genes can share a name -- DGCR5 is both an lncRNA and a transcribed pseudogene -- and the
    later row wins, which is the biotype the picker has always filtered on.
    """
    biotypes = {}
    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "gene" or fields[0] not in CANONICAL_CHROMS:
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
        genes = []
        gtf_path = os.path.join(db_dir, "GRCh38.gtf")
        if os.path.exists(gtf_path):
            biotypes = gene_biotypes(gtf_path)
            genes = sorted(name for name, biotype in biotypes.items() if biotype in TARGETABLE_BIOTYPES)
            print(f"Filtered {len(biotypes)} genes to {len(genes)} targetable ones.")
        else:
            print("WARNING: no GTF to read genes from.")

        # A shorter list is the intent; an empty one is an outage, because the picker then offers
        # nothing at all. An unreadable or unexpected GTF falls back to the annotation database:
        # it is the slow path this script exists to avoid, but it is always there.
        if not genes:
            print("WARNING: no targetable genes from the GTF; falling back to the annotation database.")
            from tauso.genome.read_human_genome import get_locus_to_data_dict

            genes = sorted(get_locus_to_data_dict(include_introns=False).keys())

        output_path = os.path.join(db_dir, "available_genes.json")
        with open(output_path, "w") as f:
            json.dump(genes, f)

        print(f"Successfully cached {len(genes)} genes to {output_path}")

    except Exception as e:
        print(f"CRITICAL ERROR caching genes: {e}")


if __name__ == "__main__":
    main()
