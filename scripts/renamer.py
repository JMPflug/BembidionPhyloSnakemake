# Version: 31 Jan 2021 v2
#   Can exclude loci based on taxon occupancy
#       -specified in config.yaml

from Bio import SeqIO

def load_names(namefile):
    names = {}
    with open(namefile) as f:
        for line in f:
            split = line.split()
            # print("{} {}".format(split[0], split[3]))
            names[split[0]] = split[3]
    return names


try:
    excluded = [f for f in snakemake.config["Excluded_taxa"].split()]
except:
    excluded = [None]
print(excluded)


def rename(input, names):
    with open(input) as f:
        records = []
        for record in SeqIO.parse(f, "fasta"):
            record.description = ""
            if record.id not in excluded and names[record.id] not in excluded:
                if record.id in names:
                    # print("Replacing {} with {}".format(record.id, names[record.id]))
                    record.id = record.id.replace(record.id, names[record.id])
                    records.append(record)
                else:
                    records.append(record)
            else:
                print("Excluded taxon: {} or {}".format(record.id, names[record.id]))

    SeqIO.write(records, snakemake.output[0], "fasta")

# print("RENAMER INPUT:")
# print(snakemake.input[0])
rename(snakemake.input[0], load_names("resources/OldToNewNames.txt"))
