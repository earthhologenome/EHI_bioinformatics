## Raphael Eisenhofer 2023: I've rewritten this to work with Checkm2 output

import sys

if len(sys.argv) == 2:
    print("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tcoding_density\tsize")
else:
    print("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tcoding_density\tsize\tbinner")

for line in open(sys.argv[1]):
    if line.startswith("Name"):
        continue  # Skip the header line
    cols = line.strip().split("\t")
    name = cols[0]
    completeness = float(cols[1])
    contamination = float(cols[2])
    gc = float(cols[9])
    lineage = cols[11]
    n50 = int(cols[6])
    size = int(cols[8])
    coding_density = float(cols[5])


    if len(sys.argv) == 2:
        print("\t".join([name, str(completeness)[:5], str(contamination)[:5], str(gc)[:5], lineage, str(n50), str(coding_density)[:5], str(size)]))
    else:
        binner=sys.argv[2]
        print("\t".join([name, str(completeness)[:5], str(contamination)[:5], str(gc)[:5], lineage, str(n50), str(coding_density)[:5], str(size), binner]))
