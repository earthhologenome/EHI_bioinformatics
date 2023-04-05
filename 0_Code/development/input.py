input_file = "abb_input.tsv"

# Read in input file and create a list of unique combinations of ID and EHI_number
ids = set()
ehis = set()
with open(input_file) as f:
    for line in f:
        if line.startswith("ID"):
            continue
        id, batch, ehi = line.strip().split("\t")
        ids.add(id)
        ehis.add((id, ehi))
        
combinations = [(id, ehi) for id in ids for ehi in ehis if ehi[0] == id]

rule all:
    input:
        expand("output/{id}_{ehi}.txt", id=[c[0] for c in combinations], ehi=[c[1] for c in combinations])

rule run_pipeline:
    input:
        "input/{id}/{ehi}.fastq"
    output:
        "output/{id}_{ehi}.txt"
    shell:
        "touch --output {output}"