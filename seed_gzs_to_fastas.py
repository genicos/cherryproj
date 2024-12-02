import os, sys, gzip
from pathlib import Path

usage = f"usage: {sys.argv[0]} <seed gz dir> <output dir>"


seed_dir = sys.argv[1]
output_dir = sys.argv[2]
seed_files = []
for root,dirs,files in os.walk(seed_dir):
    for name in files:
        if name.endswith(('gz')): # if fasta already generated
            seed_files.append(os.path.join(root,name))

Path(output_dir).mkdir(parents=True, exist_ok=True)
for seed_file in seed_files:
    fasta_name = str(Path(os.path.basename(seed_file)).stem.split(".")[0]) + ".fasta"
    with open(os.path.join(output_dir,fasta_name), 'w') as F:
        with gzip.open(seed_file ,'rt') as S:
            for line in S:
                if not line.startswith("#"):
                    if len(line.strip().split()) == 2:
                        id, aln = line.strip().split()
                        aln = aln.replace('.','-')
                        F.write(f">{id}\n")
                        F.write(f"{aln}\n")
