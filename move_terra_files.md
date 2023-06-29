```
#!/bin/bash
#SBATCH --job-name=terra_transfer
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=main
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=16
#SBATCH --output=%x.%j.log
#SBATCH --time=3:00:00

conda activate terra-env

python3 ~/progs/terra_scripts/pull_terra_table.py \
  --workspace hprc_polishing \
  --workspace-namespace human-pangenome-ucsc \
  --table-name phasing_homozygous \
  --dir /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/terra_phasing_homozygous \
  --threads 1
```
