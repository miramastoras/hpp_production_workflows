


Debug taylors script:
```
mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --logDebug --jobStore ./big_store --batchSystem slurm --batchLogsDir ./logs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/workflows/dipcall_happy_eval.wdl /private/groups/patenlab/tswon/100523_downsampling/dipcall_happy_output/dipcall_happy_downsample_30x.json -o dipcall_output_35x -m toil_dipcall_35x.json 2>&1 | tee toil_dipcall_35x_log.txt
```


```
bcftools isec /private/groups/patenlab/tswon/100523_downsampling/deepPolisher_output/DP_output_35x/polisher_output.vcf.gz /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG002_y2_DCv1.2_PHv5_DPmodel5/HG002_y2_DCv1.2_PHv5_DPmodel5_polisher_output.vcf.gz -p /private/groups/patenlab/mira/hprc_polishing/debug
```
Run dipcall with polished assemblies

```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/dipcall.wdl
```

```
{
  "runDipcall.dipcall.referenceFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runDipcall.dipcall.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runDipcall.dipcall.threadCount": 64,
  "runDipcall.dipcall.assemblyFastaMat": "/private/groups/patenlab/tswon/100523_downsampling/applyPolish_output/AP_output_hap2_35x/HG002_y2_PHv5_DeepPolisher_model5_35x_hap2.polished.fasta",
  "runDipcall.dipcall.isMaleSample": true,
  "runDipcall.dipcall.memSizeGB": 64,
  "runDipcall.dipcall.assemblyFastaPat": "/private/groups/patenlab/tswon/100523_downsampling/applyPolish_output/AP_output_hap1_35x/HG002_y2_PHv5_DeepPolisher_model5_35x_hap1.polished.fasta",
  "runDipcall.dipcall.referenceIsHS38": true
}
```


```
mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --jobStore ./big_store --batchSystem slurm --batchLogsDir ./logs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/dipcall.wdl dipcall_35x_wdl_inputs.json -o toil_dipcall_35x_wdl_out -m toil_dipcall_35x_wdl.json 2>&1 | tee dipcall_35x_log.txt
```

Bedtools intersect
```
bedtools intersect -a /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -b HG002_y2_PHv5_DeepPolisher_model5_35x_hap1.polished.dipcall.bed > 35x_dipcall_GIABconf.bed
```

Hap.py
```
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz /private/groups/patenlab/mira/hprc_polishing/debug/toil_dipcall_35x_wdl_out/HG002_y2_PHv5_DeepPolisher_model5_35x_hap1.polished.dipcall.vcf.gz -r /private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta -f /private/groups/patenlab/mira/hprc_polishing/debug/toil_dipcall_35x_wdl_out/35x_dipcall_GIABconf.bed -o /private/groups/patenlab/mira/hprc_polishing/debug/toil_dipcall_35x_wdl_out/35x_happy --pass-only --no-roc --no-json --engine=vcfeval --threads=32

```

Rerun whole pipeline, using taylors json
```
/private/groups/patenlab/tswon/100523_downsampling/dipcall_happy_output/dipcall_happy_downsample_35x.json

mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --jobStore ./big_store --batchSystem slurm --batchLogsDir ./logs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/workflows/dipcall_happy_eval.wdl dipcall_happy_downsample_35x_taylor.json -o toil_dipcall_happy_downsample_35x_taylor_out -m toil_dipcall_happy_downsample_35x_taylor.json 2>&1 | tee toil_dipcall_happy_downsample_35x_taylor_log.txt

```


```
{
  "dipcall_happy.dipcall_t.threadCount": 128,
  "dipcall_happy.happy_t.memSizeGB": 256,
  "dipcall_happy.happy_t.sample": "HG002",
  "dipcall_happy.dipcall_t.assemblyFastaMat": "/private/groups/patenlab/mira/hprc_polishing/wdl_tests/full_hprc_WDL_mixed_aligners/toil_hprc_DP_wdl_out/HG002.y2.polished.full_hprc_WDL_mixed_aligners.hap2.fasta",
  "dipcall_happy.happy_t.threadCount": 128,
  "dipcall_happy.happy_t.truthVCF": "/private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
  "dipcall_happy.dipcall_t.assemblyFastaPat": "/private/groups/patenlab/mira/hprc_polishing/wdl_tests/full_hprc_WDL_mixed_aligners/toil_hprc_DP_wdl_out/HG002.y2.polished.full_hprc_WDL_mixed_aligners.hap1.fasta",
  "dipcall_happy.confidenceBedFile": "/private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
  "dipcall_happy.bedtoolsIntersect.threadCount": 32,
  "dipcall_happy.referenceFasta": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "dipcall_happy.dipcall_t.memSizeGB": 256,
  "dipcall_happy.referenceFastaFai": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "dipcall_happy.dipcall_t.isMaleSample": true,
  "dipcall_happy.bedtoolsIntersect.memSizeGB": 16
}
```
This failed... `toil_dipcall_happy_downsample_35x_taylor_log.txt`

Rerun Happy with exact same input data as the manual run, but with the happy wdl

```
docker run -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz /private/groups/patenlab/mira/hprc_polishing/debug/toil_dipcall_35x_wdl_out/HG002_y2_PHv5_DeepPolisher_model5_35x_hap1.polished.dipcall.vcf.gz -r /private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta -f /private/groups/patenlab/mira/hprc_polishing/debug/toil_dipcall_35x_wdl_out/35x_dipcall_GIABconf.bed -o /private/groups/patenlab/mira/hprc_polishing/debug/toil_dipcall_35x_wdl_out/35x_happy --pass-only --no-roc --no-json --engine=vcfeval --threads=32

docker run --rm -it -u `id -u`:`id -g` -v /private/groups/patenlab/mira:/private/groups/patenlab/mira jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py /private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz /private/groups/patenlab/mira/hprc_polishing/debug/toil_dipcall_35x_wdl_out/HG002_y2_PHv5_DeepPolisher_model5_35x_hap1.polished.dipcall.vcf.gz -r /private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta -f /private/groups/patenlab/mira/hprc_polishing/debug/toil_dipcall_35x_wdl_out/35x_dipcall_GIABconf.bed -o /private/groups/patenlab/mira/hprc_polishing/debug/toil_dipcall_35x_wdl_out/35x_happy --pass-only --engine=vcfeval --threads=32
```

```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/hapDotPy.wdl
```

```
{
  "runhapDotPy.hapDotPy.assemblyIndex": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai",
  "runhapDotPy.hapDotPy.bedRegions": "/private/groups/patenlab/mira/hprc_polishing/debug/toil_dipcall_35x_wdl_out/35x_dipcall_GIABconf.bed",
  "runhapDotPy.hapDotPy.truthVCF": "/private/groups/patenlab/mira/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
  "runhapDotPy.hapDotPy.queryVCF": "/private/groups/patenlab/mira/hprc_polishing/debug/toil_dipcall_35x_wdl_out/HG002_y2_PHv5_DeepPolisher_model5_35x_hap1.polished.dipcall.vcf.gz",
  "runhapDotPy.hapDotPy.threadCount": 32,
  "runhapDotPy.hapDotPy.assembly": "/private/groups/patenlab/mira/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
  "runhapDotPy.hapDotPy.sample": "HG002"
}
```

```
mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --jobStore ./big_store --batchSystem slurm --batchLogsDir ./logs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/hapDotPy.wdl happy_wdl_docker.json -o toil_happy_wdl_docker_out -m toil_happy_wdl_docker.json 2>&1 | tee toil_happy_wdl_docker_log.txt

mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --jobStore ./big_store2 --batchSystem slurm --batchLogsDir ./logs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/hapDotPy.wdl happy_new_docker.json -o toil_happy_new_docker_out -m toil_happy_new_docker.json 2>&1 | tee toil_happy_new_docker_log.txt

```

Rerun winnowmap alignments with more coverage files

HG002

```
{
  "longReadAlignmentScattered.preset": "map-pb",
  "longReadAlignmentScattered.sampleSuffix": "DCv1.2.full_coverage.winnowmapv2.03",
  "longReadAlignmentScattered.aligner": "winnowmap",
  "longReadAlignmentScattered.kmerSize": 15,
  "longReadAlignmentScattered.readFiles": ["/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/m64011_190714_120746.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/m64011_190728_111204.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/m64011_190830_220126.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/m64011_190901_095311.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/m64012_190920_173625.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/m64012_190921_234837.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/m64015_190920_185703.dc.q20.fastq.gz","/private/groups/patenlab/masri/hprc/polishing/HG002/reads/HiFi_DC_v1.2/m64015_190922_010918.dc.q20.fastq.gz"],
  "longReadAlignmentScattered.dockerImage": "mobinasri/long_read_aligner:v0.3.3",
  "longReadAlignmentScattered.alignment.suffix": "all2dip.winnowmap",
  "longReadAlignmentScattered.sampleName": "HG002",
  "longReadAlignmentScattered.assembly": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "longReadAlignmentScattered.options": "--cs --eqx -L -Y -I8g"
}
```
```
mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --jobStore ./big_store --batchSystem slurm --batchLogsDir ./logs ~/progs/hpp_production_workflows/QC/wdl/tasks/long_read_aligner_scattered_PhaseHom.wdl long_read_alignment_scattered.json -o toil_winnow_full_cov_out -m toil_winnow_full_cov.json 2>&1 | tee toil_winnow_full_cov_log.txt

```

Run PHARAOH on HG002 80x hifi  alignments

```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/workflows/PHARAOH.wdl
```


```
{
  "PHARAOH.maternalFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa",
  "PHARAOH.allONTToMatBam": "/private/groups/patenlab/masri/hprc/polishing/HG002/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/UL_alignments/mat/slurm_run/HG002.mat.trio_hifiasm_0.19.5.DC_1.2.UL_R941_Guppy6_40x.winnowmap_2.0.3.bam",
  "PHARAOH.allONTToMatBai": "/private/groups/patenlab/masri/hprc/polishing/HG002/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/UL_alignments/mat/slurm_run/HG002.mat.trio_hifiasm_0.19.5.DC_1.2.UL_R941_Guppy6_40x.winnowmap_2.0.3.bam.bai",
  "PHARAOH.diploidFaGz": "/private/groups/patenlab/masri/hprc/polishing/HG002/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/assembly/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.gz",
  "PHARAOH.paternalFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa",
  "PHARAOH.allHifiToDiploidBai": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/HG002.DCv1.2.full_coverage.winnowmapv2.03.bam.bai",
  "PHARAOH.maternalFastaIndex": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa.fai",
  "PHARAOH.allHifiToDiploidBam": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/alignments/hifi_winnowmap/toil_winnow_full_cov_out/HG002.DCv1.2.full_coverage.winnowmapv2.03.bam",
  "PHARAOH.allONTToPatBai": "/private/groups/patenlab/masri/hprc/polishing/HG002/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/UL_alignments/pat/slurm_run/HG002.pat.trio_hifiasm_0.19.5.DC_1.2.UL_R941_Guppy6_40x.winnowmap_2.0.3.bam.bai",
  "PHARAOH.allONTToPatBam": "/private/groups/patenlab/masri/hprc/polishing/HG002/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/UL_alignments/pat/slurm_run/HG002.pat.trio_hifiasm_0.19.5.DC_1.2.UL_R941_Guppy6_40x.winnowmap_2.0.3.bam",
  "PHARAOH.sampleName": "HG002",
  "PHARAOH.paternalFastaIndex": "/private/groups/patenlab/mira/hprc_polishing/data/HG002_y2_polishing/HG002.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa.fai",
  "PHARAOH.useMargin": false,
  "PHARAOH.PharaohAligner":"winnowmap",
  "PHARAOH.PharaohKmerSize":"15",
  "PHARAOH.PharaohHiFiPreset":"map-pb",
  "PHARAOH.pafAligner":"minimap2"
}
```

```
cd /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG002_y2_DCv1.2_80x_PHv5

mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --jobStore ./big_store --batchSystem slurm --batchLogsDir ./logs ~/progs/hpp_production_workflows/QC/wdl/workflows/PHARAOH.wdl PHARAOH.json -o toil_pharaoh_out -m toil_pharaoh.json 2>&1 | tee toil_pharaoh_log.txt
```

HG005

```
{
  "longReadAlignmentScattered.preset": "map-pb",
  "longReadAlignmentScattered.sampleSuffix": "DCv1.2.full_coverage.winnowmapv2.03",
  "longReadAlignmentScattered.aligner": "winnowmap",
  "longReadAlignmentScattered.kmerSize": 15,
  "longReadAlignmentScattered.readFiles": ["https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191118_150849.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191120_193948.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191122_184406.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191124_055423.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191126_155613.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191127_220906.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191129_043425.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191202_204405.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191204_164321.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191205_225630.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191207_052215.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191209_211903.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191211_182504.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191213_003759.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191214_070352.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191216_194629.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191218_164535.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191219_225837.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_191221_052416.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200107_170917.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200108_232219.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200112_090459.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200723_190224.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200730_190124.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200801_011415.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64017_200802_073944.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64109_200210_210230.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64109_200304_195708.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64109_200309_192110.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64109_200311_013444.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64109_200805_204709.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64109_200807_075817.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64109_200808_191025.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64109_200810_062248.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64109_200813_162416.dc.q20.fastq.gz","https://storage.googleapis.com/brain-genomics/kishwar/share/deepconsensus/v1.2/m64109_200815_033514.dc.q20.fastq.gz"],
  "longReadAlignmentScattered.dockerImage": "mobinasri/long_read_aligner:v0.3.3",
  "longReadAlignmentScattered.alignment.suffix": "all2dip.winnowmap",
  "longReadAlignmentScattered.sampleName": "HG005",
  "longReadAlignmentScattered.assembly": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "longReadAlignmentScattered.options": "--cs --eqx -L -Y -I8g"
}
```
```
mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --jobStore ./big_store --batchSystem slurm --batchLogsDir ./logs ~/progs/hpp_production_workflows/QC/wdl/tasks/long_read_aligner_scattered_PhaseHom.wdl long_read_aligner_scattered.json -o toil_winnow_full_cov_out -m toil_winnow_full_cov.json 2>&1 | tee toil_winnow_full_cov_log.txt
```

Get coverage of HG005 bamfile
```
#!/bin/bash
#SBATCH --job-name=coverage
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=main
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=16
#SBATCH --output=%x.%j.log
#SBATCH --time=5:00:00

samtools coverage /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_winnowmap/toil_winnow_full_cov_out/HG005.DCv1.2.full_coverage.winnowmapv2.03.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_winnowmap/samtools_coverage.txt
```

```
# Get averages of sequencing depths
# multiple coverage for each contig by number of bases

awk '{ print $3 * $7 }' samtools_coverage.txt | awk '{sum+=$1} END { print sum}'   
399689567823

divide by number of bases
awk '{sum+=$3} END {print sum}' samtools_coverage.txt
399689567823 / 2930142973 = 136x

60/130 = 0.4615

0.45
```

```
#!/bin/bash
#SBATCH --job-name=downsampling_60x
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=main
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --time=3:00:00

samtools view -s 0.45 -b -h -@ 64 /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_winnowmap/toil_winnow_full_cov_out/HG005.DCv1.2.full_coverage.winnowmapv2.03.bam > /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_winnowmap/downsample_60x/HG005.DCv1.2.60x.winnowmapv2.03.bam

samtools index /private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_winnowmap/downsample_60x/HG005.DCv1.2.60x.winnowmapv2.03.bam
```
Run PHARAOH with 60x coverage bamfile
```
java -jar /private/home/mmastora/progs/womtool-85.jar inputs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/workflows/PHARAOH.wdl
```

```
{
  "PHARAOH.maternalFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa",
  "PHARAOH.PharaohAligner": "winnowmap",
  "PHARAOH.PharaohHiFiPreset": "map-pb",
  "PHARAOH.allONTToMatBam": "/private/groups/patenlab/masri/hprc/polishing/HG005/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/UL_alignments/mat/slurm_run/HG005.mat.trio_hifiasm_0.19.5.DC_1.2.winnowmap_2.03.UL_R941_Guppy5_40x.bam",
  "PHARAOH.diploidFaGz": "/private/groups/patenlab/masri/hprc/polishing/HG005/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/assembly/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa.gz",
  "PHARAOH.pafAligner": "winnowmap",
  "PHARAOH.PharaohKmerSize": "15",
  "PHARAOH.allONTToMatBai": "/private/groups/patenlab/masri/hprc/polishing/HG005/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/UL_alignments/mat/slurm_run/HG005.mat.trio_hifiasm_0.19.5.DC_1.2.winnowmap_2.03.UL_R941_Guppy5_40x.bam.bai",
  "PHARAOH.paternalFasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa",
  "PHARAOH.allHifiToDiploidBai": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_winnowmap/downsample_60x/HG005.DCv1.2.60x.winnowmapv2.03.bam.bai",
  "PHARAOH.maternalFastaIndex": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.mat.fa.fai",
  "PHARAOH.allHifiToDiploidBam": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_winnowmap/downsample_60x/HG005.DCv1.2.60x.winnowmapv2.03.bam",
  "PHARAOH.allONTToPatBai": "/private/groups/patenlab/masri/hprc/polishing/HG005/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/UL_alignments/pat/slurm_run/HG005.pat.trio_hifiasm_0.19.5.DC_1.2.winnowmap_2.03.UL_R941_Guppy5_40x.bam.bai",
  "PHARAOH.allONTToPatBam": "/private/groups/patenlab/masri/hprc/polishing/HG005/HPRC_Y2/trio_hifiasm_0.19.5_DC_1.2/UL_alignments/pat/slurm_run/HG005.pat.trio_hifiasm_0.19.5.DC_1.2.winnowmap_2.03.UL_R941_Guppy5_40x.bam",
  "PHARAOH.sampleName": "HG005",
  "PHARAOH.paternalFastaIndex": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.pat.fa.fai"
}
```

```
mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --jobStore ./big_store --batchSystem slurm --batchLogsDir ./logs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/workflows/PHARAOH.wdl PHARAOH_wdl_inputs.json -o toil_PH_out -m toil_PH.json 2>&1 | tee toil_PH_log.txt
```

DeepPolisher model 5 winnowmap on HG005 high coverage bam file for
```
/private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_60x_PHv5/toil_PH_out
```
```
{
  "runDeepPolisher.sampleName": "HG005",
  "runDeepPolisher.ModelFilesTarGZ": "/private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-220.tar.gz",
  "runDeepPolisher.Bai": "/private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_60x_PHv5/toil_PH_out/HG005.DCv1.2.60x.winnowmapv2.03.PHARAOH.bam.bai",
  "runDeepPolisher.Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "runDeepPolisher.dockerImage": "google/deepconsensus:polisher_v0.0.7_080323",
  "runDeepPolisher.Bam": "/private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_60x_PHv5/toil_PH_out/HG005.DCv1.2.60x.winnowmapv2.03.PHARAOH.bam"
}
```
```
mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --jobStore ./big_store --batchSystem slurm --batchLogsDir ./logs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/DeepPolisher.wdl DeepPolisher_inputs.json -o DP_output -m DP_outputs.json  2>&1 | tee toil_DP_log.txt

```
DeepPolisher model 5 winnowmap on HG005 full coverage 130x bam, original before pharaoh

```
{
  "runDeepPolisher.sampleName": "HG005",
  "runDeepPolisher.ModelFilesTarGZ": "/private/groups/patenlab/mira/hprc_polishing/data/DeepPolisher_models/checkpoint-220.tar.gz",
  "runDeepPolisher.Bai": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_winnowmap/toil_winnow_full_cov_out/HG005.DCv1.2.full_coverage.winnowmapv2.03.bam.bai",
  "runDeepPolisher.Fasta": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/HG005.trio_hifiasm_0.19.5.DC_1.2_40x.dip.fa",
  "runDeepPolisher.dockerImage": "google/deepconsensus:polisher_v0.0.7_080323",
  "runDeepPolisher.Bam": "/private/groups/patenlab/mira/hprc_polishing/data/HG005_y2_polishing/alignments/HiFi_DCv1.2_winnowmap/toil_winnow_full_cov_out/HG005.DCv1.2.full_coverage.winnowmapv2.03.bam"
}
```
```
cd /private/groups/patenlab/mira/hprc_polishing/hprc_deepPolisher_wf_runs/HG005_y2_DCv1.2_130x_DPmodel5

mkdir -p logs && time SINGULARITY_CACHEDIR=`pwd`/outputs/cache/.singularity/cache MINIWDL__SINGULARITY__IMAGE_CACHE=`pwd`/outputs/cache/.cache/miniwdl toil-wdl-runner --jobStore ./big_store --batchSystem slurm --batchLogsDir ./logs /private/home/mmastora/progs/hpp_production_workflows/QC/wdl/tasks/DeepPolisher.wdl DeepPolisher_inputs.json -o DP_output -m DP_outputs.json  2>&1 | tee toil_DP_log.txt
```
