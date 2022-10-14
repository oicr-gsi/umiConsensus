version 1.0

import "imports/pull_hsMetrics.wdl" as hsMetrics

struct InputGroup {
  File fastqR1
  File fastqR2
}

workflow umiConsensus {
  input {
    Array[InputGroup]? inputGroups
    File? sortedBam
    File? sortedBai
    String outputFileNamePrefix
    String intervalFile
    String reference
  }

  parameter_meta {
    inputGroups: "Array of fastq files to concatenate if a top-up"
    sortedBam: "Bam file from bwamem"
    sortedBai: "Bai file from bwamem"
    outputFileNamePrefix: "Prefix to use for output file"
    intervalFile: "interval file to subset variant calls"
    reference: "the reference build of the genome"
  }
    if (reference == "hg19") {
        String hg19inputRefDict = "$HG19_ROOT/hg19_random.dict"
        String hg19inputRefFasta = "$HG19_ROOT/hg19_random.fa" 
        String hg19inputHSMetricsModules = "picard/2.21.2 hg19/p13"
        String hg19alignModules = "consensus-cruncher/5.0 data-hg19-consensus-cruncher/1.0 hg19-bwa-index/0.7.12 samtools/1.9"
        String hg19bwaref = "$HG19_BWA_INDEX_ROOT/hg19_random.fa"
        String hg19blist = "$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/IDT_duplex_sequencString ing_barcodes.list"
        String hg19consensusModules = "consensus-cruncher/5.0 data-hg19-consensusString -cruncher/1.0 hg19-bwa-index/0.7.12 samtools/1.9" 
        String hg19genome = "hg19" 
        String hg19cytoband = "$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/hg19_cytoBand.tString xt" 
    }
    if (reference == "hg38") {
        String hg38inputRefDict  = "$HG38_ROOT/hg38_random.dict"
        String hg38inputRefFasta = "$HG38_ROOT/hg38_random.fa"
        String hg38inputHSMetricsModules = "picard/2.21.2 hg38/p12"
        String hg38alignModules = "consensus-cruncher/5.0 data-hg38-consensus-cruncher/1.0 hg38-bwa-index-with-alt/0.7.12 samtools/1.9"
        String hg38bwaref = "$HG38_BWA_INDEX_WITH_ALT_ROOT/hg38_random.fa"
        String hg38blist = "$DATA_HG38_CONSENSUS_CRUNCHER_ROOT/IDT_duplex_sequencing_barcodes.list"
        String hg38consensusModules = "consensus-cruncher/5.0 data-hg38-consensus-cruncher/1.0 hg38-bwa-index-with-alt/0.7.12 samtools/1.9"
        String hg38genome = "hg38" 
        String hg38cytoband = "$DATA_HG38_CONSENSUS_CRUNCHER_ROOT/hg38_cytoBand.txt"
    }

    String inputRefDict = select_first([hg19inputRefDict, hg38inputRefDict])
    String inputRefFasta = select_first([hg19inputRefFasta, hg38inputRefFasta])
    String inputHSMetricsModules = select_first([hg19inputHSMetricsModules, hg38inputHSMetricsModules])
    String alignModules = select_first([hg19alignModules, hg38alignModules])
    String bwaref = select_first([hg19bwaref, hg38bwaref])
    String blist = select_first([hg19blist, hg38blist])
    String consensusModules = select_first([hg19consensusModules, hg38consensusModules])
    String genome = select_first([hg19genome, hg38genome])
    String cytoband = select_first([hg19cytoband, hg38cytoband])
    

if (!(defined(sortedBam)) && defined(inputGroups)) {
  Array[InputGroup] inputs = select_first([inputGroups])
  scatter (ig in inputs) {
    File read1s       = ig.fastqR1
    File read2s       = ig.fastqR2
  }
}
  
  if (!(defined(sortedBam)) && defined(inputGroups)) {
    call concat {
      input:
        read1s = select_first([read1s]),
        read2s = select_first([read2s]),
        outputFileNamePrefix = outputFileNamePrefix
    }
  }

  if (!(defined(sortedBam)) && defined(read1s) && defined(read2s)) {
    call align {
      input:
        fastqR1 = select_first([concat.fastqR1]),
        fastqR2 = select_first([concat.fastqR2]),
        outputFileNamePrefix = outputFileNamePrefix,
        modules = alignModules,
        bwaref = bwaref,
        blist = blist
    }
  }

  call consensus {
    input:
      inputBam = select_first([sortedBam, align.sortedBam]),
      inputBai = select_first([sortedBai, align.sortedBai]),
      basePrefix = outputFileNamePrefix,
      modules = consensusModules,
      genome = genome,
      cytoband = cytoband
  }


  call hsMetrics.hsMetrics as hsMetricsRunDCSSC {
    input: 
      inputBam = consensus.dcsScBam,
      outputFileNamePrefix = "~{outputFileNamePrefix}.dcsSc-hsMetrics",
      baitBed = intervalFile, 
      targetBed = intervalFile,
      collectHSmetrics_modules = inputHSMetricsModules,
      collectHSmetrics_refFasta = inputRefFasta,
      bedToBaitIntervals_refDict = inputRefDict,
      bedToBaitIntervals_modules = inputHSMetricsModules,
      bedToTargetIntervals_refDict = inputRefDict,
      bedToTargetIntervals_modules = inputHSMetricsModules
  }

  call hsMetrics.hsMetrics as hsMetricsRunSSCSSC {
    input: 
      inputBam = consensus.sscsScBam,
      outputFileNamePrefix = "~{outputFileNamePrefix}.sscsSc-hsMetrics",
      baitBed = intervalFile, 
      targetBed = intervalFile,
      collectHSmetrics_modules = inputHSMetricsModules,
      collectHSmetrics_refFasta = inputRefFasta,
      bedToBaitIntervals_refDict = inputRefDict,
      bedToBaitIntervals_modules = inputHSMetricsModules,
      bedToTargetIntervals_refDict = inputRefDict,
      bedToTargetIntervals_modules = inputHSMetricsModules
  }

  call hsMetrics.hsMetrics as hsMetricsRunAllUnique {
    input: 
      inputBam = consensus.allUniqueBam,
      outputFileNamePrefix = "~{outputFileNamePrefix}.allUnique-hsMetrics",
      baitBed = intervalFile, 
      targetBed = intervalFile,
      collectHSmetrics_modules = inputHSMetricsModules,
      collectHSmetrics_refFasta = inputRefFasta,
      bedToBaitIntervals_refDict = inputRefDict,
      bedToBaitIntervals_modules = inputHSMetricsModules,
      bedToTargetIntervals_refDict = inputRefDict,
      bedToTargetIntervals_modules = inputHSMetricsModules
  }

  meta {
    author: "Alexander Fortuna, Rishi Shah and Gavin Peng"
    email: "alexander.fortuna@oicr.on.ca, rshah@oicr.on.ca and gpeng@oicr.on.ca"
    description: "Workflow to run extract UMIs from fastq and generate consensus Bams as well as run it thru mutect2 task and combinevariants task"
    dependencies: [
     {
      name: "hg19-bwa-index/0.7.12",
      url: "http://bio-bwa.sourceforge.net/"
     },
     {
      name: "samtools/1.9",
      url: "http://www.htslib.org/"
     },
     {
      name: "python/3.6",
      url: "https://www.python.org/downloads/"
     },
     {
      name: "picard/2.21.2",
      url: "https://broadinstitute.github.io/picard/"
     },
     {
      name: "rstats/3.6",
      url: "https://www.r-project.org/"
     },
     {
      name: "consensuscruncer-5.0",
      url: "https://github.com/pughlab/ConsensusCruncher"
     }
    ]
    output_meta: {
      rawBam: "aligned bam file",
      rawBamIndex: "aligned bam index",
      dcsScBam: "DCS generated from SSCS + SC",
      dcsScBamIndex: "Index for DCS SC Bam",
      allUniqueBam: "DCS (from SSCS + SC) + SSCS_SC_Singletons + remaining singletons",
      allUniqueBamIndex: "Index for All Unique Bam",
      sscsScBam: "SSCS combined with corrected singletons (from both rescue strategies)",
      sscsScBamIndex: "Index for SSCS SC Bam",
      ccFolder: "output folder containing files not needed for downstream analysis; info on family size, QC metrics",
      outputCCStats: " Consensus sequence formation metrics",
      outputCCReadFamilies: "Family size and frequency from consensusCruncher",
      dcsScHsMetrics: "Hs Metrics for duplex consensus sequences (DCS)",
      sscsScHsMetrics: "HS Metrics for single-strand consensus sequences (SSCS)",
      allUniqueHsMetrics: "HS Metrics for AllUnique"
    }

  }
  
  output {
    File? rawBam = align.sortedBam
    File? rawBamIndex = align.sortedBai
    File dcsScBam = consensus.dcsScBam
    File dcsScBamIndex = consensus.dcsScBamIndex
    File allUniqueBam = consensus.allUniqueBam
    File allUniqueBamIndex = consensus.allUniqueBamIndex
    File sscsScBam = consensus.sscsScBam
    File sscsScBamIndex = consensus.sscsScBamIndex
    File outputCCStats = consensus.statsCCFile
    File outputCCReadFamilies = consensus.readFamiliesCCFile
    File ccFolder = consensus.ccFolder
    File dcsScHsMetrics = hsMetricsRunDCSSC.outputHSMetrics
    File sscsScHsMetrics = hsMetricsRunSSCSSC.outputHSMetrics
    File allUniqueHsMetrics = hsMetricsRunAllUnique.outputHSMetrics
  }
}

task concat {
  input {
    Array[File]+ read1s
    Array[File]+ read2s
    String outputFileNamePrefix
    Int threads = 4
    Int jobMemory = 16
    Int timeout = 72
    String modules = "tabix/0.2.6"
  }

  parameter_meta {
    read1s: "array of read1s"
    read2s: "array of read2s"
    outputFileNamePrefix: "File name prefix"
    threads: "Number of threads to request"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    modules: "Required environment modules"
  }

  command <<<
    set -euo pipefail

    zcat ~{sep=" " read1s} | gzip > ~{outputFileNamePrefix}_R1_001.fastq.gz

    zcat ~{sep=" " read2s} | gzip > ~{outputFileNamePrefix}_R2_001.fastq.gz

  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
    modules: "~{modules}"

  }

  output {
    File fastqR1 = "~{outputFileNamePrefix}_R1_001.fastq.gz"
    File fastqR2 = "~{outputFileNamePrefix}_R2_001.fastq.gz"
  }
}

task align {
  input {
    File fastqR1
    File fastqR2
    String outputFileNamePrefix
    String modules 
    String consensusCruncherPy = "$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"
    String bwa = "$BWA_ROOT/bin/bwa"
    String bwaref 
    String samtools = "$SAMTOOLS_ROOT/bin/samtools"
    String blist 
    Int threads = 4
    Int jobMemory = 16
    Int timeout = 72
  }

  parameter_meta {
    fastqR1: "Path to left fastq file"
    fastqR2: "Path to right fastq file"
    outputFileNamePrefix: "File name prefix"
    consensusCruncherPy: "Path to consensusCruncher binary"
    modules: "Names and versions of modules to load"
    bwa: "Path to bwa binary"
    bwaref: "Path to bwa index"
    samtools: "Path to samtools binary"
    blist: "Path to blacklist for barcodes"
    threads: "Number of threads to request"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail

    ~{consensusCruncherPy} fastq2bam \
         --fastq1 ~{fastqR1} \
         --fastq2 ~{fastqR2}\
         --output . \
         --bwa ~{bwa} \
         --ref ~{bwaref} \
         --samtools ~{samtools} \
         --skipcheck \
         --blist ~{blist}

    # Necessary for if bam files to be named according to merged library name
    # Additionally if ".sorted" isn't omitted here, file names from align include ".sorted" twice
    mv bamfiles/*.bam bamfiles/"~{outputFileNamePrefix}.bam"
    mv bamfiles/*.bai bamfiles/"~{outputFileNamePrefix}.bam.bai"
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    cpu:     "~{threads}"
    timeout: "~{timeout}"

  }

  output {
    File? sortedBam = "bamfiles/~{outputFileNamePrefix}.bam"
    File? sortedBai = "bamfiles/~{outputFileNamePrefix}.bam.bai"
  }
}

task consensus {
  input {
    File? inputBam
    File? inputBai
    String consensusCruncherPy = "$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"
    String basePrefix
    String samtools = "$SAMTOOLS_ROOT/bin/samtools"
    String cytoband
    String genome  
    String ccDir = basePrefix + ".consensuscruncher"
    Float cutoff  = 0.7
    Int threads = 8
    Int jobMemory = 32
    Int timeout = 72
    String modules 
  }

  parameter_meta {
    inputBam: "Bam file either from bwamem or ConsensusCruncher align."
    inputBai: "Bai file either from bwamem or ConsensusCruncher align."
    consensusCruncherPy: "Path to consensusCruncher binary"
    modules: "Names and versions of modules to load"
    samtools: "Path to samtools binary"
    threads: "Number of threads to request"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    genome: "Which genome version to use"
    cytoband: "Path to cytoband for genome"
    cutoff: "Cutoff to use to call a consenus of reads"
    ccDir: "Placeholder"
  }


  command <<<
  set -euo pipefail

   ~{consensusCruncherPy} consensus \
         --input ~{inputBam} \
         --output . \
         --samtools ~{samtools} \
         --cutoff ~{cutoff} \
         --genome ~{genome} \
         --bedfile ~{cytoband} \
         --bdelim '|'

   tar cf - ~{basePrefix} | gzip --no-name > ~{ccDir}.tar.gz
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File dcsScBam = "~{basePrefix}/dcs_sc/~{basePrefix}.dcs.sc.sorted.bam"
    File dcsScBamIndex = "~{basePrefix}/dcs_sc/~{basePrefix}.dcs.sc.sorted.bam.bai"
    File allUniqueBam = "~{basePrefix}/dcs_sc/~{basePrefix}.all.unique.dcs.sorted.bam"
    File allUniqueBamIndex = "~{basePrefix}/dcs_sc/~{basePrefix}.all.unique.dcs.sorted.bam.bai"
    File sscsScBam = "~{basePrefix}/sscs_sc/~{basePrefix}.sscs.sc.sorted.bam"
    File sscsScBamIndex = "~{basePrefix}/sscs_sc/~{basePrefix}.sscs.sc.sorted.bam.bai"
    File statsCCFile = "~{basePrefix}/~{basePrefix}.stats.txt"
    File readFamiliesCCFile = "~{basePrefix}/~{basePrefix}.read_families.txt"
    File ccFolder = "~{ccDir}.tar.gz"
    
  }
}
