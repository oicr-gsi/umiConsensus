# umiConsensus

Workflow to run extract UMIs from fastq and generate consensus Bams as well as run it thru mutect2 task and combinevariants task

## Overview

## Dependencies

* [hg19-bwa-index 0.7.12](http://bio-bwa.sourceforge.net/)
* [samtools 1.9](http://www.htslib.org/)
* [python 3.6](https://www.python.org/downloads/)
* [picard 2.21.2](https://broadinstitute.github.io/picard/)
* [rstats 3.6](https://www.r-project.org/)
* [consensuscruncer-5.0.1](https://github.com/oicr-gsi/ConsensusCruncher)


## Usage

### Cromwell
```
java -jar cromwell.jar run umiConsensus.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`outputFileNamePrefix`|String|Prefix to use for output file
`intervalFile`|String|interval file to subset variant calls
`reference`|String|the reference build of the genome


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fastqGroups`|Array[fastqGroup]?|None|Array of fastq files to concatenate if a top-up
`sortedBam`|File?|None|Bam file from bwamem
`sortedBai`|File?|None|Bai file from bwamem


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`align.consensusCruncherPy`|String|"$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"|Path to consensusCruncher binary
`align.bwa`|String|"$BWA_ROOT/bin/bwa"|Path to bwa binary
`align.samtools`|String|"$SAMTOOLS_ROOT/bin/samtools"|Path to samtools binary
`align.threads`|Int|4|Number of threads to request
`align.jobMemory`|Int|16|Memory allocated for this job
`align.timeout`|Int|72|Hours before task timeout
`mergeBams.additionalParams`|String?|None|Additional parameters to pass to GATK MergeSamFiles.
`mergeBams.jobMemory`|Int|48|Memory allocated to job (in GB).
`mergeBams.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`mergeBams.cores`|Int|1|The number of cores to allocate to the job.
`mergeBams.timeout`|Int|8|Maximum amount of time (in hours) the task can run for.
`mergeBams.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`consensus.consensusCruncherPy`|String|"$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"|Path to consensusCruncher binary
`consensus.samtools`|String|"$SAMTOOLS_ROOT/bin/samtools"|Path to samtools binary
`consensus.ccDir`|String|basePrefix + ".consensuscruncher"|Placeholder
`consensus.cutoff`|Float|0.7|Cutoff to use to call a consenus of reads
`consensus.threads`|Int|8|Number of threads to request
`consensus.jobMemory`|Int|32|Memory allocated for this job
`consensus.timeout`|Int|72|Hours before task timeout
`hsMetricsRunDCSSC.collectHSmetrics_timeout`|Int|5|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunDCSSC.collectHSmetrics_maxRecordsInRam`|Int|250000|Specifies the N of records stored in RAM before spilling to disk. Increasing this number increases the amount of RAM needed.
`hsMetricsRunDCSSC.collectHSmetrics_coverageCap`|Int|500|Parameter to set a max coverage limit for Theoretical Sensitivity calculations
`hsMetricsRunDCSSC.collectHSmetrics_jobMemory`|Int|18|Memory allocated to job
`hsMetricsRunDCSSC.collectHSmetrics_filter`|String|"LENIENT"|Settings for picard filter
`hsMetricsRunDCSSC.collectHSmetrics_metricTag`|String|"HS"|Extension for metrics file
`hsMetricsRunDCSSC.bedToBaitIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunDCSSC.bedToBaitIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunDCSSC.bedToTargetIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunDCSSC.bedToTargetIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunSSCSSC.collectHSmetrics_timeout`|Int|5|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunSSCSSC.collectHSmetrics_maxRecordsInRam`|Int|250000|Specifies the N of records stored in RAM before spilling to disk. Increasing this number increases the amount of RAM needed.
`hsMetricsRunSSCSSC.collectHSmetrics_coverageCap`|Int|500|Parameter to set a max coverage limit for Theoretical Sensitivity calculations
`hsMetricsRunSSCSSC.collectHSmetrics_jobMemory`|Int|18|Memory allocated to job
`hsMetricsRunSSCSSC.collectHSmetrics_filter`|String|"LENIENT"|Settings for picard filter
`hsMetricsRunSSCSSC.collectHSmetrics_metricTag`|String|"HS"|Extension for metrics file
`hsMetricsRunSSCSSC.bedToBaitIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunSSCSSC.bedToBaitIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunSSCSSC.bedToTargetIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunSSCSSC.bedToTargetIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunAllUnique.collectHSmetrics_timeout`|Int|5|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunAllUnique.collectHSmetrics_maxRecordsInRam`|Int|250000|Specifies the N of records stored in RAM before spilling to disk. Increasing this number increases the amount of RAM needed.
`hsMetricsRunAllUnique.collectHSmetrics_coverageCap`|Int|500|Parameter to set a max coverage limit for Theoretical Sensitivity calculations
`hsMetricsRunAllUnique.collectHSmetrics_jobMemory`|Int|18|Memory allocated to job
`hsMetricsRunAllUnique.collectHSmetrics_filter`|String|"LENIENT"|Settings for picard filter
`hsMetricsRunAllUnique.collectHSmetrics_metricTag`|String|"HS"|Extension for metrics file
`hsMetricsRunAllUnique.bedToBaitIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunAllUnique.bedToBaitIntervals_jobMemory`|Int|16|Memory allocated to job
`hsMetricsRunAllUnique.bedToTargetIntervals_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`hsMetricsRunAllUnique.bedToTargetIntervals_jobMemory`|Int|16|Memory allocated to job


### Outputs

Output | Type | Description | Labels
---|---|---|---
`rawBam`|File?|aligned bam file|vidarr_label: rawBam
`rawBamIndex`|File?|aligned bam index|vidarr_label: rawBamIndex
`dcsScBam`|File|DCS generated from SSCS + SC|vidarr_label: dcsScBam
`dcsScBamIndex`|File|Index for DCS SC Bam|vidarr_label: dcsScBamIndex
`allUniqueBam`|File|DCS (from SSCS + SC) + SSCS_SC_Singletons + remaining singletons|vidarr_label: allUniqueBam
`allUniqueBamIndex`|File|Index for All Unique Bam|vidarr_label: allUniqueBamIndex
`sscsScBam`|File|SSCS combined with corrected singletons (from both rescue strategies)|vidarr_label: sscsScBam
`sscsScBamIndex`|File|Index for SSCS SC Bam|vidarr_label: sscsScBamIndex
`outputCCStats`|File| Consensus sequence formation metrics|vidarr_label: outputCCStats
`outputCCReadFamilies`|File|Family size and frequency from consensusCruncher|vidarr_label: outputCCReadFamilies
`ccFolder`|File|output folder containing files not needed for downstream analysis; info on family size, QC metrics|vidarr_label: ccFolder
`dcsScHsMetrics`|File|Hs Metrics for duplex consensus sequences (DCS)|vidarr_label: dcsScHsMetrics
`sscsScHsMetrics`|File|HS Metrics for single-strand consensus sequences (SSCS)|vidarr_label: sscsScHsMetrics
`allUniqueHsMetrics`|File|HS Metrics for AllUnique|vidarr_label: allUniqueHsMetrics


## Commands
 This section lists command(s) run by umiConsensus workflow
 
 * Running umiConsensus
 
 === Description here ===.
 Commands for running concat.
 
 ```
     set -euo pipefail
 
     zcat ~{sep=" " read1s} | gzip > ~{outputFileNamePrefix}_R1_001.fastq.gz
 
     zcat ~{sep=" " read2s} | gzip > ~{outputFileNamePrefix}_R2_001.fastq.gz
 
 ```
 Commands for running align
 ```
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
  ```
 Commands for running consensus:
 ```
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
   ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
