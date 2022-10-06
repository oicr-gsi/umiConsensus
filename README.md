# umiConsensus

Workflow to run extract UMIs from fastq and generate consensus Bams as well as run it thru mutect2 task and combinevariants task

## Overview

## Dependencies

* [hg19-bwa-index 0.7.12](http://bio-bwa.sourceforge.net/)
* [samtools 1.9](http://www.htslib.org/)
* [python 3.6](https://www.python.org/downloads/)
* [picard 2.21.2](https://broadinstitute.github.io/picard/)
* [rstats 3.6](https://www.r-project.org/)
* [consensuscruncer-5.0](https://github.com/pughlab/ConsensusCruncher)


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
`inputRefDict`|String|reference dictionary
`inputRefFasta`|String|reference fasta file
`inputHSMetricsModules`|String|module for HSmetrics


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`inputGroups`|Array[InputGroup]?|None|Array of fastq files to concatenate if a top-up
`sortedBam`|File?|None|Bam file from bwamem
`sortedBai`|File?|None|Bai file from bwamem


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`concat.threads`|Int|4|Number of threads to request
`concat.jobMemory`|Int|16|Memory allocated for this job
`concat.timeout`|Int|72|Hours before task timeout
`concat.modules`|String|"tabix/0.2.6"|Required environment modules
`align.modules`|String|"consensus-cruncher/5.0 data-hg19-consensus-cruncher/1.0 hg19-bwa-index/0.7.12 samtools/1.9"|Names and versions of modules to load
`align.consensusCruncherPy`|String|"$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"|Path to consensusCruncher binary
`align.bwa`|String|"$BWA_ROOT/bin/bwa"|Path to bwa binary
`align.bwaref`|String|"$HG19_BWA_INDEX_ROOT/hg19_random.fa"|Path to bwa index
`align.samtools`|String|"$SAMTOOLS_ROOT/bin/samtools"|Path to samtools binary
`align.blist`|String|"$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/IDT_duplex_sequencing_barcodes.list"|Path to blacklist for barcodes
`align.threads`|Int|4|Number of threads to request
`align.jobMemory`|Int|16|Memory allocated for this job
`align.timeout`|Int|72|Hours before task timeout
`consensus.consensusCruncherPy`|String|"$CONSENSUS_CRUNCHER_ROOT/bin/ConsensusCruncher.py"|Path to consensusCruncher binary
`consensus.samtools`|String|"$SAMTOOLS_ROOT/bin/samtools"|Path to samtools binary
`consensus.cytoband`|String|"$DATA_HG19_CONSENSUS_CRUNCHER_ROOT/hg19_cytoBand.txt"|Path to cytoband for genome
`consensus.genome`|String|"hg19"|Which genome version to use
`consensus.ccDir`|String|basePrefix + ".consensuscruncher"|Placeholder
`consensus.cutoff`|Float|0.7|Cutoff to use to call a consenus of reads
`consensus.threads`|Int|8|Number of threads to request
`consensus.jobMemory`|Int|32|Memory allocated for this job
`consensus.timeout`|Int|72|Hours before task timeout
`consensus.modules`|String|"consensus-cruncher/5.0 data-hg19-consensus-cruncher/1.0 hg19-bwa-index/0.7.12 samtools/1.9"|Names and versions of modules to load
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

Output | Type | Description
---|---|---
`rawBam`|File?|aligned bam file
`rawBamIndex`|File?|aligned bam index
`dcsScBam`|File|DCS generated from SSCS + SC
`dcsScBamIndex`|File|Index for DCS SC Bam
`allUniqueBam`|File|DCS (from SSCS + SC) + SSCS_SC_Singletons + remaining singletons
`allUniqueBamIndex`|File|Index for All Unique Bam
`sscsScBam`|File|SSCS combined with corrected singletons (from both rescue strategies)
`sscsScBamIndex`|File|Index for SSCS SC Bam
`outputCCStats`|File| Consensus sequence formation metrics
`outputCCReadFamilies`|File|Family size and frequency from consensusCruncher
`ccFolder`|File|output folder containing files not needed for downstream analysis; info on family size, QC metrics
`dcsScHsMetrics`|File|Hs Metrics for duplex consensus sequences (DCS)
`sscsScHsMetrics`|File|HS Metrics for single-strand consensus sequences (SSCS)
`allUniqueHsMetrics`|File|HS Metrics for AllUnique

## Commands
This section lists command(s) run by umiConsensus workflow

* Running umiConsensus

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
