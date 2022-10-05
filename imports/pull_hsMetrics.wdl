version 1.0

workflow hsMetrics {
input {
   Int collectHSmetrics_timeout = 5
   String collectHSmetrics_modules = "picard/2.21.2 hg19/p13"
   Int collectHSmetrics_maxRecordsInRam = 250000
   Int collectHSmetrics_coverageCap = 500
   Int collectHSmetrics_jobMemory = 18
   String collectHSmetrics_filter = "LENIENT"
   String collectHSmetrics_metricTag = "HS"
   String collectHSmetrics_refFasta = "$HG19_ROOT/hg19_random.fa"
   Int bedToBaitIntervals_timeout = 1
   String bedToBaitIntervals_modules = "picard/2.21.2 hg19/p13"
   Int bedToBaitIntervals_jobMemory = 16
   String bedToBaitIntervals_refDict = "$HG19_ROOT/hg19_random.dict"
   Int bedToTargetIntervals_timeout = 1
   String bedToTargetIntervals_modules = "picard/2.21.2 hg19/p13"
   Int bedToTargetIntervals_jobMemory = 16
   String bedToTargetIntervals_refDict = "$HG19_ROOT/hg19_random.dict"
   File    inputBam
   String  baitBed
   String  targetBed
   String outputFileNamePrefix = basename(inputBam, '.bam')
}

call bedToIntervals as bedToTargetIntervals { input: inputBed = targetBed, refDict = bedToTargetIntervals_refDict, jobMemory = bedToTargetIntervals_jobMemory, modules = bedToTargetIntervals_modules, timeout = bedToTargetIntervals_timeout }
call bedToIntervals as bedToBaitIntervals { input: inputBed = baitBed, refDict = bedToBaitIntervals_refDict, jobMemory = bedToBaitIntervals_jobMemory, modules = bedToBaitIntervals_modules, timeout = bedToBaitIntervals_timeout }

call collectHSmetrics{ input: inputBam = inputBam, baitIntervals = bedToBaitIntervals.outputIntervals, targetIntervals = bedToTargetIntervals.outputIntervals, outputPrefix = outputFileNamePrefix, refFasta = collectHSmetrics_refFasta, metricTag = collectHSmetrics_metricTag, filter = collectHSmetrics_filter, jobMemory = collectHSmetrics_jobMemory, coverageCap = collectHSmetrics_coverageCap, maxRecordsInRam = collectHSmetrics_maxRecordsInRam, modules = collectHSmetrics_modules, timeout = collectHSmetrics_timeout }

meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "HSMetrics 2.0"
  dependencies: [{
    name: "picard/2.21.2",
    url: "https://broadinstitute.github.io/picard/"
  }]
}

parameter_meta {
    collectHSmetrics_timeout: "Maximum amount of time (in hours) the task can run for."
    collectHSmetrics_modules: "Names and versions of modules needed"
    collectHSmetrics_maxRecordsInRam: "Specifies the N of records stored in RAM before spilling to disk. Increasing this number increases the amount of RAM needed."
    collectHSmetrics_coverageCap: "Parameter to set a max coverage limit for Theoretical Sensitivity calculations"
    collectHSmetrics_jobMemory: "Memory allocated to job"
    collectHSmetrics_filter: "Settings for picard filter"
    collectHSmetrics_metricTag: "Extension for metrics file"
    collectHSmetrics_refFasta: "Path to fasta reference file"
    bedToBaitIntervals_timeout: "Maximum amount of time (in hours) the task can run for."
    bedToBaitIntervals_modules: "Names and versions of modules needed"
    bedToBaitIntervals_jobMemory: "Memory allocated to job"
    bedToBaitIntervals_refDict: "Path to index of fasta reference file"
    bedToTargetIntervals_timeout: "Maximum amount of time (in hours) the task can run for."
    bedToTargetIntervals_modules: "Names and versions of modules needed"
    bedToTargetIntervals_jobMemory: "Memory allocated to job"
bedToTargetIntervals_refDict: "Path to index of fasta reference file"
}

output {
  File outputHSMetrics  = collectHSmetrics.outputHSMetrics
}

}

task bedToIntervals {
input {
   String inputBed
   String refDict = "$HG19_ROOT/hg19_random.dict"
   Int    jobMemory = 16
   String modules   = "picard/2.21.2 hg19/p13"
   Int timeout = 1
}

command <<<
 java -Xmx~{jobMemory-6}G -jar $PICARD_ROOT/picard.jar BedToIntervalList \
                              INPUT=~{inputBed} \
                              OUTPUT="~{basename(inputBed, '.bed')}.interval_list" \
                              SD="~{refDict}"
>>>

parameter_meta {
 inputBed: "Path to input bed file"
 refDict: "Path to index of fasta reference file"
 jobMemory: "Memory allocated to job"
 modules: "Names and versions of modules needed"
 timeout: "Maximum amount of time (in hours) the task can run for."
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File outputIntervals = "~{basename(inputBed, '.bed')}.interval_list"
}
}

task collectHSmetrics {
input { 
   File   inputBam
   String baitIntervals
   String targetIntervals
   String refFasta   = "$HG19_ROOT/hg19_random.fa"
   String metricTag  = "HS"
   String filter     = "LENIENT"
   String outputPrefix = "OUTPUT"
   Int   jobMemory   = 18
   Int   coverageCap = 500
   Int   maxRecordsInRam = 250000
   String modules    = "picard/2.21.2 hg19/p13"
   Int timeout = 5
}

command <<<
 java -Xmx~{jobMemory-6}G -jar $PICARD_ROOT/picard.jar CollectHsMetrics \
                              TMP_DIR=picardTmp \
                              BAIT_INTERVALS=~{baitIntervals} \
                              TARGET_INTERVALS=~{targetIntervals} \
                              R=~{refFasta} \
                              COVERAGE_CAP=~{coverageCap} \
                              MAX_RECORDS_IN_RAM=~{maxRecordsInRam} \
                              INPUT=~{inputBam} \
                              OUTPUT="~{outputPrefix}.~{metricTag}.txt" \
                              VALIDATION_STRINGENCY=~{filter}
>>>

parameter_meta {
 inputBam: "Input bam file"
 baitIntervals: "path to bed file with bait intervals"
 targetIntervals: "path to bed file with target intervals"
 refFasta: "Path to fasta reference file"
 metricTag: "Extension for metrics file"
 filter: "Settings for picard filter"
 outputPrefix: "prefix to build a name for output file"
 coverageCap: "Parameter to set a max coverage limit for Theoretical Sensitivity calculations"
 maxRecordsInRam: "Specifies the N of records stored in RAM before spilling to disk. Increasing this number increases the amount of RAM needed."
 jobMemory: "Memory allocated to job"
 modules: "Names and versions of modules needed"
 timeout: "Maximum amount of time (in hours) the task can run for."
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File outputHSMetrics = "~{outputPrefix}.~{metricTag}.txt"
}

}

