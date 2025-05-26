# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.1] - 2025-05-26
- Re-deployment to enable labels for optional outputs
- [GRD-948](https://jira.oicr.on.ca/browse/GRD-948)

## [1.1.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - Add vidarr labels to outputs (changes to medata only).

## [1.0.4] - 2023-12-21
### Changed
- [GRD-738]
  - Remove concat fastq, align multiple fastq separately to allow add lane-level readGroup.

## [1.0.3] - 2023-12-14
### Changed
- [GRD-683] - Used new version of consensusCruncher 5.0.1 which includes bug fix of readGroup in bam file header.

## [1.0.0] - 2022-10-06
### Changed
- [GRD-440](https://jira.oicr.on.ca/browse/GRD-440) 
- Reimplementation of the consensus cruncher workflow.
- Removed some tasks from original workflow.

## [1.0.1] - 2022-10-20
### Added
- Added set default modules.
