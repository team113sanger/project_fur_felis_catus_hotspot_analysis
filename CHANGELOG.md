# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
## [1.0.1] - 2025-03-31
### Changed
- Reorganised CLI to use a single entrypoint under `fur_hotspot_mutations`
  - e.g. `fur_hotspot_mutations --version`
- Previous scripts now have a CLI sub command:
  - `fur_hotspot_mutations extract_hotspot_mutations --help`
  - `fur_hotspot_mutations mpileup_variant_filter --help`
  - `fur_hotspot_mutations update_maf_variants --help`
  - `fur_hotspot_mutations check_mpileup_samples_from_maf --help`
- Small changes to make tests more portable
- Logging level can now be controlled from the CLI

### Added
- More portable public Docker image as `Dockerfile.public` with supporting CICD changes

## [1.0.0] - 2025-03-12
### Added
- Initial publication prep for Hotspot insepction tool used in Felis Catus FUR project.

## [0.1.0] - 2024-11-13
### Added
- Initial project setup using the template from [example-python-cicd](https://gitlab.internal.sanger.ac.uk/team113sanger/common/example-python-cicd).
