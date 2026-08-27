# Changelog

## 1.0.0 (2026-08-27)


### Features

* add coverage and metagene coverage rules to QC workflow ([a00539c](https://github.com/jsimonas/parts-seq/commit/a00539cd829a335d2159860ed9c7a8660472f9aa))
* enable mirtrace feature in configuration and update related rules ([cf35901](https://github.com/jsimonas/parts-seq/commit/cf3590150a8e401d0ee50dd5d6b4ef7b06e76e6f))


### Bug Fixes

* add mirna_gtf parameter to aggregate_mirtop_counts rule for improved functionality ([b271982](https://github.com/jsimonas/parts-seq/commit/b2719822273ece1e8d34ceb083882e891a43643e))
* add spacing for readability and maintain consistency in Snakefile and rules ([63e36e4](https://github.com/jsimonas/parts-seq/commit/63e36e42de620ce45c6e7a47c2238a56b29fe714))
* adjust mirtop species parameter handling to improve command execution ([2f8d04c](https://github.com/jsimonas/parts-seq/commit/2f8d04ce7e7cac526b8b1ed7694aae6663b0b17a))
* correct indentation in mirtrace output expansion for clarity ([8174589](https://github.com/jsimonas/parts-seq/commit/8174589996eec2cb0b06c766c1e665e787501531))
* correct min_len parameter assignment in trim_reads rule ([44d2354](https://github.com/jsimonas/parts-seq/commit/44d2354b615a684c78a536395a63d80f68dbf654))
* correct mirtrace output expansion in all rule ([e00260f](https://github.com/jsimonas/parts-seq/commit/e00260f360035777aa36c7f7c2fb70ae8f6edf36))
* ensure mirtrace enabled parameter is correctly evaluated as boolean ([7f98c41](https://github.com/jsimonas/parts-seq/commit/7f98c41b65d4937b65881bce6e170e59b05fff26))
* improve mirtrace output expansion formatting for consistency ([a1265df](https://github.com/jsimonas/parts-seq/commit/a1265dfb674fe1f34bd9bbbd3a6da1a9a9237973))
* refactor starsolo_to_multiqc script to streamline file handling and remove unused glob import ([e344fb4](https://github.com/jsimonas/parts-seq/commit/e344fb4bbb04e755bc06dfad23356314643753c8))
* remove redundant mirtrace input expansion in multiqc rule ([70f1cd8](https://github.com/jsimonas/parts-seq/commit/70f1cd8fc2e21de02d6c40b8bde929fe57587a01))
* replace config star features with PRIMARY_STAR_FEATURE for consistency in mirtop rules ([f3c2ea0](https://github.com/jsimonas/parts-seq/commit/f3c2ea06c0f1cf19de4e8b9300557653a4b60a6b))
* streamline file path handling in QC rules for consistency ([fe6386a](https://github.com/jsimonas/parts-seq/commit/fe6386a4cfd3054434b773d8e1298c51d4b10028))
* streamline mirtop_counts_per_barcode rule by removing unnecessary log directory creation ([ee8701b](https://github.com/jsimonas/parts-seq/commit/ee8701b9faae16b9a80bfc6603d81fc3bd72383b))
* update format_starsolo rule to use PRIMARY_STAR_FEATURE for consistency ([27a1f6e](https://github.com/jsimonas/parts-seq/commit/27a1f6eeae9d6ade93e4abbeb6fbeb6b01e8ce45))
* update metagene coverage section in multiqc config and adjust CSV header handling in metagene coverage script ([45c18ef](https://github.com/jsimonas/parts-seq/commit/45c18ef42d506c54afcc923722c9a4eaf2a281b0))
* update mirtop species parameter handling to use sps_flag ([7fbf530](https://github.com/jsimonas/parts-seq/commit/7fbf5307d6b380d7c43eb3e47ce1b0d832b20695))
* update sps_flag parameter handling for mirtop rule to improve robustness ([1cb24cc](https://github.com/jsimonas/parts-seq/commit/1cb24ccf1d98718a43080886adf0c1ee839251dc))
* update stage_fastqs docstring and restore merge_fastq rule ([6dd67da](https://github.com/jsimonas/parts-seq/commit/6dd67da9ba15ab47d4dbc003f231c5ddc2e6f24f))
* update starsolo rule to enhance multi-mapper handling with EM option ([06752f5](https://github.com/jsimonas/parts-seq/commit/06752f54e7473ed2edc8325b608e8d2c3b00b808))
* update starsolo_to_multiqc script to use feature parameter for file pattern matching ([985d596](https://github.com/jsimonas/parts-seq/commit/985d59638668c3eaab584c186e21ca2e7e3e053a))
* update write_mqc_linegraph to output YAML format and adjust header structure ([8409b1c](https://github.com/jsimonas/parts-seq/commit/8409b1c0968252626b38c9aaf6a040bd437ade2d))
