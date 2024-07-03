# MSA and SNP/Variation Analysis

## Overview

The Multiple Sequence Alignment (MSA) and Single Nucleotide Polymorphism (SNP)/Variation Analysis Service allows users to choose an alignment algorithm to align sequences selected from: a search result, a FASTA file saved to the workspace, or through simply cutting and pasting. The service can also be used for variation and SNP analysis with feature groups, FASTA files, aligned FASTA files, and user input FASTA records. If a single alignment file is given, then only the variation analysis is run. If multiple inputs are given, the program concatenates all sequence records and aligns them. If a mixture of protein and nucleotides are given, then nucleotides are converted to proteins.

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

There is one application service specifications defined here:

1. [MSA](app_specs/MSA.md): MSA service Performs a multiple sequence alignment with entropy and variation analysis.

The code in this module provides the BV-BRC application service wrapper scripts for the MSA service as well
as some backend utilities:

| Script name | Purpose |
| ----------- | ------- |
| [App-MSA.pl](service-scripts/App-MSA.pl) | App script for the [MSA and SNP/Variation Analysis Service](https://www.bv-brc.org/docs/quick_references/services/msa_snp_variation_service.html) |

## See also

* [MSA and SNP/Variation Analysis Service](https://www.bv-brc.org/app/MSA)
* [Quick Reference](https://www.bv-brc.org/docs/quick_references/services/msa_snp_variation_service.html)
* [MSA and SNP/Variation Analysis Service Tutorial](https://www.bv-brc.org/docs/tutorial/msa_snp_variation/msa_snp_variation.html)
* [MSA and SNP/Variation Analysis Service Instructional Video](https://www.youtube.com/watch?v=ea6GboAZPQs&feature=youtu.be)


## References

Katoh K., Misawa K., Kuma K., Miyata T. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res. 2002;30:3059–3066.

Katoh K, Rozewicki J, Yamada KD. MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization. Brief Bioinform. 2019 Jul 19;20(4):1160-1166. doi: 10.1093/bib/bbx108. PMID: 28968734; PMCID: PMC6781576.

Edgar, Robert C. (2004), MUSCLE: multiple sequence alignment with high accuracy and high throughput, Nucleic Acids Research 32(5), 1792-97.

Edgar, Robert C (2004), MUSCLE: a multiple sequence alignment method with reduced time and space complexity. BMC Bioinformatics, 5(1):113.

Darling AE, Mau B, Perna NT (2010) progressiveMauve: Multiple Genome Alignment with Gene Gain, Loss and Rearrangement. PLoS ONE 5(6): e11147. https://doi.org/10.1371/journal.pone.0011147
