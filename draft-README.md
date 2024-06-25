# MSA and SNP/Variation Analysis Service

## Overview

The Multiple Sequence Alignment (MSA) and Single Nucleotide Polymorphism (SNP)/Variation Analysis Service allows users to choose an alignment algorithm to align sequences selected from: a search result, a FASTA file saved to the workspace, or through simply cutting and pasting. The service can also be used for variation and SNP analysis with feature groups, FASTA files, aligned FASTA files, and user input FASTA records. If a single alignment file is given, then only the variation analysis is run. If multiple inputs are given, the program concatenates all sequence records and aligns them. If a mixture of protein and nucleotides are given, then nucleotides are converted to proteins. 



## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

This module provides the following application specfication(s):
* [MSA](app_specs/MSA.md)


## See also

* [MSA and SNP/Variation Analysis Service](https://www.bv-brc.org/docs/https://bv-brc.org/app/MSA.html)
* [MSA and SNP/Variation Analysis Service Tutorial](https://www.bv-brc.org/docs//tutorial/msa_snp_variation/msa_snp_variation.html)
* [MSA Viewer](https://www.bv-brc.org/docs//quick_references/other/msa_viewer.html)



## References

1.	Katoh K., Misawa K., Kuma K., Miyata T. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res. 2002;30:3059–3066. 
2.	Katoh K, Rozewicki J, Yamada KD. MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization. Brief Bioinform. 2019 Jul 19;20(4):1160-1166. doi: 10.1093/bib/bbx108. PMID: 28968734; PMCID: PMC6781576.
3.	Edgar, Robert C. (2004), MUSCLE: multiple sequence alignment with high accuracy and high throughput, Nucleic Acids Research 32(5), 1792-97. 
4.	Edgar, Robert C (2004), MUSCLE: a multiple sequence alignment method with reduced time and space complexity. BMC Bioinformatics, 5(1):113.
5.	Darling AE, Mau B, Perna NT (2010) progressiveMauve: Multiple Genome Alignment with Gene Gain, Loss and Rearrangement. PLoS ONE 5(6): e11147. https://doi.org/10.1371/journal.pone.0011147
