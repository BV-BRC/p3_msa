
# Application specification: MSA

This is the application specification for service with identifier MSA.

The backend script implementing the application is [App-MSA.pl](../service-scripts/App-MSA.pl).

The raw JSON file for this specification is [MSA.json](MSA.json).

This service performs the following task:   Compute the multiple sequence alignment and analyze SNP/variance.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| input_status | The input status. | enum  |  |  |
| input_type | The type of input selected. | enum  |  |  |
| fasta_files |  | group  |  |  |
| select_genomegroup | Genome groups | WS: genome_group  |  |  |
| feature_groups | Feature groups | WS: feature_group  |  |  |
| feature_list | feature list | string  |  |  |
| genome_list | genome list | string  |  |  |
| aligner | multiple sequence aligner | enum  |  | Muscle |
| alphabet | sequence alphabet | enum  | :heavy_check_mark: | dna |
| fasta_keyboard_input | fasta keyboard input | string  |  |  |
| ref_type | Reference sequence type | enum  |  | none |
| strategy |  |   |  | auto |
| ref_string | reference sequence identity | string  |  |  |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |

