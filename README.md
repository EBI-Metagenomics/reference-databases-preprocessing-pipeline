# Reference databases preprocessing pipeline

Nextflow workflow used to preprocess various reference database files in use by MGnify analysis pipelines.

## Databases included

- [SILVA](https://www.arb-silva.de/) (SSU + LSU)
- [PR2](https://pr2-database.org/) (SSU)
- [UNITE](https://unite.ut.ee/) (ITS)
- [ITSoneDB](https://itsonedb.cloud.ba.infn.it/) (ITS)
- [Rfam](https://rfam.org/) (rRNA covariance models)
- [UniRef90](https://www.uniprot.org/help/uniref) (protein sequences)
- [Rhea](https://www.rhea-db.org/)
- [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) (taxdump)
- [KOFAM](https://www.genome.jp/tools/kofamkoala/)

## How to run

Run this command with the flags for generation of the databases you need:

```
nextflow run main.nf \
    --generate_amplicon_db \
    --generate_uniref90_db \
    --generate_kofam_db
```

`--generate_amplicon_db` enable preprocessing of ITSone, PR2, RFAM, SILVA-LSU, SILVA-SSU and UNITE databases, required for taxonomic analysis in MGnify amplicon pipeline. It generates output like this, one subdirectory for each database:

```
├── ITSone
│   └── 1.141
│       ├── ITSone.fasta
│       ├── ITSone.otu
│       └── ITSone-tax.txt
├── PR2
│   └── 5.0.0
│       ├── PR2.fasta
│       ├── PR2.otu
│       └── PR2-tax.txt
├── RFAM
│   └── 14.10
│       ├── ribo.clan_info
│       └── ribo.cm
├── SILVA-LSU
│   └── 138.1
│       ├── SILVA-LSU.fasta
│       ├── SILVA-LSU.otu
│       └── SILVA-LSU-tax.txt
├── SILVA-SSU
│   └── 138.1
│       ├── SILVA-SSU.fasta
│       ├── SILVA-SSU.otu
│       └── SILVA-SSU-tax.txt
└── UNITE
    └── 9.0
        ├── UNITE.fasta
        ├── UNITE.otu
        └── UNITE-tax.txt
```

`--generate_rhea_tax_db` enable preprocessing of UniRef90 and RHEA database to produce MGnify custom DB for RHEA reactions annotation, and also UniRef90 and NCBI taxonomy to generate CAT_pack database for taxonomic classification of contigs. It generates the following output structure:
```
├── uniref90_rhea
│   └── 2024_05
│       ├── rhea_chebi_mapping_135.tsv
│       └── uniref90_rhea_2024_05_2024-07-31.dmnd
└── uniref90_taxonomy
    └── 2024_05
        ├── db
        │   ├── catpack.dmnd
        │   ├── catpack.fastaid2LCAtaxid
        │   └── catpack.taxids_with_multiple_offspring
        └── tax
            ├── names.dmp
            └── nodes.dmp
```
`--generate_kofam_db` enable preprocessing of KOFAM Koala and generation of the database of HMM profiles for `hmmscan`. The output structure is as follows:
```
└── kofam
    └── 2025-01-22
        ├── kofam_modified.h3f
        ├── kofam_modified.h3i
        ├── kofam_modified.h3m
        └── kofam_modified.h3p
```
You will need to modify some values in the `nextflow.config` file to update a new version of the database, change output directory, etc. 
