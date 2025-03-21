include { UNIREF90_RHEA_FILTER      } from '../../modules/local/uniref90_rhea_filter/main.nf'
include { UNIREF90_TAXONOMY_FILTER } from '../../modules/local/uniref90_taxonomy_filter/main.nf'
include { REFORMAT_RHEA_CHEBI       } from '../../modules/local/reformat_rhea_chebi/main.nf'
include { EXTRACT_TAXDUMP           } from '../../modules/local/extract_taxdump/main.nf'
include { DIAMOND_MAKEDB            } from '../../modules/nf-core/diamond/makedb/main.nf'
include { CATPACK_PREPARE           } from '../../modules/nf-core/catpack/prepare/main.nf'

workflow RHEA_AND_TAXONOMY_GENERATION {
    take:
        uniref90_fasta
        uniprot_rhea_mapping
        txt_rhea_chebi_mapping
        ncbi_taxdump

    main:
        UNIREF90_RHEA_FILTER(uniref90_fasta, uniprot_rhea_mapping)

        UNIREF90_RHEA_FILTER.out.filtered_proteins
            .map { filepath -> 
                [[id: "uniref90_rhea_${params.uniref90_version}_${params.uniprotKB_access_date}"], filepath]
            }
            .set { uniref90_with_rhea_ch }
        
        DIAMOND_MAKEDB(uniref90_with_rhea_ch, [], [], [])

        Channel.fromPath(uniref90_fasta)
            .splitFasta(by: params.uniref90_batch_size, file: 'uniref90')
            .set { uniref90_batches_ch }

        UNIREF90_TAXONOMY_FILTER(uniref90_batches_ch)

        UNIREF90_TAXONOMY_FILTER.out.filtered_proteins
            .collectFile(name: 'uniref90_filtered.fasta')
            .map { filepath ->
                [[id: "${params.uniref90_version}"], filepath]
            }
            .set { uniref90_filtered_ch }

        UNIREF90_TAXONOMY_FILTER.out.protid2taxid
            .collectFile(name: 'uniref90_filtered.protid2taxid')
            .set { uniref90_filtered_mapping_ch }

        EXTRACT_TAXDUMP(ncbi_taxdump)

        CATPACK_PREPARE(
            uniref90_filtered_ch,
            EXTRACT_TAXDUMP.out.tax_names,
            EXTRACT_TAXDUMP.out.tax_nodes,
            uniref90_filtered_mapping_ch
        )

        REFORMAT_RHEA_CHEBI(txt_rhea_chebi_mapping)

    emit:
        rhea_db            = DIAMOND_MAKEDB.out.db
        cat_diamond_db     = CATPACK_PREPARE.out.db
        cat_taxonomy_db    = CATPACK_PREPARE.out.taxonomy
        rhea_chebi_mapping = REFORMAT_RHEA_CHEBI.out.tsv_rhea_chebi_mapping
}
