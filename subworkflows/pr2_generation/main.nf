
include { UNCOMPRESS_FILE as UNCOMPRESS_FASTA} from '../../modules/local/uncompress_file/main.nf'
include { UNCOMPRESS_FILE as UNCOMPRESS_TAX} from '../../modules/local/uncompress_file/main.nf'
include { PR2_PROCESS_TAX } from '../../modules/local/pr2_process_tax/main.nf'
include { MAKE_OTU_FILE } from '../../modules/local/make_otu_file/main.nf'
include { CLEAN_FASTA } from '../../modules/local/clean_fasta/main.nf'
include { GENERATE_MSCLUSTER } from '../../modules/local/generate_mscluster/main.nf'

workflow PR2_GENERATION {

    main:

        dummy_fasta = file(params.dummy_fasta)
        pr2_fasta = file(params.pr2_download_fasta)
        pr2_tax = file(params.pr2_download_tax)

        UNCOMPRESS_FASTA(
            pr2_fasta,
            "PR2.fasta"
        )

        UNCOMPRESS_TAX(
            pr2_tax,
            "PR2.tax"
        )

        PR2_PROCESS_TAX(
            UNCOMPRESS_TAX.out.uncmp_file,
            params.pr2_tax_header,
            params.pr2_version,
            params.pr2_label
        )

        MAKE_OTU_FILE(
            PR2_PROCESS_TAX.out.tax,
            params.empty_file,
            params.pr2_version,
            params.pr2_label
        )

        CLEAN_FASTA(
            UNCOMPRESS_FASTA.out.uncmp_file,
            PR2_PROCESS_TAX.out.tax,
            params.pr2_version,
            params.pr2_label
        )

        GENERATE_MSCLUSTER(
            dummy_fasta,
            UNCOMPRESS_FASTA.out.uncmp_file,
            PR2_PROCESS_TAX.out.tax,
            params.pr2_version,
            params.pr2_label
        )

    emit:
        fasta = UNCOMPRESS_FASTA.out.uncmp_file
        tax = PR2_PROCESS_TAX.out.tax
        otu = MAKE_OTU_FILE.out.otu

}