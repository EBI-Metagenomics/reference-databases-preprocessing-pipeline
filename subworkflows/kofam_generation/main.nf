include { ADD_METADATA_TO_HMM } from '../../modules/local/add_metadata_to_hmm/main.nf'
include { HMMER_HMMPRESS      } from '../../modules/nf-core/hmmer/hmmpress/main.nf'

workflow KOFAM_GENERATION {
    take:
    ko_hmm_dir
    ko_list

    main:
    ADD_METADATA_TO_HMM(ko_hmm_dir, ko_list)

    ADD_METADATA_TO_HMM.out.kofam_modified
        .map { db_file ->
            [[id: "KOFAM_db"], db_file]
        }
        .set { hmmer_input}

    HMMER_HMMPRESS(hmmer_input)
}