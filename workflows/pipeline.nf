
// nextflow run -resume -profile slurm main.nf

include { SILVA_GENERATION as SILVA_SSU_GENERATION } from '../subworkflows/silva_generation/main.nf'
include { SILVA_GENERATION as SILVA_LSU_GENERATION } from '../subworkflows/silva_generation/main.nf'
include { PR2_GENERATION                           } from "../subworkflows/pr2_generation/main.nf"
include { UNITE_GENERATION                         } from "../subworkflows/unite_generation/main.nf"
include { ITSONEBD_GENERATION                      } from "../subworkflows/itsonedb_generation/main.nf"
include { RFAM_GENERATION                          } from "../subworkflows/rfam_generation/main.nf"
include { RHEA_AND_TAXONOMY_GENERATION             } from "../subworkflows/rhea_and_taxonomy_generation/main.nf"
include { KOFAM_GENERATION                         } from "../subworkflows/kofam_generation/main.nf"
include { UNIREF_PROTNLM_GENERATION                } from "../subworkflows/uniref_protnlm_generation/main.nf"

include { METADATA_LOG                             } from '../modules/local/metadata_log/main.nf'

workflow REFERENCE_DATABASES_PREPROCESSING {
    if (params.generate_amplicon_db) {
        silva_taxdump = file(params.silva_download_taxdump, checkIfExists: true)

        SILVA_SSU_GENERATION(
            "SSU",
            silva_taxdump,
            params.silva_ssu_download_fasta,
            params.silva_ssu_label
        )

        SILVA_LSU_GENERATION(
            "LSU",
            silva_taxdump,
            params.silva_lsu_download_fasta,
            params.silva_lsu_label
        )
        
        PR2_GENERATION()

        UNITE_GENERATION()

        ITSONEBD_GENERATION()

        RFAM_GENERATION()

        METADATA_LOG()
    }
    
    if (params.generate_rhea_tax_db) {
        RHEA_AND_TAXONOMY_GENERATION(
            file(params.uniref90_fasta_link, checkIfExists: true),
            file(params.uniprot_rhea_mapping, checkIfExists: true),
            file(params.rhea_chebi_mapping_link, checkIfExists: true),
            file(params.ncbi_taxdump_link, checkIfExists: true)
        )
    }

    if (params.generate_kofam_db) {
        KOFAM_GENERATION(
            file(params.kofam_ko_hmm_link, checkIfExists: true),
            file(params.kofam_ko_list_link, checkIfExists: true),
        )
    }

    if (params.generate_uniref_protnlm_db) {
        uniprot_links = [params.trembl_dat_link, params.uniref90_fasta_link, params.uniref50_fasta_link]

        // UniProt only serves the current release, so check it is the one we are labelling the dbs with
        uniprot_links.each { link ->
            def release = file(link).parent.resolve('RELEASE.metalink').text.find(/<version>([^<]+)</) { _match, version -> version }
            if (release != params.uniref_protnlm_version) {
                error("${link} is UniProt release ${release}, but --uniref_protnlm_version is ${params.uniref_protnlm_version}")
            }
        }

        UNIREF_PROTNLM_GENERATION(
            file(params.trembl_dat_link, checkIfExists: true),
            file(params.uniref90_fasta_link, checkIfExists: true),
            file(params.uniref50_fasta_link, checkIfExists: true)
        )
    }
}