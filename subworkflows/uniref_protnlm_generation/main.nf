include { PROTNLM_GETIDS                                } from '../../modules/local/protnlm_getids/main.nf'
include { UNIREF_PROTNLM_FILTER                         } from '../../modules/local/uniref_protnlm_filter/main.nf'
include { DIAMOND_MAKEDB as DIAMOND_MAKEDB_UNIREF_PROTNLM } from '../../modules/nf-core/diamond/makedb/main.nf'

// UniRef90/UniRef50 DIAMOND dbs without the clusters whose representative name was predicted by Google ProtNLM
workflow UNIREF_PROTNLM_GENERATION {
    take:
        trembl_dat
        uniref90_fasta
        uniref50_fasta

    main:
        PROTNLM_GETIDS(trembl_dat)

        channel.of([[id: "uniref90"], uniref90_fasta], [[id: "uniref50"], uniref50_fasta])
            .set { uniref_ch }

        UNIREF_PROTNLM_FILTER(uniref_ch, PROTNLM_GETIDS.out.protnlm_ids.first())

        DIAMOND_MAKEDB_UNIREF_PROTNLM(UNIREF_PROTNLM_FILTER.out.filtered_fasta, [], [], [])

    emit:
        protnlm_ids = PROTNLM_GETIDS.out.protnlm_ids
        stats       = UNIREF_PROTNLM_FILTER.out.stats
        db          = DIAMOND_MAKEDB_UNIREF_PROTNLM.out.db
}
