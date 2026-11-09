process PROTNLM_GETIDS {
    label 'process_single'
    label 'process_long'

    input:
    path trembl_dat

    output:
    path "protnlm_ids.txt", emit: protnlm_ids
    path "versions.yml"   , emit: versions

    script:
    """
    set -o pipefail

    # IDs of TrEMBL entries whose name was predicted by Google ProtNLM.
    # Only the first top-level name counts, as that is the name UniRef uses for the cluster.
    gzip -cdf ${trembl_dat} \\
        | awk '/^ID   /{id=\$2; named=0}
               /^DE   (RecName|SubName): Full=/ && !named++ && /Google:ProtNLM/{print id}' > protnlm_ids.txt

    # An empty list would silently keep every cluster, e.g. if UniProt changes the evidence string
    [ -s protnlm_ids.txt ] || { echo "No ProtNLM entries found in ${trembl_dat}" >&2; exit 1; }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        uniprot_release: ${params.uniref_protnlm_version}
    END_VERSIONS
    """
}
