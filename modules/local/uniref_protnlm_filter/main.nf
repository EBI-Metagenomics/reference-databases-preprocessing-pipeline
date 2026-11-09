process UNIREF_PROTNLM_FILTER {
    tag "${meta.id}"
    label 'process_single'
    label 'process_long'

    input:
    tuple val(meta), path(uniref_fasta)
    path protnlm_ids

    output:
    tuple val(meta), path("${meta.id}_protnlm_filtered.fasta"), emit: filtered_fasta
    path "protnlm_filter_stats.txt"                           , emit: stats
    path "VERSION.txt"                                        , emit: version_file
    path "versions.yml"                                       , emit: versions

    script:
    """
    set -o pipefail

    # Drop clusters whose representative name was predicted by Google ProtNLM.
    # The cluster name is the representative's protein name, so checking RepID is enough.
    gzip -cdf ${uniref_fasta} \\
        | awk 'FILENAME==ARGV[1]{bad[\$1]; next}
               /^>/{n++; match(\$0, /RepID=[^ ]+/); keep = !(substr(\$0, RSTART+6, RLENGTH-6) in bad); if (!keep) d++}
               keep
               END{printf "clusters_total\\t%d\\nclusters_removed\\t%d\\nclusters_kept\\t%d\\n", n, d, n-d > "protnlm_filter_stats.txt"}' \\
            ${protnlm_ids} - > ${meta.id}_protnlm_filtered.fasta

    printf "trembl_protnlm_entries\\t%d\\n" \$(wc -l < ${protnlm_ids}) >> protnlm_filter_stats.txt

    echo "${params.uniref_protnlm_version}" > VERSION.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${meta.id}: ${params.uniref_protnlm_version}
    END_VERSIONS
    """
}
