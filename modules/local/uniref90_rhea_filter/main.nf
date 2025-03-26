process UNIREF90_RHEA_FILTER {
    label 'process_single'
    container 'quay.io/biocontainers/biopython:1.78'

    input:
    path uniref90_fasta
    path uniprot_rhea_mapping

    output:
    path 'uniref90_with_rhea.fasta', emit: filtered_proteins
    path "versions.yml"            , emit: versions

    script:
    """
    uniref90_rhea_filter.py ${uniref90_fasta} ${uniprot_rhea_mapping} uniref90_with_rhea.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
        uniref90_db: ${params.uniref90_version}
        uniprotKB_access_date: ${params.uniprotKB_access_date}
    END_VERSIONS
    """
}