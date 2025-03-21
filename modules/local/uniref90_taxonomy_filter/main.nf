process UNIREF90_TAXONOMY_FILTER {
    label 'process_single'
    container 'community.wave.seqera.io/library/pip_biopython_pyfastx_taxoniq:c440284a61b91ed0'

    input:
    path uniref90_fasta

    output:
    path '*.filtered.fasta', emit: filtered_proteins
    path '*.protid2taxid'  , emit: protid2taxid
    path 'versions.yml'    , emit: versions

    script:
    """
    uniref90_taxonomy_filter.py ${uniref90_fasta} uniref90

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
        taxoniq: \$(python -c "import importlib.metadata; print(importlib.metadata.version('taxoniq'))")
        pyfastx: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pyfastx'))")
        uniref90_db: ${params.uniref90_version}
    END_VERSIONS
    """
}