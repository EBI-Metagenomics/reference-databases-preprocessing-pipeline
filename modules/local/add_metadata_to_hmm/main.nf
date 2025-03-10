process ADD_METADATA_TO_HMM {
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    path ko_hmm_dir
    path ko_list

    output:
    path "kofam_modified", emit: kofam_modified
    path "versions.yml"  , emit: versions

    script:
    def is_compressed = ko_hmm_dir.getExtension() == "gz" ? true : false
    def ko_hmm_input = is_compressed ? ko_hmm_dir.getName().replace(".tar.gz", "") : ko_hmm_dir

    """
    if [ "${is_compressed}" == "true" ]; then
        tar -xzf ${ko_hmm_dir}
    fi

    add_metadata_to_hmm.py ${ko_list} ${ko_hmm_input}

    cat *.modified.hmm > kofam_modified

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch kofam_modified

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}