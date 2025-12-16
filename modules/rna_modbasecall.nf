// Basecalling for m6A_DRACH modification
process RNA_MOD_M6A_DRACH {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
    containerOptions '--gpus all'
    tag "${sample_id}"
    publishDir "${params.output_dir}/basecalling_rna", mode: 'copy'

    input:
    tuple val(sample_id), path(pod5)

    output:
    tuple val(sample_id), path("${sample_id}_q${params.min_qscore}_m6a.bam")

    script:
    """
    dorado download --model ${params.model_sup}
    dorado basecaller ${params.model_sup} ${pod5} --min-qscore ${params.min_qscore} --modified-bases m6A_DRACH > ${sample_id}_q${params.min_qscore}_m6a.bam
    """
}

// Basecalling for m5C modification
process RNA_MOD_M5C {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
    containerOptions '--gpus all'
    tag "${sample_id}"
    publishDir "${params.output_dir}/basecalling_rna", mode: 'copy'

    input:
    tuple val(sample_id), path(pod5)

    output:
    tuple val(sample_id), path("${sample_id}_q${params.min_qscore}_m5c.bam")

    script:
    """
    dorado download --model ${params.model_hac}
    dorado basecaller ${params.model_hac} ${pod5} --min-qscore ${params.min_qscore} --modified-bases m5C > ${sample_id}_q${params.min_qscore}_m5c.bam
    """
}

// Basecalling for pseU modification
process RNA_MOD_PSEU {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
    containerOptions '--gpus all'
    tag "${sample_id}"
    publishDir "${params.output_dir}/basecalling_rna", mode: 'copy'

    input:
    tuple val(sample_id), path(pod5)

    output:
    tuple val(sample_id), path("${sample_id}_q${params.min_qscore}_pseu.bam")

    script:
    """
    dorado download --model ${params.model_hac}
    dorado basecaller --min-qscore ${params.min_qscore} --modified-bases pseU ${params.model_hac} ${pod5} > ${sample_id}_q${params.min_qscore}_pseu.bam
    """
}

// Basecalling for inosine_m6A modification
process RNA_MOD_INOSINE {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
    containerOptions '--gpus all'
    tag "${sample_id}"
    publishDir "${params.output_dir}/basecalling_rna", mode: 'copy'

    input:
    tuple val(sample_id), path(pod5)

    output:
    tuple val(sample_id), path("${sample_id}_q${params.min_qscore}_inosine.bam")

    script:
    """
    dorado download --model ${params.model_hac}
    dorado basecaller ${params.model_hac} ${pod5} --min-qscore ${params.min_qscore} --modified-bases inosine_m6A > ${sample_id}_q${params.min_qscore}_inosine.bam
    """
}

// Basecalling for m5C_2OmeC modification
process RNA_MOD_C_2OME {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
    containerOptions '--gpus all'
    tag "${sample_id}"
    publishDir "${params.output_dir}/basecalling_rna", mode: 'copy'

    input:
    tuple val(sample_id), path(pod5)

    output:
    tuple val(sample_id), path("${sample_id}_q${params.min_qscore}_m5C_2OmeC.bam")

    script:
    """
    dorado download --model ${params.model_sup}
    dorado basecaller ${params.model_sup} ${pod5} --min-qscore ${params.min_qscore} --modified-bases m5C_2OmeC > ${sample_id}_q${params.min_qscore}_m5C_2OmeC.bam
    """
}

// Basecalling for inosine_m6A_2OmeA modification
process RNA_MOD_A_2OME {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
    containerOptions '--gpus all'
    tag "${sample_id}"
    publishDir "${params.output_dir}/basecalling_rna", mode: 'copy'

    input:
    tuple val(sample_id), path(pod5)

    output:
    tuple val(sample_id), path("${sample_id}_q${params.min_qscore}_inosine_m6A_2OmeA.bam")

    script:
    """
    dorado download --model ${params.model_sup}
    dorado basecaller ${params.model_sup} ${pod5} --min-qscore ${params.min_qscore} --modified-bases inosine_m6A_2OmeA > ${sample_id}_q${params.min_qscore}_inosine_m6A_2OmeA.bam
    """
}

// Basecalling for pseU_2OmeU modification
process RNA_MOD_U_2OME {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
    containerOptions '--gpus all'
    tag "${sample_id}"
    publishDir "${params.output_dir}/basecalling_rna", mode: 'copy'

    input:
    tuple val(sample_id), path(pod5)

    output:
    tuple val(sample_id), path("${sample_id}_q${params.min_qscore}_pseU_2OmeU.bam")

    script:
    """
    dorado download --model ${params.model_sup}
    dorado basecaller ${params.model_sup} ${pod5} --min-qscore ${params.min_qscore} --modified-bases pseU_2OmeU > ${sample_id}_q${params.min_qscore}_pseU_2OmeU.bam
    """
}

// Basecalling for 2OmeG modification
process RNA_MOD_G_2OME {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
    containerOptions '--gpus all'
    tag "${sample_id}"
    publishDir "${params.output_dir}/basecalling_rna", mode: 'copy'

    input:
    tuple val(sample_id), path(pod5)

    output:
    tuple val(sample_id), path("${sample_id}_q${params.min_qscore}_2OmeG.bam")

    script:
    """
    dorado download --model ${params.model_sup}
    dorado basecaller ${params.model_sup} ${pod5} --min-qscore ${params.min_qscore} --modified-bases 2OmeG > ${sample_id}_q${params.min_qscore}_2OmeG.bam
    """
}

// callall modifications
process RNA_MOD_ALL {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
    containerOptions '--gpus all'
    tag "${sample_id}"
    publishDir "${params.output_dir}/basecalling_rna", mode: 'copy'

    input:
    tuple val(sample_id), path(pod5)

    output:
    tuple val(sample_id), path("${sample_id}_all_mods.bam")

    script:
    def model = params.model_sup ?: "rna004_130bps_sup@v5.2.0"
    def min_qscore = params.min_qscore ?: 10
    // Models targeting all 4 canonical bases: C, A, U, G
    def mod_models = "${model}_m5C_2OmeC@v1,${model}_inosine_m6A_2OmeA@v1,${model}_pseU_2OmeU@v1,${model}_2OmeG@v1"
    """
    set -e
    set -o pipefail

    echo "Downloading model: ${model}"
    dorado download --model ${model}

    echo "Downloading modification models..."
    dorado download --model ${model}_m5C_2OmeC@v1
    dorado download --model ${model}_inosine_m6A_2OmeA@v1
    dorado download --model ${model}_pseU_2OmeU@v1
    dorado download --model ${model}_2OmeG@v1

    echo "Starting basecalling for sample: ${sample_id}"
    echo "POD5 path: ${pod5}"
    echo "Model: ${model}"
    echo "Modification models: ${mod_models}"
    echo "Min Q score: ${min_qscore}"

    dorado basecaller ${model} ${pod5} \
        --min-qscore ${min_qscore} \
        --modified-bases-models ${mod_models} \
        --modified-bases-threshold 0.2 \
        --emit-moves > ${sample_id}_all_mods.bam

    # Check if output file was created and has content
    if [ ! -s "${sample_id}_all_mods.bam" ]; then
        echo "ERROR: Output BAM file is empty or was not created!"
        exit 1
    fi

    ls -lh ${sample_id}_all_mods.bam
    echo "Basecalling completed successfully"
    """
}