version 1.0

workflow numbat_workflow {
  input {
    String output_directory

    # pileup/phase parameters
    String label
    String sample_name
    File bam_file
    File bam_index_file
    File barcodes_file
    File count_matrix
    File ref_count_matrix
    File ref_annotations

    # parameters
    String transition = '1e-5'
    String gamma = '20'
    String min_cells = '50'
    String max_iter = '2'
    String min_LLR = '5'
    String max_entropy = '0.5'
    String init_k = '3'
    String tau = '0.3'
    String plot = 'TRUE'

    #general parameters
    Int preemptible = 2
    Int disk_space = 150
    Int num_cpu = 16
    Int memory = 128
    String docker = "mparikhbroad/numbat:latest"
  }

  call pileup_and_phase {
    input:
        output_directory = sub(output_directory, "/+$", ""),
        label = label,
        sample_name = sample_name,
        bam_file = bam_file,
        bam_index_file = bam_index_file,
        barcodes_file = barcodes_file,
        preemptible = preemptible,
        disk_space = disk_space,
        docker = docker,
        num_cpu = num_cpu,
        memory = memory,
  }

  call run_numbat {
    input:
      label = label,
      output_directory = sub(output_directory, "/+$", ""),
      allele_counts = pileup_and_phase.allele_counts,
      count_matrix = count_matrix,
      ref_count_matrix = ref_count_matrix,
      ref_annotations = ref_annotations,

      transition = transition,
      gamma = gamma,
      min_cells = min_cells,
      max_iter = max_iter,
      min_LLR = min_LLR,
      max_entropy = max_entropy,
      init_k = init_k,
      tau = tau,
      plot = plot,

      preemptible = preemptible,
      disk_space = disk_space,
      docker = docker,
      num_cpu = num_cpu,
      memory = memory,
  }

  output {
    File pileup_and_phase_output = pileup_and_phase.allele_counts
    File numbat_output = run_numbat.output_zip
  }
}

task pileup_and_phase {
  input {
    String output_directory
    String label
    String sample_name
    File bam_file
    File bam_index_file
    File barcodes_file
    Int preemptible
    Int disk_space
    Int num_cpu
    Int memory
    Int? memory_override
    String docker
  }

  command <<<

    set -e

    mkdir -p 'inputs/~{label}'
    mkdir -p 'outputs/~{label}'

    R --no-save  <<RSCRIPT

    library(glue)

    label = "~{label}"
    sample = "~{sample_name}"

    bam_file = "~{bam_file}"
    bam_index_file  = "~{bam_index_file}"
    barcode_file = "~{barcodes_file}"

    # system(glue('cp {bam_file} inputs/{label}/{sample}.bam'))
    # bam_file <- glue('{getwd()}/inputs/{label}/{sample}.bam')
    # all_bam_files <- c(all_bam_files, bam_file)

    # system(glue('cp {bam_index_file} inputs/{label}/{sample}.bam.bai'))
    # bam_index_file <- glue('inputs/{label}/{sample}.bam.bai')

    if (endsWith(barcode_file[[1]], 'gz')) {
        system(glue('cp {barcode_file} inputs/{label}/{sample}.tsv.gz'))
        system(glue('gunzip inputs/{label}/{sample}.tsv.gz'))
    } else {
        system(glue('cp {barcode_file} inputs/{label}/{sample}.tsv'))
    }
    barcode_file <- glue('{getwd()}/inputs/{label}/{sample}.tsv')

    args = c('--label', label,
             '--samples', sample,
             '--bams', bam_file,
             '--barcodes', barcode_file,
             '--outdir', glue('{getwd()}/outputs/{label}/'), 
             '--gmap', '/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz', 
             '--snpvcf', '/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf', 
             '--paneldir', '/data/1000G_hg38', 
             '--ncores', '~{num_cpu}')


    arg_vals <- paste(args, collapse = ' ')
    print(glue("Rscript /numbat/inst/bin/pileup_and_phase.R {arg_vals}"))

    system(paste0(glue("Rscript /numbat/inst/bin/pileup_and_phase.R {arg_vals}")))

    RSCRIPT

    gsutil -m cp -r "outputs/~{label}/" ~{output_directory}/

  >>>

  output {
    File allele_counts = "outputs/~{label}/~{sample_name}_allele_counts.tsv.gz"
  }

  runtime {
    preemptible: preemptible
    bootDiskSizeGb: 12
    disks: "local-disk ${disk_space} HDD"
    docker: docker
    cpu: num_cpu
    memory: select_first([memory_override, memory]) + "GB"
  }
}

task run_numbat {
  input {
    String label
    String output_directory
    File allele_counts

    File count_matrix
    File ref_count_matrix
    File ref_annotations

    String transition
    String gamma
    String min_cells
    String max_iter
    String min_LLR
    String max_entropy
    String init_k
    String tau
    String plot

    Int preemptible
    Int disk_space
    Int num_cpu
    Int memory
    Int? memory_override
    String docker
  }

  command <<<

    set -e

    mkdir -p 'inputs/~{label}'
    mkdir -p 'outputs/~{label}'

    R --no-save  <<RSCRIPT

    library(glue)
    library(numbat)

    label = '~{label}'
    print(glue('label = {label}'))

    # allele_counts
    allele_counts_file = "~{allele_counts}"
    allele_counts = read.table(allele_counts_file, sep='\t')

    # count_mat
    count_mat_file = "~{count_matrix}"
    count_mat = readRDS(count_mat_file)

    # ref_count_mat
    ref_count_mat_file = "~{ref_count_matrix}"
    ref_count_mat = readRDS(ref_count_mat_file)

    # ref_annot
    ref_annot_file = "~{ref_annotations}"
    ref_annot = read.table(ref_annot_file, sep=",", header=TRUE) #fix this

    ref_internal = aggregate_counts(ref_count_mat, ref_annot, normalized=TRUE)

    out = run_numbat(
      count_mat,
      ref_internal,
      allele_counts,
      t = ~{transition},
      gamma = ~{gamma},
      min_cells = ~{min_cells},
      max_iter = ~{max_iter},
      min_LLR = ~{min_LLR},
      max_entropy = ~{max_entropy},
      init_k = ~{init_k},
      tau = ~{tau},
      plot = ~{plot},
      n_cores = ~{num_cpu}
      out_dir = "outputs/~{label}"
    )

    RSCRIPT

    gsutil -m cp -r "outputs/~{label}/" "~{output_directory}/"

    zip -r outputs.zip outputs/~{label}/*

  >>>

  output {
    File output_zip = "outputs.zip"
  }

  runtime {
    preemptible: preemptible
    bootDiskSizeGb: 12
    disks: "local-disk ${disk_space} HDD"
    docker: docker
    cpu: num_cpu
    memory: select_first([memory_override, memory]) + "GB"
  }
}