version 1.0

workflow prs_aml_filter_vars {
    input {
        File sum_stats_file
        File geno_pgen_file
        File geno_psam_file
        File geno_pvar_file
        Int max_num_vars
    }

    call prs_aml_filter_vars_task {
        input:
            sum_stats_file = sum_stats_file,
            geno_pgen_file = geno_pgen_file,
            geno_psam_file = geno_psam_file,
            geno_pvar_file = geno_pvar_file,
            max_num_vars = max_num_vars
    }

    output {
        File dosage_table = prs_aml_filter_vars_task.dosage_table
        File filtered_vars_json = prs_aml_filter_vars_task.filtered_vars_json
        File filtered_vars_meta = prs_aml_filter_vars_task.filtered_vars_meta
    }

    meta {
        description: "Preprocess genotypes for AutoML_PRS"
    }
}

task prs_aml_filter_vars_task {
    input {
        File sum_stats_file
        File geno_pgen_file
        File geno_psam_file
        File geno_pvar_file
        Int max_num_vars
    }

    command <<<
        echo "Filtering variants by p-value and window size"

        python3 /home/AutoML_PRS/data_preprocessing/filter_vars_by_pval.py \
            --sum-stats-file ~{sum_stats_file} \
            -p 1e-5 1e-8 1e-11 1e-14 1e-17 1e-24 \
            -w 0 2500 5000 10000 25000 100000 250000 \ 
            -m ~{max_num_vars}

        echo "Filtering variants by p-value and window size complete"

        echo "Filter and convert to Python readable format"

        plink2 \
            --pgen ~{geno_pgen_file} \
            --psam ~{geno_psam_file} \
            --pvar ~{geno_pvar_file} \
            --extract filtered_vars_all.txt \
            --export A \
            --out filtered_doseage_table

        python3 /home/AutoML_PRS/data_preprocessing/raw_to_input_parquet.py \
            -f filtered_vars_raw.json \
            -r filtered_doseage_table.raw
        >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/nonlin_prs_prs_aml_filter_vars:latest"
    }

    output {
        File dosage_table = "filtered_vars.parquet"
        File filtered_vars_json = "filtered_vars.json"
        File filtered_vars_meta = "filtered_vars_meta.json"
    }
}