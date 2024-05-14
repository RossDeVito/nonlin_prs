version 1.0

workflow prs_aml_filter_vars {
    input {
        File clumps_file
        File best_pval_file
        File geno_pgen_file
        File geno_psam_file
        File geno_pvar_file
    }

    call prs_aml_filter_vars_task {
        input:
            clumps_file = clumps_file,
            best_pval_file = best_pval_file,
            geno_pgen_file = geno_pgen_file,
            geno_psam_file = geno_psam_file,
            geno_pvar_file = geno_pvar_file,
    }

    output {
        File dosage_table = prs_aml_filter_vars_task.dosage_table
        File filtered_vars_json = prs_aml_filter_vars_task.filtered_vars_json
        File filtered_vars_meta = prs_aml_filter_vars_task.filtered_vars_meta
    }

    meta {
        description: "Preprocess genotypes for AutoML_PRS using clump leads"
    }
}

task prs_aml_filter_vars_task {
    input {
        File clumps_file
        File best_pval_file
        File geno_pgen_file
        File geno_psam_file
        File geno_pvar_file
    }

    command <<<
        echo "Filtering variants to clump leads"

        python3 /home/AutoML_PRS/data_preprocessing/filter_vars_by_clump.py \
            --clumps-file ~{clumps_file} \
            --pval-thresh-file ~{best_pval_file}

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