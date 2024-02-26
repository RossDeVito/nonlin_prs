version 1.0

workflow gwas_plink2 {
    input {
        File geno_pgen_file
        File covar_file
        File pheno_file
        File split_file
    }

    call gwas_plink2_task {
        input:
            geno_pgen_file = geno_pgen_file,
            covar_file = covar_file,
            pheno_file = pheno_file,
            split_file = split_file
    }

    output {
        Array[File] glm_linear_output = gwas_plink2_task.glm_linear_output
        File runtime_json = gwas_plink2_task.runtime_json
    }

    meta {
        description: "Run GWAS using plink2"
    }
}

task gwas_plink2_task {
    input {
        File geno_pgen_file
        File covar_file
        File pheno_file
        File split_file
    }

    command <<<
        # Time GWAS
        START_TIME=$(date +%s)

        plink2 --glm 'hide-covar' \
            --pfile ~{geno_pgen_file} \
            --keep ~{split_file} \
            --covar ~{covar_file} \
            --pheno ~{pheno_file} \
            --out gwas_plink2

        # Save runtime as JSON as 'runtime_seconds' key
        END_TIME=$(date +%s)
        ELAPSED_TIME=$((END_TIME-START_TIME))
        echo "{\"runtime_seconds\": $ELAPSED_TIME}" > runtime.json
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/nonlin_prs_gwas_plink2:latest"
    }

    output {
        Array[File] glm_linear_output = glob("gwas_plink2.*.glm.linear")
        File runtime_json = "runtime.json"
    }
}