version 1.0

workflow prs_prsice2 {
    input {
        File geno_bed_file
        File geno_bim_file
        File geno_fam_file
        File sum_stats_file
        File pheno_file
        File covar_file
        File keep_file
    }

    call prs_prsice2_task {
        input:
            geno_bed_file = geno_bed_file,
            geno_bim_file = geno_bim_file,
            geno_fam_file = geno_fam_file,
            sum_stats_file = sum_stats_file,
            pheno_file = pheno_file,
            covar_file = covar_file,
            keep_file = keep_file
    }

    output {
        Array[File] prsice2_output = prs_prsice2_task.prsice2_output
        File runtime_json = prs_prsice2_task.runtime_json
    }

    meta {
        description: "Run PRSice2 C+T PRS"
    }
}

task prs_prsice2_task {
    input {
        File geno_bed_file
        File geno_bim_file
        File geno_fam_file
        File sum_stats_file
        File pheno_file
        File covar_file
        File keep_file
    }

    command <<<
        # Get common prefix for BED files
        BED_PREFIX=$(basename ~{geno_bed_file} .bed)
        BIM_PREFIX=$(basename ~{geno_bim_file} .bim)
        FAM_PREFIX=$(basename ~{geno_fam_file} .fam)

        echo "BED_PREFIX: $BED_PREFIX"
        echo "BIM_PREFIX: $BIM_PREFIX"
        echo "FAM_PREFIX: $FAM_PREFIX"

        # Assert that all files have the same prefix
        if [ "$BED_PREFIX" != "$BIM_PREFIX" ] || [ "$BED_PREFIX" != "$FAM_PREFIX" ]; then
            echo "BED, BIM, and FAM files must have the same prefix"
            exit 1
        fi

        # Set number of threads
        N_THREADS=$(lscpu | grep "^CPU(s):" | awk '{print $2}')

        # Time PRSice2
        START_TIME=$(date +%s)

        PRSice_linux \
            --base ~{sum_stats_file} \
            --A1 A1 \
            --beta \
            --bp POS \
            --chr "#CHROM" \
            --snp ID \
            --pvalue P \
            --target ${BED_PREFIX} \
            --nonfounders \
            --keep ~{keep_file} \
            --ignore-fid \
            --pheno ~{pheno_file} \
            --cov ~{covar_file} \
            --thread ${N_THREADS} \
            --out prs_prsice2

        # Save runtime as JSON as 'runtime_seconds' key
        END_TIME=$(date +%s)
        ELAPSED_TIME=$((END_TIME-START_TIME))
        echo "{\"runtime_seconds\": $ELAPSED_TIME}" > runtime.json
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/nonlin_prs_prs_prsice2:latest"
    }

    output {
        Array[File] prsice2_output = glob("prs_prsice2*")
        File runtime_json = "runtime.json"
    }
}