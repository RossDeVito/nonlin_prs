version 1.0

workflow prs_prsice2 {
    input {
        File geno_bed_file_val
        File geno_bim_file_val
        File geno_fam_file_val
        File geno_bed_file_all
        File geno_bim_file_all
        File geno_fam_file_all
        File sum_stats_file
        File pheno_file
        File covar_file
        File keep_file
        File pred_file
    }

    call prs_prsice2_task {
        input:
            geno_bed_file_val = geno_bed_file_val,
            geno_bim_file_val = geno_bim_file_val,
            geno_fam_file_val = geno_fam_file_val,
            geno_bed_file_all = geno_bed_file_all,
            geno_bim_file_all = geno_bim_file_all,
            geno_fam_file_all = geno_fam_file_all,
            sum_stats_file = sum_stats_file,
            pheno_file = pheno_file,
            covar_file = covar_file,
            keep_file = keep_file,
            pred_file = pred_file
    }

    output {
        Array[File] prsice2_logging_output = prs_prsice2_task.logging_output
        File prsice2_score_output = prs_prsice2_task.score_output
        File prsice2_runtime_json = prs_prsice2_task.runtime_json
        File val_preds = prs_prsice2_task.val_preds
        File test_preds = prs_prsice2_task.test_preds
        File prsice2_clumps = prs_prsice2_task.clumps
        File prsice2_best_thresh = prs_prsice2_task.best_thresh
    }

    meta {
        description: "Run PRSice2 C+T PRS"
    }
}

task prs_prsice2_task {
    input {
        File geno_bed_file_val
        File geno_bim_file_val
        File geno_fam_file_val
        File geno_bed_file_all
        File geno_bim_file_all
        File geno_fam_file_all
        File sum_stats_file
        File pheno_file
        File covar_file
        File keep_file
        File pred_file
    }

    command <<<
        # Get common prefix for BED files for use with PRSice2
        BED_DIR_VAL=$(dirname ~{geno_bed_file_val})
        BED_FNAME_VAL=$(basename ~{geno_bed_file_val} .bed)
        BED_PREFIX_VAL=${BED_DIR_VAL}/${BED_FNAME_VAL}

        BIM_DIR_VAL=$(dirname ~{geno_bim_file_val})
        BIM_FNAME_VAL=$(basename ~{geno_bim_file_val} .bim)
        BIM_PREFIX_VAL=${BIM_DIR_VAL}/${BIM_FNAME_VAL}

        FAM_DIR_VAL=$(dirname ~{geno_fam_file_val})
        FAM_FNAME_VAL=$(basename ~{geno_fam_file_val} .fam)
        FAM_PREFIX_VAL=${FAM_DIR_VAL}/${FAM_FNAME_VAL}

        echo "BED_PREFIX_VAL: $BED_PREFIX_VAL"
        echo "BIM_PREFIX_VAL: $BIM_PREFIX_VAL"
        echo "FAM_PREFIX_VAL: $FAM_PREFIX_VAL"

        # Assert that all files have the same prefix
        if [ "$BED_PREFIX_VAL" != "$BIM_PREFIX_VAL" ] || [ "$BED_PREFIX_VAL" != "$FAM_PREFIX_VAL" ]; then
            echo "BED, BIM, and FAM files used to fit PRSice-2 must have the same prefix"
            exit 1
        fi

        # Set number of threads
        N_THREADS=$(lscpu | grep "^CPU(s):" | awk '{print $2}')

        # Run and time PRSice2
        START_TIME=$(date +%s)

        PRSice_linux \
            --base ~{sum_stats_file} \
            --A1 A1 \
            --beta \
            --bp POS \
            --chr "#CHROM" \
            --snp ID \
            --pvalue P \
            --target ${BED_PREFIX_VAL} \
            --nonfounders \
            --keep ~{keep_file} \
            --ignore-fid \
            --pheno ~{pheno_file} \
            --cov ~{covar_file} \
            --thread ${N_THREADS} \
            --out prs_prsice2

        END_TIME=$(date +%s)
        PRSICE_RUNTIME=$((END_TIME-START_TIME))

        # Get best p-value threshold from PRSice2 output and make input
        # file for plink2 scoring
        BEST_P_THRESH=$(tail -n 1 prs_prsice2.summary | awk '{print $3}')
        echo "best_p 0.0 $BEST_P_THRESH" > prsice2_best_p_thresh.txt

        # Use plink2 to clump variants using parameters used by PRSice2
        # and the validation set only variants
        plink2 \
            --bfile ${BED_PREFIX_VAL} \
            --clump-p1 1 \
            --clump-r2 0.1 \
            --clump-kb 250 \
            --clump ~{sum_stats_file} \
            --out prs_prsice2_clump

        # Score all samples with plink2 using the clumped variants and
        # the best p-value threshold from PRSice2
        plink2 \
            --bed ~{geno_bed_file_all} \
            --bim ~{geno_bim_file_all} \
            --fam ~{geno_fam_file_all} \
            --score ~{sum_stats_file} 3 7 12 header \
            --extract prs_prsice2_clump.clumps \
            --q-score-range prsice2_best_p_thresh.txt \
                ~{sum_stats_file} 3 15 header \
            --out prs_prsice2_score

        # Fit wrapper with python3 script. Saves runtime
        CURRENT_DIR=$(pwd)

        START_TIME=$(date +%s)

        python3 /home/fit_wrapper.py \
            --pheno-file ~{pheno_file} \
            --covar-file ~{covar_file} \
            --score-file prs_prsice2_score.best_p.sscore \
            --val-iids ~{keep_file} \
            --test-iids ~{pred_file} \
            --out-dir $CURRENT_DIR

        END_TIME=$(date +%s)
        FIT_WRAPPER_RUNTIME=$((END_TIME-START_TIME))

        # Save runtime as JSON with total time as 'runtime_seconds' key,
        # and individual times as 'prsice_runtime_seconds' and
        # 'wrapper_runtime_seconds' keys
        TOTAL_RUNTIME=$((PRSICE_RUNTIME+FIT_WRAPPER_RUNTIME))
        echo "{\"runtime_seconds\": $TOTAL_RUNTIME, \"prsice_runtime_seconds\": $PRSICE_RUNTIME, \"wrapper_runtime_seconds\": $FIT_WRAPPER_RUNTIME}" > runtime.json
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/nonlin_prs_prs_prsice2:latest"
    }

    output {
        Array[File] logging_output = glob("*.log")
        File score_output = "prs_prsice2_score.best_p.sscore"
        File val_preds = "val_preds.csv"
        File test_preds = "test_preds.csv"
        File runtime_json = "runtime.json"
        File clumps = "prs_prsice2_clump.clumps"
        File best_thresh = "prsice2_best_p_thresh.txt"
    }
}