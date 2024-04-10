version 1.0

workflow prs_aml {
    input {
        File geno_parquet
        File var_subset_json
        File pheno_file
        File covar_file
        File train_ids
        File val_ids
        File test_ids
        String train_config_path
    }

    call prs_aml_task {
        input:
            geno_parquet = geno_parquet,
            var_subset_json = var_subset_json,
            pheno_file = pheno_file,
            covar_file = covar_file,
            train_ids = train_ids,
            val_ids = val_ids,
            test_ids = test_ids
            train_config_path = train_config_path
    }

    output {
        File best_model_config = prs_aml_task.best_model_config
        File best_model = prs_aml_task.best_model
        File fit_log = prs_aml_task.fit_log
        File training_config = prs_aml_task.training_config
        File learning_curve = prs_aml_task.learning_curve
        File val_preds = prs_aml_task.val_preds
        File test_preds = prs_aml_task.test_preds
        File runtime_json = prs_aml_task.runtime_json
    }

    meta {
        description: "Fit AutoML_PRS"
    }
}

task prs_aml_task {
    input {
        File geno_parquet
        File var_subset_json
        File pheno_file
        File covar_file
        File train_ids
        File val_ids
        File test_ids
        String train_config_path
    }

    command <<<
        fit_automl_prs \
            --training-config ~{train_config_path} \
            --geno-parquet ~{geno_parquet} \
            --var-subsets ~{var_subset_json} \
            --pheno ~{pheno_file} \
            --covars ~{covar_file} \
            --train-ids ~{train_ids} \
            --val-ids ~{val_ids} \
            --test-ids ~{test_ids} 
            
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/nonlin_prs_prs_aml:latest"
    }

    output {
        File best_model_config = "best_model_config.json"
        File best_model = "best_model.pkl"
        File fit_log = "fit.log"
        File training_config = "training_config.json"
        File learning_curve = "learning_curve.png"
        File val_preds = "val_preds.csv"
        File test_preds = "test_preds.csv"
        File runtime_json = "runtime.json"
    }
}