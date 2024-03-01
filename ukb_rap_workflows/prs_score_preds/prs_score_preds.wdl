version 1.0

workflow prs_score_preds {
	input {
		File val_preds
		File test_preds
		File pheno_file
		File test_wb_samples
	}

	call score_preds {
		input:
			val_preds = val_preds,
			test_preds = test_preds,
			pheno_file = pheno_file,
			test_wb_samples = test_wb_samples
	}

	output {
		File scores_json = score_preds.scores_json
		Array[File] plots = score_preds.plots
	}

	meta {
		description: "Score and plot a model's predictions."
	}
}

task score_preds {
	input {
		File val_preds
		File test_preds
		File pheno_file
		File test_wb_samples
	}

	command <<<
		CURRENT_DIR=$(pwd)

		python3 /home/score_preds.py \
			--val-preds ~{val_preds} \
			--test-preds ~{test_preds} \
			--pheno-file ~{pheno_file} \
			--wb ~{test_wb_samples} \
			--out-dir $CURRENT_DIR
	>>>

	runtime {
		docker: "gcr.io/ucsd-medicine-cast/nonlin_prs_prs_score_preds:latest"
	}

	output {
		File scores_json = "scores.json"
		Array[File] plots = glob("*.png")
	}
}