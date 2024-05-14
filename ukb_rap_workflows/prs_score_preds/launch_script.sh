# MODEL_TYPE=prsice
MODEL_TYPE=basil_lasso

# MODEL_TYPE=automl_max35000_v2_lgbm_v0h24
# MODEL_TYPE=automl_max35000_v2_npart_elastic_net_v0h24

PHENO=standing_height_50
# PHENO=body_fat_percentage_23099
# PHENO=platelet_count_30080
# PHENO=glycated_haemoglobin_30750

WB_ONLY=false

if [ "$WB_ONLY" = true ]; then
	WB_ONLY_FLAG="--wb"
else
	WB_ONLY_FLAG=""
fi

# Launch GWAS workflow
python launcher.py -m ${MODEL_TYPE} -p ${PHENO} ${WB_ONLY_FLAG}