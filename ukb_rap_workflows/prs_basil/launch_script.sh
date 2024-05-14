# PHENO=standing_height_50
PHENO=body_fat_percentage_23099		# n=45, n=49 wb
# PHENO=platelet_count_30080
# PHENO=glycated_haemoglobin_30750

WB_ONLY=false

MODEL_TYPE=lasso
# MODEL_TYPE=ridge
# MODEL_TYPE=elastic_net_0_1
# MODEL_TYPE=elastic_net_0_5
# MODEL_TYPE=elastic_net_0_9

if [ "$WB_ONLY" = true ]; then
	WB_ONLY_FLAG="--wb"
else
	WB_ONLY_FLAG=""
fi
# Launch BASIL workflow
python launcher.py -p ${PHENO} -m $MODEL_TYPE ${WB_ONLY_FLAG} -n 45