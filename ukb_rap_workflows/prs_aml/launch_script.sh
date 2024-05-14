PHENO=standing_height_50
# PHENO=body_fat_percentage_23099
# PHENO=platelet_count_30080
# PHENO=glycated_haemoglobin_30750

WB_ONLY=false

CLUMPS=false
BASIL_LASSO=true

MODEL_CONFIG=lgbm_v0h24
# MODEL_CONFIG=elastic_net_v0h24
# MODEL_CONFIG=npart_elastic_net_v0h24
# MODEL_CONFIG=sgd_elastic_net_v0h24

if [ "$CLUMPS" = true ]; then
	DATA_DESC=clumps
	MODEL_CONFIG="${MODEL_CONFIG}_clumps"
	INSTANCE_FLAG="-l"
elif [ "$BASIL_LASSO" = true ]; then
	DATA_DESC=basil_lasso
	MODEL_CONFIG="${MODEL_CONFIG}_basil"
	INSTANCE_FLAG="-l"
else
	DATA_DESC=max30000_v2
	INSTANCE_FLAG=""
fi

if [ "$WB_ONLY" = true ]; then
	WB_ONLY_FLAG="--wb"
else
	WB_ONLY_FLAG=""
fi

# Launch AutoML-PRS workflow
python launcher.py \
	-p ${PHENO} ${WB_ONLY_FLAG} ${INSTANCE_FLAG} \
	-d ${DATA_DESC} \
	-m ${MODEL_CONFIG}