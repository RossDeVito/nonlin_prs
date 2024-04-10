# PHENO=standing_height_50
PHENO=body_fat_percentage_23099
# PHENO=platelet_count_30080
# PHENO=glycated_haemoglobin_30750

WB_ONLY=true
DATA_DESC=max35000_v2

# MODEL_CONFIG=lgbm_v0h24_p1e-8_w10000
# MODEL_CONFIG=lgbm_v0h24_p1e-11_w250000
MODEL_CONFIG=elastic_net_v0h24_p1e-11_w250000
# MODEL_CONFIG=npart_elastic_net_v0h2_p1e-8_w0

if [ "$WB_ONLY" = true ]; then
	WB_ONLY_FLAG="--wb"
else
	WB_ONLY_FLAG=""
fi

# Launch AutoML-PRS workflow
python launcher.py \
	-p ${PHENO} ${WB_ONLY_FLAG} \
	-d ${DATA_DESC} \
	-m ${MODEL_CONFIG}