# PHENO=standing_height_50
# PHENO=body_fat_percentage_23099
PHENO=platelet_count_30080
# PHENO=glycated_haemoglobin_30750

WB_ONLY=true

if [ "$WB_ONLY" = true ]; then
	WB_ONLY_FLAG="--wb"
else
	WB_ONLY_FLAG=""
fi

# Launch GWAS workflow
python launcher.py -p ${PHENO} ${WB_ONLY_FLAG}