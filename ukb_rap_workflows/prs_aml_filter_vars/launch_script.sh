# PHENO=standing_height_50
# PHENO=body_fat_percentage_23099
# PHENO=platelet_count_30080
PHENO=glycated_haemoglobin_30750

WB_ONLY=false

if [ "$WB_ONLY" = true ]; then
	WB_ONLY_FLAG="--wb"
else
	WB_ONLY_FLAG=""
fi

# Launch GWAS workflow
python launcher.py \
	-p ${PHENO} ${WB_ONLY_FLAG} \
	--max-variants 30000 \
	--out-desc _v2

# _v0 	-p 1e-5 1e-8 1e-11 1e-14 1e-17 
#		-w 0 5000 25000 100000

# _v1 	-p 1e-5 1e-8 1e-11 1e-14 1e-17 1e-24
#		-w 0 5000 10000 25000 100000 250000

# _v2 	-p 1e-5 1e-8 1e-11 1e-14 1e-17 1e-24
#		-w 0 2500 5000 10000 25000 100000 250000