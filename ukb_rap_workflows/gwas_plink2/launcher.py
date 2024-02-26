"""Launch plink2 GWAS workflow on the UKB RAP.

Uses PGEN input genotype file.

Required args:

* -p, --pheno-name: Name of the phenotype to use for the GWAS.

Optional args:

* --wb: Flag to use just the white British subset of the UKB data.
	False when not provided.
* --geno-dir: Directory containing the genotype files. Default: 
	'/project/rdevito/nonlin_prs/data/geno_data/qced_common/pgen'
* --splits-dir: Directory containing train/val/test splits in
	the form of list of sample IDs. Default: 
	'/rdevito/nonlin_prs/data/sample_data/splits'
* --default_covar_

"""

import argparse
import dxpy
import sys

WORKFLOW_ID = 'workflow-GgKVPXjJv7B4G4gXXfK6kKjG'