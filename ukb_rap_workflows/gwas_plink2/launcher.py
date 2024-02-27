"""Launch plink2 GWAS workflow on the UKB RAP.

Uses PGEN input genotype file.

Example usage:
```
python launcher.py -p platelet_count_30080 --wb
```

Required args:

* -p, --pheno-name: Name of the phenotype to use for the GWAS.

Optional args:

* --wb: Flag to use just the white British subset of the UKB data.
	False when not provided.
* --dev: Flag to only include the development set of the train split.
	False when not provided. Will be reflected in the output directory name.
* --pheno-dir: Directory containing the phenotype files. Default:
	'/rdevito/nonlin_prs/data/pheno_data/pheno'
* --pheno-metadata-file: File containing the phenotype metadata.
	Default: '../../data/pheno_metadata.json'
* --geno-dir: Directory containing the genotype files. Default: 
	'/rdevito/nonlin_prs/data/geno_data/qced_common/pgen'
* --splits-dir: Directory containing train/val/test splits in
	the form of list of sample IDs. Default: 
	'/rdevito/nonlin_prs/data/sample_data/splits'
* --covar-dir: Directory containing the covariate files. Default:
	'/rdevito/nonlin_prs/data/covar_data/tsv'
* --default_covar_set: Default covariate set file name (w/o '.tsv') 
	if no specific set of covariates is provided by the phenotype
	metadata file. If in metadata file, the key should be 'covar_set'.
	Default: 'covar_std_v1'
* --output-dir: Directory in which a folder will be created to store
	the output of the GWAS workflow. Default: 
	'/rdevito/nonlin_prs/gwas/gwas_output'. Final output directory
	will be of the form: {output_dir}/{pheno_name}_glm[_wb][_dev]
"""

import argparse
import json

import dxpy


WORKFLOW_ID = 'workflow-GgKXpk8Jv7BBFbGbPJzp97by'
DEFAULT_INSTANCE = 'mem2_ssd1_v2_x32'


def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'-p', '--pheno-name',
		required=True,
		help='Name of the phenotype to use for the GWAS.'
	)
	parser.add_argument(
		'--wb',
		action='store_true',
		help='Flag to use just the white British subset of the UKB data.'
	)
	parser.add_argument(
		'--dev',
		action='store_true',
		help='Flag to only include the development set of the train split.'
	)
	parser.add_argument(
		'--pheno-dir',
		default='/rdevito/nonlin_prs/data/pheno_data/pheno',
		help='Directory containing the phenotype files.'
	)
	parser.add_argument(
		'--pheno-metadata-file',
		default='../../data/pheno_metadata.json',
		help='File containing the phenotype metadata.'
	)
	parser.add_argument(
		'--geno-dir',
		default='/rdevito/nonlin_prs/data/geno_data/qced_common/pgen',
		help='Directory containing the genotype files.'
	)
	parser.add_argument(
		'--splits-dir',
		default='/rdevito/nonlin_prs/data/sample_data/splits',
		help='Directory containing train/val/test splits in the form of '
			 'list of sample IDs.'
	)
	parser.add_argument(
		'--covar-dir',
		default='/rdevito/nonlin_prs/data/covar_data/tsv',
		help='Directory containing the covariate files.'
	)
	parser.add_argument(
		'--default_covar_set',
		default='covar_std_v1',
		help='Default covariate set file name (w/o \'.tsv\') if no specific '
			 'set of covariates is provided by the phenotype metadata file. '
			 'If in metadata file, the key should be \'covar_set\'.'
	)
	parser.add_argument(
		'--output-dir',
		default='/rdevito/nonlin_prs/gwas/gwas_output',
		help='Directory in which a folder will be created to store the output '
			 'of the GWAS workflow. Final output directory will be of the form: '
			 '{output_dir}/{pheno_name}_glm[_wb][_dev]'
	)
	return parser.parse_args()


def launch_gwas_workflow(
	geno_file,
	covar_file,
	pheno_file,
	split_file,
	output_dir,
	workflow_id=WORKFLOW_ID,
	instance_type=DEFAULT_INSTANCE,
	name='gwas_plink2'
):
	"""Launch the GWAS workflow on the UKB RAP.
	
	Args:
		geno_file: Path to the PGEN genotype file in UKB RAP storage.
			The filename should exclude the .pgen/.psam/.pvar extensions.
		covar_file: Path to the covariate file in UKB RAP storage.
		pheno_file: Path to the phenotype file in UKB RAP storage.
		split_file: Path to the train split file in UKB RAP storage.
		output_dir: Path to the output directory in UKB RAP storage.
		workflow_id: ID of the plink2 GWAS workflow on the UKB RAP.
			Defaults to WORKFLOW_ID constant.
		instance_type: Instance type to use for the workflow. Defaults
			to DEFAULT_INSTANCE constant.
		name: Name of the job. Defaults to 'gwas_plink2'.
	"""

	# Get workflow
	workflow = dxpy.dxworkflow.DXWorkflow(dxid=workflow_id)

	# Get data links for inputs
	geno_pgen_link = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=geno_file.split('/')[-1] + '.pgen',
			folder='/'.join(geno_file.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	geno_psam_link = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=geno_file.split('/')[-1] + '.psam',
			folder='/'.join(geno_file.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	geno_pvar_link = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=geno_file.split('/')[-1] + '.pvar',
			folder='/'.join(geno_file.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	covar_link = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=covar_file.split('/')[-1],
			folder='/'.join(covar_file.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	pheno_link = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=pheno_file.split('/')[-1],
			folder='/'.join(pheno_file.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	split_link = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=split_file.split('/')[-1],
			folder='/'.join(split_file.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)

	# Set up workflow input
	prefix = 'stage-common.'
	workflow_input = {
		f'{prefix}geno_pgen_file': geno_pgen_link,
		f'{prefix}geno_psam_file': geno_psam_link,
		f'{prefix}geno_pvar_file': geno_pvar_link,
		f'{prefix}covar_file': covar_link,
		f'{prefix}pheno_file': pheno_link,
		f'{prefix}split_file': split_link
	}

	# Run workflow
	analysis = workflow.run(
		workflow_input,
		folder=output_dir,
		name=name,
		instance_type=instance_type,
	)
	print("Started analysis %s (%s)\n"%(analysis.get_id(), name))

	return analysis


if __name__ == '__main__':
	# Parse args
	args = parse_args()
	print(f'Phenotype: {args.pheno_name}')
	
	# Open phenotype metadata file
	with open(args.pheno_metadata_file, 'r') as f:
		pheno_metadata = json.load(f)

	# Get the covariate set to use
	if 'covar_set' in pheno_metadata[args.pheno_name]:
		covar_set = pheno_metadata[args.pheno_name]['covar_set']
	else:
		covar_set = args.default_covar_set

	print(f'Using covariate set {covar_set}.')

	# Set genotype file and split file names based on if only using WB
	if args.wb:
		geno_fname = 'allchr_wbqc'
		split_fname = 'train_wb.txt'
	else:
		geno_fname = 'allchr_allqc'
		split_fname = 'train_all.txt'

	# Add dev suffix to split file name if using dev set
	if args.dev:
		split_fname = split_fname.replace('.txt', '_dev.txt')

	print(f'Genotype file: {geno_fname}')
	print(f'Split file: {split_fname}')

	# Set the output directory
	output_dir = f'{args.output_dir}/{args.pheno_name}_glm'
	if args.wb:
		output_dir += '_wb'
	if args.dev:
		output_dir += '_dev'
	print(f'Output directory: {output_dir}')

	# Launch the GWAS workflow
	job_name = f'gwas_plink2_{args.pheno_name}'
	if args.wb:
		job_name += '_wb'
	if args.dev:
		job_name += '_dev'

	print(f'Launching GWAS workflow with name: {job_name}')
	launch_gwas_workflow(
		geno_file=f'{args.geno_dir}/{geno_fname}',
		covar_file=f'{args.covar_dir}/{covar_set}.tsv',
		pheno_file=f'{args.pheno_dir}/{args.pheno_name}.pheno',
		split_file=f'{args.splits_dir}/{split_fname}',
		output_dir=output_dir,
		name=job_name
	)