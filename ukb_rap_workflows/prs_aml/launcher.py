"""Launch AutoML-PRS fitting.

Required args:

* -p, --pheno-name: Phenotype name. Must be the file name of a 
	phenotype file in --pheno-dir without the '.pheno' extension.
* -m, --model-config: Filename without '.json' extension of a model
	configuration file already on docker image in the --model-config-dir
	directory.
* -d, --data-version-desc: Description of the data version to be used.
	Corresponds to suffix of directory names in --geno-dir. For example,
	'body_fat_percentage_23099_max50000_v0' would have the data version
	'max50000_v0'.

Optional args:

* --wb: Flag to use just the white British subset of the UKB data.
	False when not provided.
* -l, --large-instance: Flag to use a larger instance type. False when not
	provided.
* --pheno-dir: Directory containing the phenotype files. Default:
	'/rdevito/nonlin_prs/data/pheno_data/pheno'
* --pheno-metadata-file: File containing the phenotype metadata.
	Default: '../../data/pheno_metadata.json'
* --geno-dir: Directory containing subdirs that contain the genotype files.
	Default: '/rdevito/nonlin_prs/automl_prs/prepro_data'. Subdirs are
	of the form '{pheno-name}[_wb]_{data-version-desc}'.
* --splits-dir: Directory containing train/val/test splits in
	the form of list of sample IDs. Default: 
	'/rdevito/nonlin_prs/data/sample_data/splits'
* --covar-dir: Directory containing the covariate files. Default:
	'/rdevito/nonlin_prs/data/covar_data/tsv'
* --default_covar_set: Default covariate set file name (w/o '.tsv') 
	if no specific set of covariates is provided by the phenotype
	metadata file. If in metadata file, the key should be 'covar_set'.
	Default: 'covar_std_v1'
* --model-config-dir: Directory containing model configuration files on
	the docker image. Default: '/home/model_configs'.
* --output-dir: Directory in which a folder will be created to store
	the output of the GWAS workflow. Default:
		'/rdevito/nonlin_prs/automl_prs/output'.
	Final output directory will be of the form: 
		{output_dir}/{pheno-name}[_{wb}][_{data-version-desc}][_{model-config}]
"""

import argparse
import json

import dxpy


WORKFLOW_ID = 'workflow-Gj6yq88Jv7BBG3K4J3K5kv1Q'
DEFAULT_INSTANCE = 'mem3_ssd1_v2_x64'
LARGE_INSTANCE = 'mem3_ssd1_v2_x96'


def parse_args():
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter
	)
	parser.add_argument(
		'-p', '--pheno-name',
		required=True,
		help='Phenotype name. Must be the file name of a phenotype file in '
			'--pheno-dir without the \'.pheno\' extension.'
	)
	parser.add_argument(
		'-m', '--model-config',
		required=True,
		help='Filename without \'.json\' extension of a model configuration '
			'file already on docker image in the --model-config-dir directory.'
	)
	parser.add_argument(
		'-d', '--data-version-desc',
		required=True,
		help='Description of the data version to be used. Corresponds to '
			'suffix of directory names in --geno-dir. For example, '
			'\'body_fat_percentage_23099_max50000_v0\' would have the data '
			'version \'max50000_v0\'.'
	)
	parser.add_argument(
		'--wb',
		action='store_true',
		help='Flag to use just the white British subset of the UKB data.'
	)
	parser.add_argument(
		'-l', '--large-instance',
		action='store_true',
		help='Flag to use a larger instance type. False when not provided.'
	)
	parser.add_argument(
		'--pheno-dir',
		default='/rdevito/nonlin_prs/data/pheno_data/pheno',
		help='Directory containing the phenotype files. Default: '
			'\'/rdevito/nonlin_prs/data/pheno_data/pheno\''
	)
	parser.add_argument(
		'--pheno-metadata-file',
		default='../../data/pheno_metadata.json',
		help='File containing the phenotype metadata. Default: '
			'\'../../data/pheno_metadata.json\''
	)
	parser.add_argument(
		'--geno-dir',
		default='/rdevito/nonlin_prs/automl_prs/prepro_data',
		help='Directory containing the genotype files. Default: '
			'\'/rdevito/nonlin_prs/automl_prs/prepro_data\'.'
	)
	parser.add_argument(
		'--splits-dir',
		default='/rdevito/nonlin_prs/data/sample_data/splits',
		help='Directory containing train/val/test splits in the form of list '
			'of sample IDs. Default: \'/rdevito/nonlin_prs/data/sample_data/splits\''
	)
	parser.add_argument(
		'--covar-dir',
		default='/rdevito/nonlin_prs/data/covar_data/tsv',
		help='Directory containing the covariate files. Default: '
			'\'/rdevito/nonlin_prs/data/covar_data/tsv\'.'
	)
	parser.add_argument(
		'--default_covar_set',
		default='covar_std_v1',
		help='Default covariate set file name (w/o \'.tsv\') if no specific '
			'set of covariates is provided by the phenotype metadata file. If '
			'in metadata file, the key should be \'covar_set\'. Default: '
			'\'covar_std_v1\''
	)
	parser.add_argument(
		'--model-config-dir',
		default='/home/model_configs',
		help='Directory containing model configuration files on the docker '
			'image. Default: \'/home/model_configs\'.'
	)
	parser.add_argument(
		'--output-dir',
		default='/rdevito/nonlin_prs/automl_prs/output',
		help='Directory in which a folder will be created to store the output '
			'of the GWAS workflow. Default: \'/rdevito/nonlin_prs/automl_prs/output\'.'
	)

	return parser.parse_args()


def get_dxlink_from_path(path_to_link):
	"""Get dxlink from path."""
	print(f'Finding data object for {path_to_link}', flush=True)
	return dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=path_to_link.split('/')[-1],
			folder='/'.join(path_to_link.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)


def launch_automl_prs_workflow(
	geno_parquet,
	var_subset_json,
	pheno_file,
	covar_file,
	train_samp_file,
	val_samp_file,
	test_samp_file,
	train_config_path,
	output_dir,
	job_name='fit_prs_automl',
	instance_type=DEFAULT_INSTANCE,
):
	"""Launch AutoML-PRS fitting.

	Args:
		geno_parquet: Path to the genotype parquet file.
		var_subset_json: Path to the variant subset JSON file.
		pheno_file: Path to the phenotype file.
		covar_file: Path to the covariate file.
		train_samp_file: Path to the train sample IDs file.
		val_samp_file: Path to the validation sample IDs file.
		test_samp_file: Path to the test sample IDs file.
		train_config_path: Path to the training configuration file.
	"""

	# Get workflow
	workflow = dxpy.dxworkflow.DXWorkflow(dxid=WORKFLOW_ID)

	# Get data object IDs
	geno_parquet_dxlink = get_dxlink_from_path(geno_parquet)
	var_subset_json_dxlink = get_dxlink_from_path(var_subset_json)
	pheno_file_dxlink = get_dxlink_from_path(pheno_file)
	covar_file_dxlink = get_dxlink_from_path(covar_file)
	train_samp_file_dxlink = get_dxlink_from_path(train_samp_file)
	val_samp_file_dxlink = get_dxlink_from_path(val_samp_file)
	test_samp_file_dxlink = get_dxlink_from_path(test_samp_file)

	# Set workflow input
	prefix = 'stage-common.'
	workflow_input = {
		f'{prefix}geno_parquet': geno_parquet_dxlink,
		f'{prefix}var_subset_json': var_subset_json_dxlink,
		f'{prefix}pheno_file': pheno_file_dxlink,
		f'{prefix}covar_file': covar_file_dxlink,
		f'{prefix}train_ids': train_samp_file_dxlink,
		f'{prefix}val_ids': val_samp_file_dxlink,
		f'{prefix}test_ids': test_samp_file_dxlink,
		f'{prefix}train_config_path': train_config_path
	}

	# Run workflow
	analysis = workflow.run(
		workflow_input,
		folder=output_dir,
		name=job_name,
		instance_type=instance_type,
		priority='high',
		ignore_reuse=True
	)
	print(f'Started analysis {analysis.get_id()} ({job_name})', flush=True)

	return analysis


if __name__ == '__main__':

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

	# Set genotype parquet file
	geno_subdir = args.pheno_name
	if args.wb:
		geno_subdir += '_wb'
	geno_subdir += f'_{args.data_version_desc}'

	geno_parquet = f'{args.geno_dir}/{geno_subdir}/filtered_vars.parquet'
	var_ss_json = f'{args.geno_dir}/{geno_subdir}/filtered_vars.json'
	print(f'Genotype parquet: {geno_parquet}')
	print(f'Variant subset JSON: {var_ss_json}')

	# Set sample split ID paths
	if args.wb:
		train_samp_fname = f'{args.splits_dir}/train_wb.txt'
		val_samp_fname = f'{args.splits_dir}/val_wb.txt'
	else:
		train_samp_fname = f'{args.splits_dir}/train_all.txt'
		val_samp_fname = f'{args.splits_dir}/val_all.txt'
	test_samp_fname = f'{args.splits_dir}/test_all.txt'

	print(f'Train sample IDs file: {train_samp_fname}')
	print(f'Val sample IDs file: {val_samp_fname}')
	print(f'Test sample IDs file: {test_samp_fname}')

	# Set the output directory
	desc = args.pheno_name
	if args.wb:
		desc += '_wb'
	desc += f'_{args.data_version_desc}_{args.model_config}'

	output_dir = f'{args.output_dir}/{desc}'

	print(f'Output directory: {output_dir}')

	# Set instance type
	if args.large_instance:
		instance_type = LARGE_INSTANCE
	else:
		instance_type = DEFAULT_INSTANCE
	print(f'Instance type: {instance_type}')

	# Launch workflow
	job_name = f'prs_automl_{desc}'
	print(f'Launching AutoML-PRS workflow with name: {job_name}')

	launch_automl_prs_workflow(
		geno_parquet=geno_parquet,
		var_subset_json=var_ss_json,
		pheno_file=f'{args.pheno_dir}/{args.pheno_name}.pheno',
		covar_file=f'{args.covar_dir}/{covar_set}.tsv',
		train_samp_file=train_samp_fname,
		val_samp_file=val_samp_fname,
		test_samp_file=test_samp_fname,
		train_config_path=f'{args.model_config_dir}/{args.model_config}.json',
		output_dir=output_dir,
		job_name=job_name,
		instance_type=instance_type
	)

	print()