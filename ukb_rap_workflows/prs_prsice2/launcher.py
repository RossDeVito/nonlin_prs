"""Launch PRSice2 C+T PRS workflow.

Uses BED input genotype file.

Required args:

* -p, --pheno-name: Phenotype name. Must be the file name of a 
	phenotype file in --pheno-dir without the '.pheno' extension.

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
	'/rdevito/nonlin_prs/data/geno_data/qced_common/bed'
* --splits-dir: Directory containing train/val/test splits in
	the form of list of sample IDs. Default: 
	'/rdevito/nonlin_prs/data/sample_data/splits'
* --covar-dir: Directory containing the covariate files. Default:
	'/rdevito/nonlin_prs/data/covar_data/tsv'
* --default_covar_set: Default covariate set file name (w/o '.tsv') 
	if no specific set of covariates is provided by the phenotype
	metadata file. If in metadata file, the key should be 'covar_set'.
	Default: 'covar_std_v1'
* --sum-stats-dir: Directory containing directories containing summary
	statistics files from gwas_plink2.py. Default:
	'/rdevito/nonlin_prs/gwas/gwas_output'. Directory where the
	appropriate summary statistics file will be found is of the form:
	{sum_stats_dir}/{pheno_name}_glm[_wb][_dev]/. The filename for the
	summary statistics file is assumed to be 
	'gwas_plink2.{pheno_name}.glm.linear'.
* --output-dir: Directory in which a folder will be created to store
	the output of the GWAS workflow. Default: 
	'/rdevito/nonlin_prs/sum_stats_prs/PRSice2/prsice2_output'. Final
	output directory will be of the form: {output_dir}/{pheno_name}[_wb][_dev]
"""

import argparse
import json

import dxpy


WORKFLOW_ID = 'workflow-GjQv0BjJv7B90q8GpqyKYzx9'
DEFAULT_INSTANCE = 'mem3_ssd1_v2_x64'


def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'-p', '--pheno-name',
		required=True,
		help='Phenotype name. Must be the file name of a phenotype file in '
		'--pheno-dir without the \'.pheno\' extension.'
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
		default='/rdevito/nonlin_prs/data/geno_data/qced_common/bed',
		help='Directory containing the genotype files.'
	)
	parser.add_argument(
		'--splits-dir',
		default='/rdevito/nonlin_prs/data/sample_data/splits',
		help='Directory containing train/val/test splits in the form of list '
		'of sample IDs.'
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
		'set of covariates is provided by the phenotype metadata file.'
	)
	parser.add_argument(
		'--sum-stats-dir',
		default='/rdevito/nonlin_prs/gwas/gwas_output',
		help='Directory containing directories containing summary statistics '
		'files from gwas_plink2.py. Directory where the appropriate summary '
		'statistics file will be found is of the form: '
		'{sum_stats_dir}/{pheno_name}_glm[_wb][_dev]/. The filename for the '
		'summary statistics file is assumed to be '
		'\'gwas_plink2.{pheno_name}.glm.linear\'.'
	)
	parser.add_argument(
		'--output-dir',
		default='/rdevito/nonlin_prs/sum_stats_prs/PRSice2/prsice2_output',
		help='Directory in which a folder will be created to store the output '
		'of the GWAS workflow.'
	)
	return parser.parse_args()


def launch_prsice2_workflow(
	geno_prefix_val,
	geno_prefix_all,
	sum_stats_file,
	pheno_file,
	covar_file,
	keep_file,
	pred_file,
	output_dir,
	workflow_id=WORKFLOW_ID,
	instance_type=DEFAULT_INSTANCE,
	name='prs_prsice2'
):
	"""Launch PRSice2 C+T PRS workflow on UKB RAP.
	
	Args:
		geno_prefix_val: Path to the genotype file prefix in UKB RAP storage
			for the validation set the PRSice2 will fit on. The filename
			should exclude the .bed/.bim/.fam extensions.
		geno_prefix_all: Path to the genotype file prefix in UKB RAP storage
			for the full dataset. Used to compute PRS scores for all samples
			before wrapper is applied. The filename should exclude the
			.bed/.bim/.fam extensions.
		sum_stats_file: Path to the summary statistics file in UKB RAP storage.
		pheno_file: Path to the phenotype file in UKB RAP storage.
		covar_file: Path to the covariate file in UKB RAP storage.
		keep_file: Path to the keep file in UKB RAP storage. Samples in this
			file will be used to fit the PRSice2 model.
		pred_file: Path to the pred file in UKB RAP storage. Samples in this
			file will have scores predicted for them, but will not be fit on.
		output_dir: Path to the output directory in UKB RAP storage.
		workflow_id: ID of the PRSice2 C+T PRS workflow.
		instance_type: Instance type to use for the workflow.
		name: Name of the workflow.
	"""

	# Get workflow
	workflow = dxpy.dxworkflow.DXWorkflow(dxid=workflow_id)

	# Get data links for inputs
	geno_bed_link_val = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=geno_prefix_val.split('/')[-1] + '.bed',
			folder='/'.join(geno_prefix_val.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	geno_bim_link_val = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=geno_prefix_val.split('/')[-1] + '.bim',
			folder='/'.join(geno_prefix_val.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	geno_fam_link_val = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=geno_prefix_val.split('/')[-1] + '.fam',
			folder='/'.join(geno_prefix_val.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	geno_bed_link_all = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=geno_prefix_all.split('/')[-1] + '.bed',
			folder='/'.join(geno_prefix_all.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	geno_bim_link_all = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=geno_prefix_all.split('/')[-1] + '.bim',
			folder='/'.join(geno_prefix_all.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	geno_fam_link_all = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=geno_prefix_all.split('/')[-1] + '.fam',
			folder='/'.join(geno_prefix_all.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	sum_stats_link = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=sum_stats_file.split('/')[-1],
			folder='/'.join(sum_stats_file.split('/')[:-1]),
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
	covar_link = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=covar_file.split('/')[-1],
			folder='/'.join(covar_file.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	keep_link = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=keep_file.split('/')[-1],
			folder='/'.join(keep_file.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)
	pred_link = dxpy.dxlink(
		list(dxpy.find_data_objects(
			name=pred_file.split('/')[-1],
			folder='/'.join(pred_file.split('/')[:-1]),
			project=dxpy.PROJECT_CONTEXT_ID
		))[0]['id']
	)

	# Set up workflow input
	prefix = 'stage-common.'
	workflow_input = {
		f'{prefix}geno_bed_file_val': geno_bed_link_val,
		f'{prefix}geno_bim_file_val': geno_bim_link_val,
		f'{prefix}geno_fam_file_val': geno_fam_link_val,
		f'{prefix}geno_bed_file_all': geno_bed_link_all,
		f'{prefix}geno_bim_file_all': geno_bim_link_all,
		f'{prefix}geno_fam_file_all': geno_fam_link_all,
		f'{prefix}sum_stats_file': sum_stats_link,
		f'{prefix}pheno_file': pheno_link,
		f'{prefix}covar_file': covar_link,
		f'{prefix}keep_file': keep_link,
		f'{prefix}pred_file': pred_link,
	}

	# Run workflow
	analysis = workflow.run(
		workflow_input,
		folder=output_dir,
		name=name,
		instance_type=instance_type,
		ignore_reuse=True
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

	# Set genotype data paths
	if args.wb:
		if args.dev:
			bed_fname_val = 'allchr_wbqc_dev'
			bed_fname_all = 'allchr_wbqc_dev'
		else:
			bed_fname_val = 'allchr_wbqc_val'
			bed_fname_all = 'allchr_wbqc'
	else:
		if args.dev:
			bed_fname_val = 'allchr_allqc_dev'
			bed_fname_all = 'allchr_allqc_dev'
		else:
			bed_fname_val = 'allchr_allqc_val'
			bed_fname_all = 'allchr_allqc'

	geno_prefix_val = f'{args.geno_dir}/{bed_fname_val}'
	geno_prefix_all = f'{args.geno_dir}/{bed_fname_all}'
	print(f'Genotype file prefix (val): {bed_fname_val}')

	# Set path to summary statistics file from GWAS
	if args.wb:
		if args.dev:
			sum_stats_dir = f'{args.sum_stats_dir}/{args.pheno_name}_glm_wb_dev'
		else:
			sum_stats_dir = f'{args.sum_stats_dir}/{args.pheno_name}_glm_wb'
	else:
		if args.dev:
			sum_stats_dir = f'{args.sum_stats_dir}/{args.pheno_name}_glm_dev'
		else:
			sum_stats_dir = f'{args.sum_stats_dir}/{args.pheno_name}_glm'
	
	sum_stats_file = f'{sum_stats_dir}/gwas_plink2.{args.pheno_name}.glm.linear'
	print(f'Summary statistics file: {sum_stats_file}')

	# Set path to keep file containing sample IDs to keep
	# and pred file containing sample IDs to predict for but not fit on
	if args.dev:
		if args.wb:
			keep_fname = 'val_wb_dev.txt'
			pred_fname = 'test_wb_dev.txt'
		else:
			keep_fname = 'val_all_dev.txt'
			pred_fname = 'test_all_dev.txt'
	else:
		if args.wb:
			keep_fname = 'val_wb.txt'
		else:
			keep_fname = 'val_all.txt'
		pred_fname = 'test_all.txt'

	keep_file = f'{args.splits_dir}/{keep_fname}'
	pred_file = f'{args.splits_dir}/{pred_fname}'
	print(f'Keep file: {keep_file}')
	print(f'Pred file: {pred_file}')

	# Set the output directory
	output_dir = f'{args.output_dir}/{args.pheno_name}'
	if args.wb:
		output_dir += '_wb'
	if args.dev:
		output_dir += '_dev'
	print(f'Output directory: {output_dir}')

	# Launch the PRSice2 workflow
	job_name = f'prs_prsice2_{args.pheno_name}'
	if args.wb:
		job_name += '_wb'
	if args.dev:
		job_name += '_dev'

	print(f'Launching PRSice2 workflow with name: {job_name}')
	launch_prsice2_workflow(
		geno_prefix_val=geno_prefix_val,
		geno_prefix_all=geno_prefix_all,
		sum_stats_file=sum_stats_file,
		pheno_file=f'{args.pheno_dir}/{args.pheno_name}.pheno',
		covar_file=f'{args.covar_dir}/{covar_set}.tsv',
		keep_file=keep_file,
		pred_file=pred_file,
		output_dir=output_dir,
		name=job_name,
	)

