"""Launch variant filtering before AutoML PRS can be fit.

Uses PGEN input genotype file.

Required args:

* -p, --pheno-name: Phenotype name. Must be the prefix of a directory in
	--basil-out-dir.

Optional args:

* --basil-out-dir: Directory containing directories with the BASIL output
	files. Default: '/rdevito/nonlin_prs/batch_iterative_prs/output'.
* --basil-desc: String signifying BASIL regularization option. Default
	is 'lasso'. The directory containing the BASIL output files should be
	'{pheno_name}[_wb]_{basil_desc}'.
* --wb: Flag to use clumps and QC from just the white British
	subset of the UKB data. False when not provided.
* -g, --geno-dir: Directory containing the genotype files. Default: 
	'/rdevito/nonlin_prs/data/geno_data/qced_common/pgen'
* -o, --out-dir: Directory in which to save the filtered variant file
	and the JSON file with information on what variants are includes
	under which p-value and window thresholds. Default:
	'/rdevito/nonlin_prs/automl_prs/prepro_data'.
"""

import argparse

import dxpy


WORKFLOW_ID = 'workflow-GjV2jkjJv7BPKQgvkZJVFjZ7'
DEFAULT_INSTANCE = 'mem3_ssd1_v2_x64'


def parse_args():
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter
	)
	parser.add_argument(
		'-p', '--pheno-name',
		required=True,
		help='Phenotype name. Must be the name of a directory in '
			'--prsice-out-dir excluding an optional \'_wb\' suffix.'
	)
	parser.add_argument(
		'--basil-out-dir',
		default='/rdevito/nonlin_prs/batch_iterative_prs/output',
		help='Directory containing directories with the BASIL output '
			'files. Default: \'/rdevito/nonlin_prs/batch_iterative_prs/output\'.'
	)
	parser.add_argument(
		'--basil-desc',
		default='lasso',
		help='String signifying BASIL regularization option. Default '
			'is \'lasso\'. The directory containing the BASIL output files should be '
			'\'{pheno_name}[_wb]_{basil_desc}\'.'
	)
	parser.add_argument(
		'--wb',
		action='store_true',
		help='Flag to use clumps and QC from just the white British '
			'subset of the UKB data. False when not provided.'
	)
	parser.add_argument(
		'-g', '--geno-dir',
		default='/rdevito/nonlin_prs/data/geno_data/qced_common/pgen',
		help='Directory containing the genotype files. Default: '
			'\'/rdevito/nonlin_prs/data/geno_data/qced_common/pgen\''
	)
	parser.add_argument(
		'-o', '--out-dir',
		default='/rdevito/nonlin_prs/automl_prs/prepro_data',
		help='Directory in which to save the filtered variant file '
			'and the JSON file with information on what variants are includes '
			'under which p-value and window thresholds. Default: '
			'\'/rdevito/nonlin_prs/automl_prs/prepro_data\'.'
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


def launch_automl_prepro_workflow(
	basil_incl_file,
	pgen_path,
	out_dir,
	name='prs_automl_prepro_basil'
):
	"""Launch genotype preprocessing for autoML on UKB RAP.
	
	Args:
		basil_incl_file (str): Path to the CSV file of BASIL included variants.
		pval_path (str): Path to the file containing the p-value threshold.
		pgen_path (str): Path to the PGEN file without extension.
		out_dir (str): Path to the output directory.
		name (str): Name of the workflow. Default: 'prs_automl_prepro'.
	"""

	# Get workflow
	workflow = dxpy.dxworkflow.DXWorkflow(dxid=WORKFLOW_ID)

	# Get data links for inputs
	basil_link = get_dxlink_from_path(basil_incl_file)
	geno_pgen_link = get_dxlink_from_path(f'{pgen_path}.pgen')
	geno_psam_link = get_dxlink_from_path(f'{pgen_path}.psam')
	geno_pvar_link = get_dxlink_from_path(f'{pgen_path}.pvar')

	# Set workflow input
	prefix = 'stage-common.'
	workflow_input = {
		f'{prefix}basil_incl_file': basil_link,
		f'{prefix}geno_pgen_file': geno_pgen_link,
		f'{prefix}geno_psam_file': geno_psam_link,
		f'{prefix}geno_pvar_file': geno_pvar_link,
	}

	# Run workflow
	analysis = workflow.run(
		workflow_input,
		folder=out_dir,
		name=name,
		instance_type=DEFAULT_INSTANCE,
		priority='high',
		ignore_reuse=True
	)
	print("Started analysis %s (%s)\n"%(analysis.get_id(), name))

	return analysis


if __name__ == '__main__':

	args = parse_args()
	
	# Set path to clumps file and best p-value file
	basil_dir = args.pheno_name
	if args.wb:
		basil_dir += '_wb'

	basil_dir = f'{args.basil_out_dir}/{basil_dir}_{args.basil_desc}'
	basil_incl_file = f'{basil_dir}/included_features.csv'

	# Set path to genotype file
	if args.wb:
		pgen_fname = 'allchr_wbqc'
	else:
		pgen_fname = 'allchr_allqc'
	pgen_path = f'{args.geno_dir}/{pgen_fname}'

	# Set output directory
	pheno_out_dir = args.pheno_name
	if args.wb:
		pheno_out_dir += '_wb'
	out_dir = f'{args.out_dir}/{pheno_out_dir}_basil_{args.basil_desc}'

	print(f'Output directory: {out_dir}', flush=True)

	# Set job name
	job_name = f'prs_automl_prepro_{pheno_out_dir}_basil_{args.basil_desc}'

	# Launch workflow
	launch_automl_prepro_workflow(
		basil_incl_file,
		pgen_path,
		out_dir,
		name=job_name
	)
