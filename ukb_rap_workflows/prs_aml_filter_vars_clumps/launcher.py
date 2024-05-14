"""Launch variant filtering before AutoML PRS can be fit.

Uses PGEN input genotype file.

Required args:

* -p, --pheno-name: Phenotype name. Must be the name of a directory in
	--prsice-out-dir excluding an optional '_wb' suffix.

Optional args:

* --prsice-out-dir: Directory containing the PRSice output files. Default:
	'/rdevito/nonlin_prs/sum_stats_prs/PRSice2/prsice2_output'.
* --wb: Flag to use clumps and QC from just the white British
	subset of the UKB data. False when not provided.
* -g, --geno-dir: Directory containing the genotype files. Default: 
	'/rdevito/nonlin_prs/data/geno_data/qced_common/pgen'
* -o, --out-dir: Directory in which to save the filtered variant file
	and the JSON file with information on what variants are includes
	under which p-value and window thresholds. Default:
	'/rdevito/nonlin_prs/automl_prs/prepro_data'.
* -d, --out-desc: String to be added to end of job name and output directory.
	Default: ''.
"""

import argparse

import dxpy


WORKFLOW_ID = 'workflow-GjQxJXQJv7B824Q9jQp9zyjP'
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
		'--prsice-out-dir',
		default='/rdevito/nonlin_prs/sum_stats_prs/PRSice2/prsice2_output',
		help='Directory containing the PRSice output files. Default: '
			'\'/rdevito/nonlin_prs/sum_stats_prs/PRSice2/prsice2_output\'.'
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
	parser.add_argument(
		'-d', '--out-desc',
		default='',
		help='String to be added to end of job name and output directory. Default: \'\'.'
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
	clumps_path,
	pval_path,
	pgen_path,
	out_dir,
	name='prs_automl_prepro'
):
	"""Launch genotype preprocessing for autoML on UKB RAP.
	
	Args:
		clumps_path (str): Path to the '.clumps' file.
		pval_path (str): Path to the file containing the p-value threshold.
		pgen_path (str): Path to the PGEN file without extension.
		out_dir (str): Path to the output directory.
		name (str): Name of the workflow. Default: 'prs_automl_prepro'.
	"""

	# Get workflow
	workflow = dxpy.dxworkflow.DXWorkflow(dxid=WORKFLOW_ID)

	# Get data links for inputs
	clumps_link = get_dxlink_from_path(clumps_path)
	pval_link = get_dxlink_from_path(pval_path)
	geno_pgen_link = get_dxlink_from_path(f'{pgen_path}.pgen')
	geno_psam_link = get_dxlink_from_path(f'{pgen_path}.psam')
	geno_pvar_link = get_dxlink_from_path(f'{pgen_path}.pvar')

	# Set workflow input
	prefix = 'stage-common.'
	workflow_input = {
		f'{prefix}clumps_file': clumps_link,
		f'{prefix}best_pval_file': pval_link,
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
	clumps_dir = args.pheno_name
	if args.wb:
		clumps_dir += '_wb'

	clumps_path = f'{args.prsice_out_dir}/{clumps_dir}/prs_prsice2_clump.clumps'
	pval_path = f'{args.prsice_out_dir}/{clumps_dir}/prsice2_best_p_thresh.txt'

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
	out_dir = f'{args.out_dir}/{pheno_out_dir}_clumps'

	# Set job name
	job_name = f'prs_automl_prepro_{pheno_out_dir}_clumps'

	# Launch workflow
	launch_automl_prepro_workflow(
		clumps_path,
		pval_path,
		pgen_path,
		out_dir,
		name=job_name
	)
