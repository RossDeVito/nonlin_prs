"""Launch variant filtering before AutoML PRS can be fit.

Uses PGEN input genotype file.

Required args:

* -p, --pheno-name: Phenotype name. Must be the prefix of a 
	directory in --sum-stats-dir before the '_glm*' part.

Optional args:

* -m, --max-variants: Maximum number of variants to include in any one
	set of variants at a threshold. Default: 80,000.
* -s, --sum-stats-dir: Directory containing the summary statistics
	files. Default: '/rdevito/nonlin_prs/gwas/gwas_output'.
* --wb: Flag to use summary statistics and QC from just the white British
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


WORKFLOW_ID = 'workflow-GjV317jJv7B9QX1qPV1zgxXB'
DEFAULT_INSTANCE = 'mem3_ssd1_v2_x64'


def parse_args():
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter
	)
	parser.add_argument(
		'-p', '--pheno-name',
		required=True,
		help='Phenotype name. Must be the prefix of a directory in '
			'--sum-stats-dir before the \'_glm*\' part.'
	)
	parser.add_argument(
		'-m', '--max-variants',
		type=int,
		default=80_000,
		help='Maximum number of variants to include in any one set of '
			'variants at a threshold. Default: 80,000.'
	)
	parser.add_argument(
		'-s', '--sum-stats-dir',
		default='/rdevito/nonlin_prs/gwas/gwas_output',
		help='Directory containing the summary statistics files. Default: '
			'\'/rdevito/nonlin_prs/gwas/gwas_output\'.'
	)
	parser.add_argument(
		'--wb',
		action='store_true',
		help='Flag to use summary statistics from just the white British '
			'subset of the UKB data. False when not provided.'
	)
	parser.add_argument(
		'-g', '--geno-dir',
		default='/rdevito/nonlin_prs/data/geno_data/qced_common/pgen',
		help='Directory containing the genotype files. Default: '
			'\'/rdevito/nonlin_prs/data/geno_data/qced_common/pgen\'.'
	)
	parser.add_argument(
		'-o', '--out-dir',
		default='/rdevito/nonlin_prs/automl_prs/prepro_data',
		help='Directory in which to save the filtered variant file and the '
			'JSON file with information on what variants are includes under '
			'which p-value and window thresholds. Default: '
			'\'/rdevito/nonlin_prs/automl_prs/prepro_data\'.'
	)
	parser.add_argument(
		'-d', '--out-desc',
		default='',
		help='String to be added to end of job name and output directory. '
			'Default: \'\''
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
	sum_stats_path,
	pgen_path,
	out_dir,
	max_num_vars,
	name='prs_automl_prepro'
):
	"""Launch genotype preprocessing for autoML on UKB RAP.
	
	Args:
		sum_stats_path (str): Path to the summary statistics file.
		pgen_path (str): Path to the PGEN file without extension.
		out_dir (str): Path to the output directory.
		name (str): Name of the workflow. Default: 'prs_automl_prepro'.
	"""

	# Get workflow
	workflow = dxpy.dxworkflow.DXWorkflow(dxid=WORKFLOW_ID)

	# Get data links for inputs
	sum_stats_link = get_dxlink_from_path(sum_stats_path)
	geno_pgen_link = get_dxlink_from_path(f'{pgen_path}.pgen')
	geno_psam_link = get_dxlink_from_path(f'{pgen_path}.psam')
	geno_pvar_link = get_dxlink_from_path(f'{pgen_path}.pvar')

	# Set workflow input
	prefix = 'stage-common.'
	workflow_input = {
		f'{prefix}sum_stats_file': sum_stats_link,
		f'{prefix}geno_pgen_file': geno_pgen_link,
		f'{prefix}geno_psam_file': geno_psam_link,
		f'{prefix}geno_pvar_file': geno_pvar_link,
		f'{prefix}max_num_vars': max_num_vars
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
	
	# Set path to sum stats file
	pheno_sum_stats_dir = args.pheno_name + '_glm'
	if args.wb:
		pheno_sum_stats_dir += '_wb'

	sum_stats_fname = f'gwas_plink2.{args.pheno_name}.glm.linear'
	sum_stats_path = f'{args.sum_stats_dir}/{pheno_sum_stats_dir}/{sum_stats_fname}'

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
	out_dir = f'{args.out_dir}/{pheno_out_dir}_max{args.max_variants}{args.out_desc}'

	# Set job name
	job_name = f'prs_automl_prepro_{pheno_out_dir}_max{args.max_variants}{args.out_desc}'

	# Launch workflow
	launch_automl_prepro_workflow(
		sum_stats_path,
		pgen_path,
		out_dir,
		max_num_vars=args.max_variants,
		name=job_name
	)
