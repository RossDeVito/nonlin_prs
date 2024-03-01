"""Launch PRS scoring and plotting workflow.

Required args:

* -m, --model-type: Model type. One of "prsice" or "basil"
* -p, --pheno-name: Phenotype name

Optional args:
* --wb: Flag for model fit on White British only
"""

import argparse

import dxpy


WORKFLOW_ID = 'workflow-GgQgb00Jv7BJqj1qP1XXfyzF'
DEFAULT_INSTANCE = 'mem1_ssd1_v2_x2'

PHENO_DIR = '/rdevito/nonlin_prs/data/pheno_data/pheno'
SPLIT_DIR = '/rdevito/nonlin_prs/data/sample_data/splits'
TEST_WB_SPLIT_FNAME = 'test_wb.txt'


def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'-m', '--model-type',
		choices=['prsice', 'basil'],
		required=True,
		help='Model type'
	)
	parser.add_argument(
		'-p', '--pheno-name',
		required=True,
		help='Phenotype name'
	)
	parser.add_argument(
		'--wb',
		action='store_true',
		help='Flag for model fit on White British only'
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


def launch_workflow(
	model_dir,
	pheno_file,
	wb_split_file,
	instance_type=DEFAULT_INSTANCE,
	name='score_prs_preds'
):
	"""Launch PRS scoring and plotting workflow."""
	
	# Get links
	val_pred_link = get_dxlink_from_path(f'{model_dir}/val_preds.csv')
	test_pred_link = get_dxlink_from_path(f'{model_dir}/test_preds.csv')
	pheno_link = get_dxlink_from_path(pheno_file)
	wb_split_link = get_dxlink_from_path(wb_split_file)

	# Set up workflow input
	prefix = 'stage-common.'
	workflow_input = {
		f'{prefix}val_preds': val_pred_link,
		f'{prefix}test_preds': test_pred_link,
		f'{prefix}pheno_file': pheno_link,
		f'{prefix}test_wb_samples': wb_split_link
	}

	# Get workflow
	workflow = dxpy.dxworkflow.DXWorkflow(dxid=WORKFLOW_ID)

	# Run workflow
	analysis = workflow.run(
		workflow_input,
		folder=model_dir,
		name=name,
		instance_type=instance_type,
	)
	print("Started analysis %s (%s)\n"%(analysis.get_id(), name))

	return analysis


if __name__ == '__main__':

	args = parse_args()

	# Set model output dir and model dir
	if args.model_type == 'prsice':
		model_out_dir = '/rdevito/nonlin_prs/sum_stats_prs/PRSice2/prsice2_output'
		model_dir = f'{model_out_dir}/{args.pheno_name}'
		if args.wb:
			model_dir = f'{model_dir}_wb'
	elif args.model_type == 'basil':
		raise NotImplementedError('BASIL model not yet supported')
	else:
		raise ValueError(f'Invalid model type: {args.model_type}')
	
	# Set pheno and split file paths
	pheno_file = f'{PHENO_DIR}/{args.pheno_name}.pheno'
	split_file = f'{SPLIT_DIR}/{TEST_WB_SPLIT_FNAME}'

	# Launch workflow
	if args.wb:
		name = f'score_prs_preds_{args.model_type}_wb_{args.pheno_name}'
	else:
		name = f'score_prs_preds_{args.model_type}_{args.pheno_name}'

	launch_workflow(
		model_dir,
		pheno_file,
		split_file,
		instance_type=DEFAULT_INSTANCE,
		name=name
	)
