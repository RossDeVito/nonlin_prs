"""Make and save table of PRS performance metrics for all models.


"""

import json
from collections import defaultdict

import dxpy
import pandas as pd


def get_dxlink_from_path(path_to_link):
	"""Get dxlink from path."""
	print(f'Finding data object for {path_to_link}', flush=True)

	data_objs = list(dxpy.find_data_objects(
		name=path_to_link.split('/')[-1],
		folder='/'.join(path_to_link.split('/')[:-1]),
		project=dxpy.PROJECT_CONTEXT_ID
	))

	if len(data_objs) == 0:
		return None
	elif len(data_objs) > 1:
		raise ValueError(f'Multiple data objects found for {path_to_link}')
	else:
		return dxpy.dxlink(data_objs[0]['id'])


def get_scores_dict(scores_json_path):
	"""Return scores.json in storage as a dictionary.
	
	Returns None if scores.json does not exist.
	"""
	# Download scores JSON to temp dir
	scores_link = get_dxlink_from_path(scores_json_path)

	if scores_link is None:
		return None
	
	dxpy.download_dxfile(
		scores_link,
		scores_fname_local
	)

	# Load scores
	with open(scores_fname_local, 'r') as f:
		scores = json.load(f)

	return scores


if __name__ == '__main__':

	# Options
	phenos = [
		'standing_height_50',
		'body_fat_percentage_23099',
		'platelet_count_30080',
		'glycated_haemoglobin_30750'
	]

	model_types = [
		'prsice',
		'basil_lasso',
		'automl_max50000_v0_lgbm_v0h2_p1e-8_w0',
		'automl_max50000_v0_lgbm_v0h24_p1e-5_w5000',
		'automl_max35000_v2_lgbm_v0h24_p1e-11_w250000',
		'automl_max50000_v0_elastic_net_v0h2_p1e-8_w0',
		'automl_max35000_v2_elastic_net_v0h24_p1e-11_w250000',
		'automl_max35000_v2_lgbm_v0h24',
		'automl_max35000_v2_npart_elastic_net_v0h24',
	]

	model_dirs = {
		'prsice': '/rdevito/nonlin_prs/sum_stats_prs/PRSice2/prsice2_output',
		'basil_lasso': '/rdevito/nonlin_prs/batch_iterative_prs/output',
		'automl': '/rdevito/nonlin_prs/automl_prs/output',
	}

	model_class_meta = {
		'prsice': {
			'model_class': 'C+T',
			'model_desc': 'PRSice-2'
		},
		'basil_lasso': {
			'model_class': 'Iterative screening',
			'model_desc': 'BASIL Lasso'
		},
		'automl_max50000_v0_lgbm_v0h2_p1e-8_w0': {
			'model_class': 'AutoML non-linear',
			'model_desc': 'lgbm_v0h2_p1e-8_w0'
		},
		'automl_max50000_v0_lgbm_v0h24_p1e-5_w5000': {
			'model_class': 'AutoML non-linear',
			'model_desc': 'lgbm_v0h24_p1e-5_w5000'
		},
		'automl_max35000_v2_lgbm_v0h24_p1e-11_w250000': {
			'model_class': 'AutoML non-linear',
			'model_desc': 'lgbm_v0h24_p1e-11_w250000'
		},
		'automl_max50000_v0_elastic_net_v0h2_p1e-8_w0': {
			'model_class': 'AutoML linear',
			'model_desc': 'elastic_net_v0h2_p1e-8_w0'
		},
		'automl_max35000_v2_elastic_net_v0h24_p1e-11_w250000': {
			'model_class': 'AutoML linear',
			'model_desc': 'elastic_net_v0h24_p1e-11_w250000'
		},
		'automl_max35000_v2_lgbm_v0h24': {
			'model_class': 'AutoML non-linear',
			'model_desc': 'lgbm_v0h24'
		},
		'automl_max35000_v2_npart_elastic_net_v0h24': {
			'model_class': 'AutoML linear',
			'model_desc': 'n-part elastic net v0h24'
		},
	}

	gwas_dir = '/rdevito/nonlin_prs/gwas/gwas_output'

	save_dir = 'scores'
	temp_dir = 'tmp'
	scores_fname_local = f'{temp_dir}/scores.json'

	# Get plink2 GWAS runtimes
	gwas_runtimes = defaultdict(dict)

	for pheno in phenos:
		print(f'\nGetting GWAS runtimes for {pheno}...', flush=True)

		desc = f'{pheno}_glm'

		gwas_runtimes[pheno]['all'] = get_scores_dict(
			f'{gwas_dir}/{desc}/runtime.json'
		)['runtime_seconds']	# type: ignore

		gwas_runtimes[pheno]['wb'] = get_scores_dict(
			f'{gwas_dir}/{desc}_wb/runtime.json'
		)['runtime_seconds']	# type: ignore

	# Get scores
	val_scores = []
	test_all_scores = []
	test_wb_scores = []
	test_nwb_scores = []

	for model_type in model_types:
		print(f'\nGetting scores for {model_type} models...', flush=True)

		if model_type.startswith('automl'):
			model_type_dir = model_dirs['automl']
		else:
			model_type_dir = model_dirs[model_type]

		# Get scores for each phenotype
		for pheno in phenos:
			print(f'\t{pheno}', flush=True)

			pheno_dir_template = "{:}"
			if model_type == 'basil_lasso':
				pheno_dir_template = '{:}_lasso'

			if model_type.startswith('automl'):
				desc = '_'.join(model_type.split('_')[1:])
				pheno_dir_template = '{:}_' + desc


			# Get scores for models trained on all samples
			all_scores = get_scores_dict(
				f'{model_type_dir}/{pheno_dir_template.format(pheno)}/scores.json'
			)
			if all_scores is not None:
				val_scores.append({
					'model_class': model_class_meta[model_type]['model_class'],
					'model_desc': model_class_meta[model_type]['model_desc'],
					'pheno': pheno,
					'white_british': False,
					**all_scores['val']
				})
				test_all_scores.append({
					'model_class': model_class_meta[model_type]['model_class'],
					'model_desc': model_class_meta[model_type]['model_desc'],
					'pheno': pheno,
					'white_british': False,
					**all_scores['test']
				})
				test_wb_scores.append({
					'model_class': model_class_meta[model_type]['model_class'],
					'model_desc': model_class_meta[model_type]['model_desc'],
					'pheno': pheno,
					'white_british': False,
					**all_scores['test_wb']
				})
				test_nwb_scores.append({
					'model_class': model_class_meta[model_type]['model_class'],
					'model_desc': model_class_meta[model_type]['model_desc'],
					'pheno': pheno,
					'white_british': False,
					**all_scores['test_nwb']
				})

				# Get runtime from runtime.json
				runtime = get_scores_dict(
					f'{model_type_dir}/{pheno_dir_template.format(pheno)}/runtime.json'
				)['runtime_seconds']	# type: ignore

				if not model_type.startswith('basil'):
					# Add GWAS runtime
					runtime += gwas_runtimes[pheno]['all']

				val_scores[-1]['runtime_seconds'] = runtime
				test_all_scores[-1]['runtime_seconds'] = runtime
				test_wb_scores[-1]['runtime_seconds'] = runtime
				test_nwb_scores[-1]['runtime_seconds'] = runtime

			else:
				print(f'\t\tNo scores found for {pheno} models trained on all samples')
			
			# Get scores for models trained on white British samples
			wb_scores = get_scores_dict(
				f'{model_type_dir}/{pheno_dir_template.format(pheno + "_wb")}/scores.json'
			)
			if wb_scores is not None:
				val_scores.append({
					'model_class': model_class_meta[model_type]['model_class'],
					'model_desc': model_class_meta[model_type]['model_desc'],
					'pheno': pheno,
					'white_british': True,
					**wb_scores['val']
				})
				test_all_scores.append({
					'model_class': model_class_meta[model_type]['model_class'],
					'model_desc': model_class_meta[model_type]['model_desc'],
					'pheno': pheno,
					'white_british': True,
					**wb_scores['test']
				})
				test_wb_scores.append({
					'model_class': model_class_meta[model_type]['model_class'],
					'model_desc': model_class_meta[model_type]['model_desc'],
					'pheno': pheno,
					'white_british': True,
					**wb_scores['test_wb']
				})
				test_nwb_scores.append({
					'model_class': model_class_meta[model_type]['model_class'],
					'model_desc': model_class_meta[model_type]['model_desc'],
					'pheno': pheno,
					'white_british': True,
					**wb_scores['test_nwb']
				})

				# Get runtime from runtime.json
				runtime = get_scores_dict(
					f'{model_type_dir}/{pheno_dir_template.format(pheno + "_wb")}/runtime.json'
				)['runtime_seconds']	# type: ignore

				if not model_type.startswith('basil'):
					# Add GWAS runtime
					runtime += gwas_runtimes[pheno]['wb']

				val_scores[-1]['runtime_seconds'] = runtime
				test_all_scores[-1]['runtime_seconds'] = runtime
				test_wb_scores[-1]['runtime_seconds'] = runtime
				test_nwb_scores[-1]['runtime_seconds'] = runtime
				
			else:
				print(f'\t\tNo scores found for {pheno} models trained on white British samples')

	# Make dataframes
	val_scores_df = pd.DataFrame(val_scores)
	test_all_scores_df = pd.DataFrame(test_all_scores)
	test_wb_scores_df = pd.DataFrame(test_wb_scores)
	test_nwb_scores_df = pd.DataFrame(test_nwb_scores)

	# Save
	val_scores_df.to_csv(f'{save_dir}/val_scores.csv', index=False)
	test_all_scores_df.to_csv(f'{save_dir}/test_all_scores.csv', index=False)
	test_wb_scores_df.to_csv(f'{save_dir}/test_wb_scores.csv', index=False)
	test_nwb_scores_df.to_csv(f'{save_dir}/test_nwb_scores.csv', index=False)
