"""Make and save table of PRS performance metrics for all models.


"""

import json

import dxpy
import pandas as pd


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


def get_scores_dict(scores_json_path):
	"""Return scores.json in storage as a dictionary."""
	# Download scores JSON to temp dir
	scores_link = get_dxlink_from_path(scores_json_path)
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

	save_dir = 'scores'
	temp_dir = 'tmp'
	scores_fname_local = f'{temp_dir}/scores.json'

	# Get scores
	val_scores = []
	test_all_scores = []
	test_wb_scores = []
	test_nwb_scores = []

	# Get scores for prsice-2
	model_type_dir = '/rdevito/nonlin_prs/sum_stats_prs/PRSice2/prsice2_output'
	
	for pheno in phenos:

		# Get scores
		all_scores = get_scores_dict(f'{model_type_dir}/{pheno}/scores.json')
		wb_scores = get_scores_dict(f'{model_type_dir}/{pheno}_wb/scores.json')

		# Append to appropriate lists after adding 'model_class',
		# 'model_desc', 'pheno', and 'white_british' keys
		val_scores.append({
			'model_class': 'C+T',
			'model_desc': 'PRSice-2',
			'pheno': pheno,
			'white_british': False,
			**all_scores['val']
		})
		val_scores.append({
			'model_class': 'C+T',
			'model_desc': 'PRSice-2',
			'pheno': pheno,
			'white_british': True,
			**wb_scores['val']
		})

		test_all_scores.append({
			'model_class': 'C+T',
			'model_desc': 'PRSice-2',
			'pheno': pheno,
			'white_british': False,
			**all_scores['test']
		})
		test_all_scores.append({
			'model_class': 'C+T',
			'model_desc': 'PRSice-2',
			'pheno': pheno,
			'white_british': True,
			**wb_scores['test']
		})

		test_wb_scores.append({
			'model_class': 'C+T',
			'model_desc': 'PRSice-2',
			'pheno': pheno,
			'white_british': False,
			**all_scores['test_wb']
		})
		test_wb_scores.append({
			'model_class': 'C+T',
			'model_desc': 'PRSice-2',
			'pheno': pheno,
			'white_british': True,
			**wb_scores['test_wb']
		})

		test_nwb_scores.append({
			'model_class': 'C+T',
			'model_desc': 'PRSice-2',
			'pheno': pheno,
			'white_british': False,
			**all_scores['test_nwb']
		})
		test_nwb_scores.append({
			'model_class': 'C+T',
			'model_desc': 'PRSice-2',
			'pheno': pheno,
			'white_british': True,
			**wb_scores['test_nwb']
		})

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
