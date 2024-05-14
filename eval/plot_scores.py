import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == '__main__':

	scores_dir = 'scores'
	
	# Load scores for test all, wb, and nwb sets
	test_all_scores = pd.read_csv(f'{scores_dir}/test_all_scores.csv')
	test_wb_scores = pd.read_csv(f'{scores_dir}/test_wb_scores.csv')
	test_nwb_scores = pd.read_csv(f'{scores_dir}/test_nwb_scores.csv')

	# Add 'eval_set' column
	test_all_scores['eval_set'] = 'All'
	test_wb_scores['eval_set'] = 'White British'
	test_nwb_scores['eval_set'] = 'Not White British'

	# Concatenate
	all_scores = pd.concat([test_all_scores, test_wb_scores, test_nwb_scores])

	# Plot
	all_scores['white_british'] = all_scores['white_british'].astype(str).replace({'True': 'True', 'False': 'False'})

	sns.set_style('whitegrid')

	g = sns.catplot(
		data=all_scores,
		row='pheno',
		col='eval_set',
		x='model_desc',
		y='r2',
		hue='white_british',
		kind='bar',
		margin_titles=True,
		dodge=True,
		sharey='row',	# type: ignore
	)

	g.set_titles(row_template="{row_name}")
        
	plt.show()
	plt.close()
