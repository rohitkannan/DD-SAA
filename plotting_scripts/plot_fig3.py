# Python script for plotting J-SAA+OLS, ER-SAA+OLS, ER-SAA+Lasso, and PP+Lasso for omega = 1
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from statistics import mean, quantiles
from scipy import stats


num_degrees = 3
degrees = [1.0, 0.5, 2.0]
num_covariate_dimensions = 2
covariate_dimensions = np.array([10, 100], dtype=np.int32)
num_data_replicates = 100
num_saa_replicates = 30
num_sample_sizes = 3
sample_sizes = np.array([[14, 16, 22], [131, 152, 202]], dtype=np.int32)
conf_level = 0.99
t_value = stats.t.ppf(conf_level, num_saa_replicates-1)
whiskers = (5.0,95.0)
font_size = 32
line_width = 4
fig_size = (10,8)
dpi_val = 1200

# First load the FI-SAA data
SAA_data = np.zeros((num_degrees, num_data_replicates, num_saa_replicates))

for deg in range(num_degrees):
	for data_rep in range(num_data_replicates):
		for saa_rep in range(num_saa_replicates):
			SAA_file = open("case1_saa/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/saa_" + str(saa_rep+1) + "/saa_obj.txt", "r")
			SAA_data[deg,data_rep,saa_rep] = float(SAA_file.readline())
			SAA_file.close()
		
SAA_mean = np.mean(SAA_data, axis=2)


# Next, load the ER-SAA+OLS data
ERSAA_ols_data = np.zeros((num_covariate_dimensions, num_degrees, num_data_replicates, num_sample_sizes, num_saa_replicates))
ERSAA_ols_ucb = np.zeros((num_covariate_dimensions, num_degrees, num_data_replicates, num_sample_sizes))
ERSAA_samp_num = np.array([[1, 2, 3], [1, 2, 3]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(2+cov_dim) + "_ersaa/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(ERSAA_samp_num[cov_dim][samp_size]+1) + "/ersaa_obj.txt") as ERSAA_ols_file:
					ERSAA_ols_list = ERSAA_ols_file.read().splitlines() 
					ERSAA_ols_data[cov_dim,deg,data_rep,samp_size,:] = [float(i) for i in ERSAA_ols_list]
					ERSAA_ols_gap = ERSAA_ols_data[cov_dim,deg,data_rep,samp_size,:] - SAA_data[deg,data_rep,:]
					ERSAA_ols_gap_mean = np.mean(ERSAA_ols_gap)
					ERSAA_ols_gap_var = np.var(ERSAA_ols_gap)
					ERSAA_ols_ucb[cov_dim,deg,data_rep,samp_size] = (100.0/np.abs(SAA_mean[deg,data_rep]))*(ERSAA_ols_gap_mean + t_value*np.sqrt(ERSAA_ols_gap_var/num_saa_replicates))
				ERSAA_ols_file.close()


# Next, load the ER-SAA+Lasso data
ERSAA_lasso_data = np.zeros((num_covariate_dimensions, num_degrees, num_data_replicates, num_sample_sizes, num_saa_replicates))
ERSAA_lasso_ucb = np.zeros((num_covariate_dimensions, num_degrees, num_data_replicates, num_sample_sizes))
ERSAA_samp_num = np.array([[1, 2, 3], [1, 2, 3]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(7+cov_dim) + "_ersaa/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/lasso/samp_" + str(ERSAA_samp_num[cov_dim][samp_size]+1) + "/ersaa_obj.txt") as ERSAA_lasso_file:
					ERSAA_lasso_list = ERSAA_lasso_file.read().splitlines() 
					ERSAA_lasso_data[cov_dim,deg,data_rep,samp_size,:] = [float(i) for i in ERSAA_lasso_list]
					ERSAA_lasso_gap = ERSAA_lasso_data[cov_dim,deg,data_rep,samp_size,:] - SAA_data[deg,data_rep,:]
					ERSAA_lasso_gap_mean = np.mean(ERSAA_lasso_gap)
					ERSAA_lasso_gap_var = np.var(ERSAA_lasso_gap)
					ERSAA_lasso_ucb[cov_dim,deg,data_rep,samp_size] = (100.0/np.abs(SAA_mean[deg,data_rep]))*(ERSAA_lasso_gap_mean + t_value*np.sqrt(ERSAA_lasso_gap_var/num_saa_replicates))
				ERSAA_lasso_file.close()


# Next, load the J-SAA+OLS data
JSAA_ols_data = np.zeros((num_covariate_dimensions, num_degrees, num_data_replicates, num_sample_sizes, num_saa_replicates))
JSAA_ols_ucb = np.zeros((num_covariate_dimensions, num_degrees, num_data_replicates, num_sample_sizes))
JSAA_samp_num = np.array([[1, 2, 3], [1, 2, 3]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(1+cov_dim) + "_jsaa/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(JSAA_samp_num[cov_dim][samp_size]+1) + "/jsaa_obj.txt") as JSAA_ols_file:
					JSAA_ols_list = JSAA_ols_file.read().splitlines() 
					JSAA_ols_data[cov_dim,deg,data_rep,samp_size,:] = [float(i) for i in JSAA_ols_list]
					JSAA_ols_gap = JSAA_ols_data[cov_dim,deg,data_rep,samp_size,:] - SAA_data[deg,data_rep,:]
					JSAA_ols_gap_mean = np.mean(JSAA_ols_gap)
					JSAA_ols_gap_var = np.var(JSAA_ols_gap)
					JSAA_ols_ucb[cov_dim,deg,data_rep,samp_size] = (100.0/np.abs(SAA_mean[deg,data_rep]))*(JSAA_ols_gap_mean + t_value*np.sqrt(JSAA_ols_gap_var/num_saa_replicates))
				JSAA_ols_file.close()


# Next, load the Point Prediction data with Lasso regression
PP_lasso_data = np.zeros((num_covariate_dimensions, num_degrees, num_data_replicates, num_sample_sizes, num_saa_replicates))
PP_lasso_ucb = np.zeros((num_covariate_dimensions, num_degrees, num_data_replicates, num_sample_sizes))
PP_samp_num = np.array([[1, 2, 3], [1, 2, 3]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(4+cov_dim) + "_pointpred/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/lasso/samp_" + str(PP_samp_num[cov_dim][samp_size]+1) + "/pointpred_obj.txt") as PP_lasso_file:
					PP_lasso_list = PP_lasso_file.read().splitlines() 
					PP_lasso_data[cov_dim,deg,data_rep,samp_size,:] = [float(i) for i in PP_lasso_list]
					PP_lasso_gap = PP_lasso_data[cov_dim,deg,data_rep,samp_size,:] - SAA_data[deg,data_rep,:]
					PP_lasso_gap_mean = np.mean(PP_lasso_gap)
					PP_lasso_gap_var = np.var(PP_lasso_gap)
					PP_lasso_ucb[cov_dim,deg,data_rep,samp_size] = (100.0/np.abs(SAA_mean[deg,data_rep]))*(PP_lasso_gap_mean + t_value*np.sqrt(PP_lasso_gap_var/num_saa_replicates))
				PP_lasso_file.close()


# Plot PP+Lasso ER-SAA+Lasso vs ER-SAA+OLS vs J-SAA+OLS
box_colors = ['black', 'limegreen', 'blue', 'red']
plot_y_max = [10.0, 12.0, 50.0]
for deg in range(num_degrees):
	for cov in range(num_covariate_dimensions):
		plt.figure(figsize=fig_size, dpi=dpi_val)
		plt.rcParams["font.weight"] = "bold"
		plt.rcParams["axes.labelweight"] = "bold"
		plt.rcParams["axes.titleweight"] = "bold"
		plt.rcParams.update({'font.size': font_size})

		for meth in range(4):
			data = np.zeros((num_data_replicates,num_sample_sizes))
			box_positions = np.zeros(num_sample_sizes)
			box_width = 0.5
			for samp_size in range(num_sample_sizes):
				if meth == 0:
					data[:,samp_size] = JSAA_ols_ucb[cov,deg,:,samp_size]
					box_positions[samp_size] = 1 + 2.5*samp_size
				elif meth == 1:
					data[:,samp_size] = ERSAA_ols_ucb[cov,deg,:,samp_size]
					box_positions[samp_size] = 1.5 + 2.5*samp_size
				elif meth == 2:
					data[:,samp_size] = ERSAA_lasso_ucb[cov,deg,:,samp_size]
					box_positions[samp_size] = 2 + 2.5*samp_size
				elif meth == 3:
					data[:,samp_size] = PP_lasso_ucb[cov,deg,:,samp_size]
					box_positions[samp_size] = 2.5 + 2.5*samp_size

			plt.boxplot(data,whis=whiskers, showfliers=False, 
						positions=box_positions, widths=box_width,
						boxprops={'linewidth':line_width,'color':box_colors[meth]}, 
						medianprops={'linewidth':line_width,'color':box_colors[meth]},
						whiskerprops={'linewidth':line_width,'color':box_colors[meth]},
						capprops={'linewidth':line_width,'color':box_colors[meth]})
			if cov == 0:
				plt.ylabel('UCB on % optimality gap')
			if deg == 0:
				plt.title(r'$d_x$' + ' = ' + str(covariate_dimensions[cov]))
			x1,x2,y1,y2 = plt.axis()  
			plt.axis((x1,2.5*(num_sample_sizes)+0.5,0,plot_y_max[deg]))


			box_positions_overall = np.zeros(5*num_sample_sizes+1)
			box_labels_overall = np.empty(5*num_sample_sizes+1,dtype=object)
			for samp_size in range(num_sample_sizes):
				box_positions_overall[5*samp_size] = 1 + 2.5*samp_size
				box_labels_overall[5*samp_size] = 'J'
				box_positions_overall[5*samp_size+1] = 1.5 + 2.5*samp_size
				box_labels_overall[5*samp_size+1] = 'O'
				box_positions_overall[5*samp_size+2] = 1.75 + 2.5*samp_size
				box_labels_overall[5*samp_size+2] = '\n' + str(sample_sizes[cov][samp_size])
				box_positions_overall[5*samp_size+3] = 2 + 2.5*samp_size
				box_labels_overall[5*samp_size+3] = 'L'
				box_positions_overall[5*samp_size+4] = 2.5 + 2.5*samp_size
				box_labels_overall[5*samp_size+4] = 'P\''
			box_positions_overall[5*num_sample_sizes] = 0.5
			box_labels_overall[5*num_sample_sizes] = '\nn ='
			
			plt.xticks(box_positions_overall, box_labels_overall)
			plt.tick_params(axis='x',
							which='both',
							bottom=False)

			if cov == 0:
				plt.subplots_adjust(left=0.16,right=0.99,bottom=0.12)
			else:
				plt.subplots_adjust(left=0.11,right=0.99,bottom=0.12)
			plt.savefig("/Users/rohitkannan/Desktop/OPRE Revision/test_cases_v2/figures/fig3_cov" + str(cov+1) + "_deg" + str(deg+1) + ".eps", format='eps')
