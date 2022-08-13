# Python script for plotting kNN-SAA, ER-SAA+kNN, ER-SAA+OLS, and PP+OLS for d_x = 10, omega = 1, and varying sigma
# Also plots Fig. 7
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from statistics import mean, quantiles
from scipy import stats


num_degrees = 3
degrees = [1.0, 0.5, 2.0]
num_data_replicates = 100
num_saa_replicates = 30
num_sample_sizes = 3
sample_sizes = np.array([55,220,1100], dtype=np.int32)
num_stdevs = 3
stdevs = np.array([5, 10, 20], dtype=np.int32)
conf_level = 0.99
t_value = stats.t.ppf(conf_level, num_saa_replicates-1)
whiskers = (5.0,95.0)
font_size = 32
line_width = 4
fig_size = (10,8)
dpi_val = 1200

# First load the FI-SAA data
SAA_data = np.zeros((num_stdevs, num_degrees, num_data_replicates, num_saa_replicates))

for std in range(num_stdevs):
	for deg in range(num_degrees):
		for data_rep in range(num_data_replicates):
			for saa_rep in range(num_saa_replicates):
				SAA_file = open("case" + str(1+std) + "_saa/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/saa_" + str(saa_rep+1) + "/saa_obj.txt", "r")
				SAA_data[std,deg,data_rep,saa_rep] = float(SAA_file.readline())
				SAA_file.close()
		
SAA_mean = np.mean(SAA_data, axis=3)


# Next, load the ER-SAA+OLS data
ERSAA_ols_data = np.zeros((num_stdevs, num_degrees, num_data_replicates, num_sample_sizes, num_saa_replicates))
ERSAA_ols_ucb = np.zeros((num_stdevs, num_degrees, num_data_replicates, num_sample_sizes))
ERSAA_samp_num = np.array([[5, 6, 7], [2, 3, 4], [2, 3, 4]], dtype=np.int32)
ERSAA_case_num = np.array([2, 9, 10], dtype=np.int32)

for std in range(num_stdevs):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(ERSAA_case_num[std]) + "_ersaa/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(ERSAA_samp_num[std][samp_size]+1) + "/ersaa_obj.txt") as ERSAA_ols_file:
					ERSAA_ols_list = ERSAA_ols_file.read().splitlines() 
					ERSAA_ols_data[std,deg,data_rep,samp_size,:] = [float(i) for i in ERSAA_ols_list]
					ERSAA_ols_gap = ERSAA_ols_data[std,deg,data_rep,samp_size,:] - SAA_data[std,deg,data_rep,:]
					ERSAA_ols_gap_mean = np.mean(ERSAA_ols_gap)
					ERSAA_ols_gap_var = np.var(ERSAA_ols_gap)
					ERSAA_ols_ucb[std,deg,data_rep,samp_size] = (100.0/np.abs(SAA_mean[std,deg,data_rep]))*(ERSAA_ols_gap_mean + t_value*np.sqrt(ERSAA_ols_gap_var/num_saa_replicates))
				ERSAA_ols_file.close()


# Next, load the ER-SAA+kNN data
ERSAA_knn_data = np.zeros((num_stdevs, num_degrees, num_data_replicates, num_sample_sizes, num_saa_replicates))
ERSAA_knn_ucb = np.zeros((num_stdevs, num_degrees, num_data_replicates, num_sample_sizes))
ERSAA_samp_num = np.array([2, 3, 4], dtype=np.int32)
ERSAA_case_num = np.array([5, 11, 12], dtype=np.int32)

for std in range(num_stdevs):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(ERSAA_case_num[std]) + "_ersaa/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/kNN/samp_" + str(ERSAA_samp_num[samp_size]+1) + "/ersaa_obj.txt") as ERSAA_knn_file:
					ERSAA_knn_list = ERSAA_knn_file.read().splitlines() 
					ERSAA_knn_data[std,deg,data_rep,samp_size,:] = [float(i) for i in ERSAA_knn_list]
					ERSAA_knn_gap = ERSAA_knn_data[std,deg,data_rep,samp_size,:] - SAA_data[std,deg,data_rep,:]
					ERSAA_knn_gap_mean = np.mean(ERSAA_knn_gap)
					ERSAA_knn_gap_var = np.var(ERSAA_knn_gap)
					ERSAA_knn_ucb[std,deg,data_rep,samp_size] = (100.0/np.abs(SAA_mean[std,deg,data_rep]))*(ERSAA_knn_gap_mean + t_value*np.sqrt(ERSAA_knn_gap_var/num_saa_replicates))
				ERSAA_knn_file.close()


# Next, load the BK data with kNN regression
BK_data = np.zeros((num_stdevs, num_degrees, num_data_replicates, num_sample_sizes, num_saa_replicates))
BK_ucb = np.zeros((num_stdevs, num_degrees, num_data_replicates, num_sample_sizes))
BK_samp_num = np.array([2, 3, 4], dtype=np.int32)
BK_case_num = np.array([2, 4, 5], dtype=np.int32)

for std in range(num_stdevs):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(BK_case_num[std]) + "_bk/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/samp_" + str(BK_samp_num[samp_size]+1) + "/bk_obj.txt") as BK_file:
					BK_list = BK_file.read().splitlines() 
					BK_data[std,deg,data_rep,samp_size,:] = [float(i) for i in BK_list]
					BK_gap = BK_data[std,deg,data_rep,samp_size,:] - SAA_data[std,deg,data_rep,:]
					BK_gap_mean = np.mean(BK_gap)
					BK_gap_var = np.var(BK_gap)
					BK_ucb[std,deg,data_rep,samp_size] = (100.0/np.abs(SAA_mean[std,deg,data_rep]))*(BK_gap_mean + t_value*np.sqrt(BK_gap_var/num_saa_replicates))
				BK_file.close()


# Next, load the Point Prediction data with OLS regression
PP_ols_data = np.zeros((num_stdevs, num_degrees, num_data_replicates, num_sample_sizes, num_saa_replicates))
PP_ols_ucb = np.zeros((num_stdevs, num_degrees, num_data_replicates, num_sample_sizes))
PP_samp_num = np.array([2, 3, 4], dtype=np.int32)
PP_case_num = np.array([2, 6, 7], dtype=np.int32)

for std in range(num_stdevs):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(PP_case_num[std]) + "_pointpred/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(PP_samp_num[samp_size]+1) + "/pointpred_obj.txt") as PP_ols_file:
					PP_ols_list = PP_ols_file.read().splitlines() 
					PP_ols_data[std,deg,data_rep,samp_size,:] = [float(i) for i in PP_ols_list]
					PP_ols_gap = PP_ols_data[std,deg,data_rep,samp_size,:] - SAA_data[std,deg,data_rep,:]
					PP_ols_gap_mean = np.mean(PP_ols_gap)
					PP_ols_gap_var = np.var(PP_ols_gap)
					PP_ols_ucb[std,deg,data_rep,samp_size] = (100.0/np.abs(SAA_mean[std,deg,data_rep]))*(PP_ols_gap_mean + t_value*np.sqrt(PP_ols_gap_var/num_saa_replicates))
				PP_ols_file.close()


# Plot PP+OLS vs ER-SAA+OLS vs ER-SAA+kNN vs kNN-SAA
box_colors = ['black', 'limegreen', 'blue', 'red']
plot_y_max = [12.0, 10.0, 30.0]
for deg in range(num_degrees):
	for std in range(num_stdevs):
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
					data[:,samp_size] = BK_ucb[std,deg,:,samp_size]
					box_positions[samp_size] = 1 + 2.5*samp_size
				elif meth == 1:
					data[:,samp_size] = ERSAA_knn_ucb[std,deg,:,samp_size]
					box_positions[samp_size] = 1.5 + 2.5*samp_size
				elif meth == 2:
					data[:,samp_size] = ERSAA_ols_ucb[std,deg,:,samp_size]
					box_positions[samp_size] = 2 + 2.5*samp_size
				elif meth == 3:
					data[:,samp_size] = PP_ols_ucb[std,deg,:,samp_size]
					box_positions[samp_size] = 2.5 + 2.5*samp_size

			plt.boxplot(data,whis=whiskers, showfliers=False, 
						positions=box_positions, widths=box_width,
						boxprops={'linewidth':line_width,'color':box_colors[meth]}, 
						medianprops={'linewidth':line_width,'color':box_colors[meth]},
						whiskerprops={'linewidth':line_width,'color':box_colors[meth]},
						capprops={'linewidth':line_width,'color':box_colors[meth]})
			if std == 0:
				plt.ylabel('UCB on % optimality gap')
			if deg == 0 or deg == 1:
				plt.title(r'$\sigma$' + ' = ' + str(stdevs[std]))
			x1,x2,y1,y2 = plt.axis()  
			plt.axis((x1,2.5*(num_sample_sizes)+0.5,0,plot_y_max[deg]))


			box_positions_overall = np.zeros(5*num_sample_sizes+1)
			box_labels_overall = np.empty(5*num_sample_sizes+1,dtype=object)
			for samp_size in range(num_sample_sizes):
				box_positions_overall[5*samp_size] = 1 + 2.5*samp_size
				box_labels_overall[5*samp_size] = 'R'
				box_positions_overall[5*samp_size+1] = 1.5 + 2.5*samp_size
				box_labels_overall[5*samp_size+1] = 'k'
				box_positions_overall[5*samp_size+2] = 1.75 + 2.5*samp_size
				box_labels_overall[5*samp_size+2] = '\n' + str(sample_sizes[samp_size])
				box_positions_overall[5*samp_size+3] = 2 + 2.5*samp_size
				box_labels_overall[5*samp_size+3] = 'O'
				box_positions_overall[5*samp_size+4] = 2.5 + 2.5*samp_size
				box_labels_overall[5*samp_size+4] = 'P'
			box_positions_overall[5*num_sample_sizes] = 0.5
			box_labels_overall[5*num_sample_sizes] = '\nn ='
			
			plt.xticks(box_positions_overall, box_labels_overall)
			plt.tick_params(axis='x',
							which='both',
							bottom=False)

			if std == 0:
				plt.subplots_adjust(left=0.16,right=0.99,bottom=0.12)
			else:
				plt.subplots_adjust(left=0.11,right=0.99,bottom=0.12)
			plt.savefig("/Users/rohitkannan/Desktop/OPRE Revision/test_cases_v2/figures/fig4_sigma" + str(std+1) + "_deg" + str(deg+1) + ".eps", format='eps')
