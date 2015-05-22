#we use this script to process datas!
#import packages
import numpy as np
data_all=np.loadtxt('results.txt',delimiter='\t',skiprows=1)
#create output file handle
fp=open('data.txt','w')
fp.write('density\t t_tilde\t energy\t energy_Variance\t Esqure\t Esquare_variance\n')
#input information from "parameters.txt"
fd=open('parameters.txt','r')
paras=fd.read()
paras_str=paras.split()
#number of bins
nbins=eval(paras_str[11])
temp=data_all.shape
N_sectors=temp[0]/nbins
data_sectors=np.array_split(data_all,N_sectors)
#now we analyze each sector
for i in range(N_sectors):
	current_sector=data_sectors[i]
	#the average and standard deviation of energy data
	E0_aver=np.mean(current_sector[:,3])
	E0_stand=np.std(current_sector[:,3])
	#the average and standard deviation of variance of energy data
	var_e0_aver=np.mean(current_sector[:,5])
	var_e0_stand=np.std(current_sector[:,5])
	t_tilde=current_sector[0,1]
	density=current_sector[0,0]
	fp.write('%0.3f\t %0.3f\t %0.6f\t %0.6f\t %0.6f\t %0.6f\n' % (density,t_tilde,E0_aver,E0_stand,var_e0_aver,var_e0_stand))
fp.close()
