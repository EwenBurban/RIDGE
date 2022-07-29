#!/shared/software/miniconda/envs/r-3.5.1/bin/Rscript
# #!/shared/home/croux/.conda/envs/R_env/bin/Rscript
# #!/usr/bin/Rscript

#################################################################################################################################
#################################################################################################################################
#####                                                                                                                       #####
#####    This file is part of Demographic Inferences with Linked Selection : DILS.                                          #####
#####                                                                                                                       #####   
#####    DILS is free software: you can redistribute it and/or modify                                                       #####
#####    it under the terms of the GNU General Public License as published by                                               #####
#####    the Free Software Foundation, either version 3 of the License, or                                                  #####
#####    (at your option) any later version.                                                                                #####
#####                                                                                                                       #####    
#####    DILS is distributed in the hope that it will be useful,                                                            #####
#####    but WITHOUT ANY WARRANTY; without even the implied warranty of                                                     #####
#####    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                      #####
#####    GNU General Public License for more details.                                                                       #####
#####                                                                                                                       #####    
#####    You should have received a copy of the GNU General Public License                                                  #####
#####    along with DILS.  If not, see <https://www.gnu.org/licenses/>.                                                     #####
#####                                                                                                                       #####    
#####    Please send bugreports with examples or suggestions to                                                             #####
#####    camille.roux@univ-lille.fr                                                                                         #####
#####                                                                                                                       #####    
#####    Or write a post on https://groups.google.com/forum/#!forum/dils---demographic-inferences-with-linked-selection     #####
#####                                                                                                                       #####
#################################################################################################################################
#################################################################################################################################

for(i in commandArgs()){
	tmp = strsplit(i, '=')
	if(tmp[[1]][1] == 'Nref'){ Nref = as.double(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'nameA'){ nameA = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nameB'){ nameB = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nMin'){ nMin = as.integer(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'sub_dir_sim'){ sub_dir_sim = tmp[[1]][2] }
	if(tmp[[1]][1] == 'nSubdir'){ nSubdir = as.integer(tmp[[1]][2]) } # number of subdirectories where simulations were ran
	if(tmp[[1]][1] == 'ncores'){ ncores = as.integer(tmp[[1]][2]) } # number of cores for the random forest
	if(tmp[[1]][1] == 'ntree'){ ntree = as.integer(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'useSFS'){ useSFS = as.integer(tmp[[1]][2]) } # 0: SFS used. 1: SFS used
	if(tmp[[1]][1] == 'bestModel'){ bestModel = tmp[[1]][2] } # name of the best model
	if(tmp[[1]][1] == 'timeStamp'){ timeStamp = tmp[[1]][2] } # name of timeStamp
	if(tmp[[1]][1] == 'nPosterior'){ nPosterior = as.integer(tmp[[1]][2]) }
	if(tmp[[1]][1] == 'binpath'){ binpath = tmp[[1]][2] } # binpath = path where to source get_parameters
	if(tmp[[1]][1] == 'path2observation'){ path2observation = tmp[[1]][2] } # path2observation = path where the observed data to fit can be found
	if(tmp[[1]][1] == 'output_dir'){ output_dir = tmp[[1]][2] } # path2observation = path where the observed data to fit can be found
	if(tmp[[1]][1] == 'mode'){ mode_var2remove = tmp[[1]][2] } # path2observation = path where the observed data to fit can be found
}

outfile = paste(output_dir, '/report_', nameA, '_', nameB, '.txt', sep='')

# colors
coul = c('#ffffcc', '#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#0c2c84')
coul = colorRampPalette(coul)

######################################################
# parameters of the best model, IM_2M_2N and SI_2N
source(paste(binpath, '/get_parameters_2pop.R', sep=''))

# SI_2N; IM_2M_2N; AM_2M_2N; SC_2M_2N
options(digits=5)
#list_models_param = c('SI_1N', 'SI_2N', 'IM_2M_2N', 'AM_2M_2N', 'SC_2M_2N')
list_models_param = c(bestModel)
#for(model_tmp in list_models_param){
#	write(paste('\n#####\n\nparameters of model using neural_network (upper lines) and random_forest (lower lines): ', model_tmp, sep=''), outfile, append=T)
	posterior = get_posterior(nameA=nameA, nameB=nameB, nSubdir=nSubdir, sub_dir_sim=sub_dir_sim, model=bestModel, sub_dir_model=bestModel, nPosterior=nPosterior, figure=T, timeStamp=timeStamp, path2observation=path2observation, useSFS=useSFS, ncores=ncores,mode_var2remove=mode_var2remove)
#	write('param\tHPD2.5%\tmedian\tHPD%97.5', outfile, append=T)
####for(i in 1:ncol(posterior[['neural_network']])){
####	param_i = colnames(posterior[['neural_network']])[i]
####	if(param_i=='N1' || param_i=='N2' || param_i=='Na'){
####		write(paste(param_i, as.numeric(quantile(posterior[['neural_network']][,i], 0.025))*Nref, as.numeric(quantile(posterior[['neural_network']][,i], 0.5))*Nref, as.numeric(quantile(posterior[['neural_network']][,i], 0.975))*Nref, sep='\t'), outfile, append=T)
####		write(paste(param_i, as.numeric(posterior[['random_forest']][[param_i]][['quantile025']])*Nref, as.numeric(posterior[['random_forest']][[param_i]][['expectation']])*Nref, as.numeric(posterior[['random_forest']][[param_i]][['quantile975']])*Nref, sep='\t'), outfile, append=T)
####	
####	}else if(param_i=='Tsplit' || param_i=='Tam' || param_i=='Tsc' || param_i=='Tdem1' || param_i=='Tdem2'){
####		write(paste(param_i, as.numeric(quantile(posterior[['neural_network']][,i], 0.025))*4*Nref, as.numeric(quantile(posterior[['neural_network']][,i], 0.5))*4*Nref, as.numeric(quantile(posterior[['neural_network']][,i], 0.975))*4*Nref, sep='\t'), outfile, append=T)
####		write(paste(param_i, as.numeric(posterior[['random_forest']][[param_i]][['quantile025']])*4*Nref, as.numeric(posterior[['random_forest']][[param_i]][['expectation']])*4*Nref, as.numeric(posterior[['random_forest']][[param_i]][['quantile975']])*4*Nref, sep='\t'), outfile, append=T)
####	
####	} else{
####		write(paste(param_i, as.numeric(quantile(posterior[['neural_network']][,i], 0.025)), as.numeric(quantile(posterior[['neural_network']][,i], 0.5)), as.numeric(quantile(posterior[['neural_network']][,i], 0.975)), sep='\t'), outfile, append=T)
####		write(paste(param_i, as.numeric(posterior[['random_forest']][[param_i]][['quantile025']]), as.numeric(posterior[['random_forest']][[param_i]][['expectation']]), as.numeric(posterior[['random_forest']][[param_i]][['quantile975']]), sep='\t'), outfile, append=T)
####	}
#	}
# }

