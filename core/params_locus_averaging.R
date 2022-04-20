#!/usr/bin/env Rscript
library(abcrf)
library(plyr)
mode_locus = F
for(arg in commandArgs()){
  tmp = strsplit(arg, '='); opt = tmp[[1]][1]; val = tmp[[1]][2]
  if(opt == "obs_dir"){obs_dir=val}
  if(opt == "timeStamp"){timeStamp=val}
  if(opt == "param_dir"){param_dir=val}
  if(opt == "weight_pattern"){weight_pattern=val}
  if(opt == "output"){output=val}
  if(opt == 'weight'){weight_file=val}
  if(opt == 'param_pattern'){param_pattern=val}
}

params_model = c('N1','N2','Na','founders1','founders2','Tdem1','Tdem2','M12_vec','M21_vec','shape_N_a','shape_N_b','shape_M12_a','shape_M12_b','shape_M21_a','shape_M21_b','Tsc','Tam','Tsplit')
MODELS_COMP = c('SC_1M_1N', 'SC_1M_2N', 'SC_2M_1N', 'SC_2M_2N', 'AM_1M_1N', 'AM_1M_2N', 'AM_2M_1N', 'AM_2M_2N', 'IM_1M_1N', 'IM_1M_2N', 'IM_2M_1N', 'IM_2M_2N', 'SI_1N', 'SI_2N')
### Param estimation ###
model_weight = read.table(weight_file,h=T)
models_used = colnames(model_weight)
list_model_param_estimation = list.files(path=param_dir,pattern=param_pattern)
param_estimation = list()
for(m in models_used){
	concerned_file =file.path(param_dir, grep(pattern=m,list_model_param_estimation,value=T))
	param_estimation[[m]] = do.call(cbind,lapply(concerned_file,read.table,h=T)) * as.numeric(model_weight[1,m]) 
}
votes_Wparams = Reduce('+',param_estimation)

col = params_model[match(colnames(votes_Wparams),params_model,nomatch=F)]

write.table(votes_Wparams,file=output,row.names=F,col.names=T,quote=F,sep='\t')
