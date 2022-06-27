#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
		y = vector()
		y[tmp[1]] = tmp[2]
		return(y)},USE.NAMES=F)
gwscan_data = read.table(args['gwscan_data'],h=T)
obs_dataset = read.table(args['obs_data'],h=T)
output = args['output']
output_list = list()
library(ggpubr)

get_pos = function(str){
    tmp=unlist(strsplit(str,split=':'))
    chr=as.integer(sub('Chr','',tmp[1]))
    tmp_bis = unlist(strsplit(tmp[2],split='-'))
    start=as.integer(tmp_bis[1])
    end=as.integer(tmp_bis[2])
    res = c('chr'=chr,'start'=start,'end'=end)
    return(res)
}

divergence_var = c('FST_avg','divAB_avg','netDivAB_avg','sf_avg')
divplot_list = list()
for(var in divergence_var){
	divplot_list[[var]] = ggviolin(gwscan_data,NULL,var,add='boxplot') + geom_hline(yintercept=obs_dataset[,var],color='red')
}
output_list[['divergence']] = ggarrange(plotlist=divplot_list,ncol=2,nrow=2)
print('divergence_done')
## diversity plots ###
diversity_plot = list()
color_vec = get_palette('npg',2)

piA = ggviolin(gwscan_data,NULL,'piA_avg',add='boxplot',color=color_vec[1]) + geom_hline(yintercept=obs_dataset[,'piA_avg'],color='red')
piB = ggviolin(gwscan_data,NULL,'piB_avg',add='boxplot',color=color_vec[2]) + geom_hline(yintercept=obs_dataset[,'piB_avg'],color='red')
print('pi done')
thetaA = ggviolin(gwscan_data,NULL,'thetaA_avg',add='boxplot',color=color_vec[1]) + geom_hline(yintercept=obs_dataset[,'thetaA_avg'],color='red')
thetaB = ggviolin(gwscan_data,NULL,'thetaB_avg',add='boxplot',color=color_vec[2]) + geom_hline(yintercept=obs_dataset[,'thetaB_avg'],color='red')
print('theta done')
DtajA = ggviolin(gwscan_data,NULL,'DtajA_avg',add='boxplot',color=color_vec[1]) + geom_hline(yintercept=obs_dataset[,'DtajA_avg'],color='red')
DtajB = ggviolin(gwscan_data,NULL,'DtajB_avg',add='boxplot',color=color_vec[2]) + geom_hline(yintercept=obs_dataset[,'DtajB_avg'],color='red')
print('dtaj done')
output_list[['diversity']] = ggarrange(plotlist=list(piA,piB,thetaA,thetaB,DtajA,DtajB),nrow=3,ncol=2,align='v')

gwscan_data$delta_pi = gwscan_data$piB_avg - gwscan_data$piA_avg
gwscan_data$delta_theta = gwscan_data$thetaB_avg - gwscan_data$thetaA_avg
gwscan_data$delta_Dtaj = gwscan_data$DtajB_avg - gwscan_data$DtajA_avg
pos = do.call(rbind,lapply(gwscan_data$dataset,get_pos))
gwscan_data = cbind(pos,gwscan_data)
all_var = c('FST_avg','divAB_avg','netDivAB_avg','sf_avg','delta_pi','delta_theta','delta_Dtaj')
chrvar_list = list()
for(var in all_var){
chrvar_list[[var]] = ggviolin(gwscan_data,'chr',var,add='boxplot',palette='npg')
}
output_list[['chr_var']] = ggarrange(plotlist=chrvar_list,nrow=2,ncol=1)

ggexport(plotlist=output_list,filename=output)
