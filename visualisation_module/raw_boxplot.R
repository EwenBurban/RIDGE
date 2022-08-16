#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
		y = vector()
		y[tmp[1]] = tmp[2]
		return(y)},USE.NAMES=F)
locus_data = read.table(args['locus_data'],h=T)
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
	divplot_list[[var]] = ggviolin(locus_data,NULL,var,add='boxplot') + geom_hline(yintercept=obs_dataset[,var],color='red')
}
output_list[['divergence']] = ggarrange(plotlist=divplot_list,ncol=2,nrow=2)
print('divergence_done')
## diversity plots ###
diversity_plot = list()
color_vec = get_palette('npg',2)

piA = ggviolin(locus_data,NULL,'piA_avg',add='boxplot',color=color_vec[1]) + geom_hline(yintercept=obs_dataset[,'piA_avg'],color='red')
piB = ggviolin(locus_data,NULL,'piB_avg',add='boxplot',color=color_vec[2]) + geom_hline(yintercept=obs_dataset[,'piB_avg'],color='red')
print('pi done')
thetaA = ggviolin(locus_data,NULL,'thetaA_avg',add='boxplot',color=color_vec[1]) + geom_hline(yintercept=obs_dataset[,'thetaA_avg'],color='red')
thetaB = ggviolin(locus_data,NULL,'thetaB_avg',add='boxplot',color=color_vec[2]) + geom_hline(yintercept=obs_dataset[,'thetaB_avg'],color='red')
print('theta done')
DtajA = ggviolin(locus_data,NULL,'DtajA_avg',add='boxplot',color=color_vec[1]) + geom_hline(yintercept=obs_dataset[,'DtajA_avg'],color='red')
DtajB = ggviolin(locus_data,NULL,'DtajB_avg',add='boxplot',color=color_vec[2]) + geom_hline(yintercept=obs_dataset[,'DtajB_avg'],color='red')
print('dtaj done')
output_list[['diversity']] = ggarrange(plotlist=list(piA,piB,thetaA,thetaB,DtajA,DtajB),nrow=3,ncol=2,align='v')

if(grepl(locus_data$dataset[1],pattern=":")){
	print('chromosom part')
	locus_data$delta_pi = locus_data$piB_avg - locus_data$piA_avg
	locus_data$delta_theta = locus_data$thetaB_avg - locus_data$thetaA_avg
	locus_data$delta_Dtaj = locus_data$DtajB_avg - locus_data$DtajA_avg

	pos = do.call(rbind,lapply(locus_data$dataset,get_pos))
	locus_data = cbind(pos,locus_data)
	all_var = c('FST_avg','divAB_avg','netDivAB_avg','sf_avg','delta_pi','delta_theta','delta_Dtaj')
	chrvar_list = list()
	for(var in all_var){
	chrvar_list[[var]] = ggviolin(locus_data,'chr',var,add='boxplot',palette='npg')
	}
	output_list[['chr_var']] = ggarrange(plotlist=chrvar_list,nrow=2,ncol=1)
}
print(output)
ggexport(plotlist=output_list,filename=output)
