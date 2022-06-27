#!/usr/bin/env Rscript
library(ggpubr)
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
            y = vector()
            y[tmp[1]] = tmp[2]
            return(y)},USE.NAMES=F)

input_file=args['input_file']
output=args['output']
## load data
data = read.table(input_file,header=T,stringsAsFactors=F)

## dataset name manipulation
data$Chr = as.numeric(sub('Chr','',sub(':.*','',data$dataset))) # gather the chr number
data$mid_pos = sapply(data$dataset,function(x){ # gather the mid position of the locus
                          str_pos = unlist(strsplit(x,':'))[2]
                          listpos = unlist(strsplit(str_pos,'-'))
                          midpoint = mean(as.numeric(listpos))
                          return(midpoint)})

data = data[order(data$Chr,data$mid_pos),]
data$order= 1:nrow(data)
row.names(data)=data$order
## plot
variable2plot = c('bialsites_avg','FST_avg','piA_avg','piB_avg','thetaA_avg','thetaB_avg','netdivAB_avg')
listplot=list()
break_pos=tapply(data$order,data$Chr,mean)
print(break_pos)
data$Chr=as.factor(data$Chr)
for (variable in variable2plot){
    listplot[[length(listplot)+1]]= ggscatter(data,'order',y=variable,
                                    color="Chr",add='loess',
                                    xlab=F,show.legend.text=F,size=0.7
    ) + scale_x_continuous(breaks=break_pos,labels=levels(data$Chr)) + theme(legend.position='n') + geom_smooth(method='loess',span=0.01,color='black',se=F)
}

ggarrange(plotlist=listplot,ncol=1,nrow=4) %>% ggexport(filename = output,height=8)
