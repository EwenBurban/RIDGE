#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
		y = vector()
		y[tmp[1]] = tmp[2]
		return(y)},USE.NAMES=F)
rho_map = read.table(args['rho_map'],sep='\t',h=T)
Nref = as.numeric(args['Nref'])
window_size = as.numeric(args['window_size'])
mu = as.numeric(args['mu'])
popfile = read.table(args['popfile'],sep=',',h=T)
nameA = args['nameA']
nameB = args['nameB']
bedfile = read.table(args['bedfile'],sep='\t',h=T)
popA = popfile[,nameA]
popB = popfile[,nameB]
size_popA = length(popA[!is.na(popA)])
size_popB = length(popB[!is.na(popB)])
output = args['output']
data = merge(bedfile,rho_map,by=c('chr','start','end'),all=F)
head(data)
colnames(data) = c('chr','start','end','rec_rate')
data$rho = data$rec_rate * 4 * Nref * window_size
data[is.na(data$rho),'rho'] = 0
data$theta = rep(mu * 4 * Nref * window_size,nrow(data))
data$locus_length = rep(window_size,nrow(data))
data$size_popA  = rep(size_popA,nrow(data))
data$size_popB  = rep(size_popB,nrow(data))
data$totpopsize = data$size_popA + data$size_popB
res = subset(data,select = c('locus_length','size_popA','size_popB','totpopsize','theta','rho'))
write.table(res,file=output,sep='\t',row.names=F)



