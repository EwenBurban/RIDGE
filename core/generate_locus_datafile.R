#!/usr/bin/env Rscript

# Goal of the script
# This script aim to produce a file defining the genomic properties of each locus defined in {bedfile}
# This information are mandatory to simulate dataset containing locus with the same properties than observed dataset.

### gather command args and define variables
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
ploidy= as.numeric(args['ploidy'])
size_popA = length(popA[!is.na(popA)]) * ploidy
size_popB = length(popB[!is.na(popB)]) * ploidy
output = args['output']

# generate a data.frame named data, containing all raw information about genomic windows selected
data = merge(bedfile,rho_map,by=c('chr','start','end'),all=F)
head(data)
colnames(data) = c('chr','start','end','rec_rate')
data$rho = data$rec_rate * 4 * Nref * window_size # calculate the rho value from local r 
data[is.na(data$rho) | data$rho < 0,'rho'] = 0 # avoid NA or aberrant values by replacing them by 0
data$theta = rep(mu * 4 * Nref * window_size,nrow(data)) # calculate theta value from {mu}
data$locus_length = rep(window_size,nrow(data)) 

#generate  vectors of the same size than data, assuming that number of sample available is the same at each loci 
data$size_popA  = rep(size_popA,nrow(data))# generate a vector of the same size than data, assuming that number of sample available is the same at each loci 
data$size_popB  = rep(size_popB,nrow(data))
data$totpopsize = data$size_popA + data$size_popB

# It keep only the necessary information for the rest of the pipeline and store it in {output} file
res = subset(data,select = c('locus_length','size_popA','size_popB','totpopsize','theta','rho'))
write.table(res,file=output,sep='\t',row.names=F)



