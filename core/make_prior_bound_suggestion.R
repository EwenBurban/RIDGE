#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
            y = vector()
            y[tmp[1]] = tmp[2]
            return(y)},USE.NAMES=F)
data=args['data']
output=args['output']
mu = as.numeric(args['mu'])
data=read.table(data,h=T)
popfile = read.table(args['popfile'],sep=',',h=T)
nameA = args['nameA']
nameB = args['nameB']
popA = popfile[,nameA]
popB = popfile[,nameB]
ploidy= as.numeric(args['ploidy'])
size_popA = length(popA[!is.na(popA)]) * ploidy
size_popB = length(popB[!is.na(popB)]) * ploidy
totpopsize=size_popA+size_popB
window_size = as.numeric(args['window_size'])
hetero_theta=args['hetero_theta']


##### Determine the bound of Ne #####
# theta = 4 * Ne * mu so Ne = theta/(4*mu)
# To obtain bound of prior, it take the minimum value of the quantile at 5% of Ne estimated across both population
# and the maximum value of the quantile at 95% of Ne estimated across both population
# Nref is the mean values of N_min and N_max
N_A=data$piA_avg / (4*mu)
N_B=data$piB_avg / (4*mu)
N_min = round(min(quantile(N_A,0.05),quantile(N_B,0.05)))
N_max = round(max(quantile(N_A,0.95),quantile(N_B,0.95)))
N_ref = round(mean(c(N_min,N_max)))
#### Determine mu from data if hetero_theta=True else mu is setup by user
# Using the watterson theta estimator, we obtain from … (TO-DO :finish commentary )
if(hetero_theta=='True'){
	mu_vec=data$bialsite_avg/(4*N_ref*sum(1/(1:(totpopsize-1)))*window_size)
	mu=mean(mu_vec,na.rm=T)
}
### Determine the bouds of Tsplit
# Da (net divergence; netDivAB) = 2 * mu * Tsplit so Tsplit=Da/(2*mu) ; this way tend to underestimate Tsplit
# Dxy (absolute divergence; divAB) = 2*mu*Tsplit + 4*Ne*mu ; this way tend to overestimate 
Tsplit_min=round(quantile(data$netDivAB_avg/(2*mu),0.05))
Tsplit_max=round(quantile((data$divAB_avg-4*N_ref*mu)/(2*mu),0.95))
### Determine the maximum proportion of Pbarrier /!\ ABANDONNED
# Pbarrier is estimated by the proportion of outlier for the following summary stat : Fst, Dxy, Da, Sf, ss, piA, piB
# To determine outlier, the tukey’s fence definition is applied
P_outlier=function(x,mode='supp'){ # calculate tukey fence and then compute the proportion of outliers 
	IQ = quantile(x,0.75)-quantile(x,0.25)
	if (mode=='supp'){fence=quantile(x,0.75) + 1.5*IQ; y=length(which(x>fence))/length(x)}
	if (mode=='inf'){fence=quantile(x,0.25) - 1.5*IQ; y=length(which(x<fence))/length(x)}
	return(y)}
P_max=max(P_outlier(data$FST_avg),P_outlier(data$divAB_avg),P_outlier(data$netDivAB_avg),P_outlier(data$sf_avg),
		  P_outlier(data$ss_avg,mode='inf'),P_outlier(data$piA_avg,mode='inf'),P_outlier(data$piB_avg,mode='inf'))


suggestion=data.frame(prior_bound=c('N_ref','N_min','N_max','Tsplit_min','Tsplit_max') , values=c(N_ref,N_min,N_max,Tsplit_min,Tsplit_max))
write.table(suggestion,file=output,row.names=F,col.names=T,quote=F)
