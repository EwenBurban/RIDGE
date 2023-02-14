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
}



######################################################

options(digits=5)
list_models_param = c(bestModel)

get_posterior<-function(nameA='spA', nameB='spB', nSubdir=10, sub_dir_sim='iteration_2', model='estimation_1', sub_dir_model=1, nPosterior=1000, figure=FALSE, timeStamp='SI_2N', path2observation='SI_2N/iteration_2/pod_0', do_nnet=T, useSFS=0, ncores=ncores){
	library(data.table)
	options(digits=5)
	###################
	# get observed data
	# observed data
	obs_ss_tmp = read.table(paste(path2observation, '/ABCstat_global.txt', sep=''), h=T)
	obs_ss_tmp = obs_ss_tmp[,-1]
	# end of remove stats
	if(useSFS == 1){
		sfs_obs = read.table(paste(path2observation, '/ABCjsfs.txt', sep=''), h=T)
		sfs_obs = sfs_obs[, -c(which(colnames(sfs_obs)=='fA0_fB0'), which(colnames(sfs_obs)=='fA1_fB0'), which(colnames(sfs_obs)=='fA0_fB1'))] # remove the singletons and the (0,0)
		ss_obs = cbind(obs_ss_tmp, sfs_obs)
	}else{
		ss_obs = obs_ss_tmp
	}

	
	# read all these simulated datasets
	ss_sim = list()
	params_sim = list()
    ss_sim_tmp= list()
    sfs_sim_tmp= list()
    params_sim_tmp = list()

	for(rep in seq(0, nSubdir-1, 1)){
        tmp_ss = read.table(file.path(timeStamp,sub_dir_sim,paste0(model,'_',rep),'ABCstat_global.txt'),h=T)
        tmp_ss = tmp_ss[order(tmp_ss[,'dataset']),]
		na_lines = unique(unlist(lapply(tmp_ss,function(x) which(is.na(x)))))
		if(length(na_lines)> 0) { tmp_ss = tmp_ss[-na_lines,]}
        ss_sim_tmp[[length(ss_sim_tmp)+1]] = tmp_ss[,-1]



        if(useSFS==1){
            sfs_sim = read.table(file.path(timeStamp,sub_dir_sim,paste0(model,'_',rep),'ABCjsfs.txt'),h=T)
            sfs_sim = sfs_sim[order(sfs_sim[,'dataset']),]
            sfs_sim_tmp[[length(sfs_sim_tmp)+1]] = sfs_sim[,-1]
        }

        tmp_params = read.table(file.path(timeStamp,sub_dir_sim,paste0(model,'_',rep),'priorfile.txt'),h=T)
        tmp_params = tmp_params[order(tmp_params[,'dataset']),]
        tmp_params = tmp_params[match(tmp_ss$dataset,tmp_params$dataset),]
        params_sim_tmp[[length(params_sim_tmp)+1]] = tmp_params[,-1]

        stopifnot(nrow(tmp_ss)==nrow(tmp_params))
        if(useSFS==1){stopifnot(nrow(sfs_sim)==nrow(tmp_params))}

	}
    ss_sim_tmp = do.call(rbind,ss_sim_tmp)
    if(useSFS==1){    sfs_sim_tmp = do.call(rbind,sfs_sim_tmp)}
    params_sim_tmp = do.call(rbind,params_sim_tmp)
	
	# statistics
	if(useSFS == 1){
		ss_sim[[model]] = cbind(ss_sim_tmp, sfs_sim_tmp)  # with SFS
	}else{
		ss_sim[[model]] = ss_sim_tmp  # without SFS
	}

	# remove the uninformative statistics
	toRemove = 1
	for(i in 2:ncol(ss_sim[[model]])){
		sd_tmp = sd(as.numeric(c(ss_obs[i], ss_sim[[model]][,i])),na.rm=T)
		mean_tmp = mean(as.numeric(c(ss_obs[i], ss_sim[[model]][,i])),na.rm=T)
		if(is.na(sd_tmp/mean_tmp) | (sd_tmp/mean_tmp) <0.1  ){
			toRemove = c(toRemove, i)}
		}
	
	toRemove = unique(toRemove)
	ss_obs = ss_obs[-toRemove]
	ss_sim[[model]] = ss_sim[[model]][, -toRemove]

	# params
	params_sim[[model]] = params_sim_tmp
	nparams = ncol(params_sim[[model]])

	### write the prior
	nPrior = nrow(params_sim[[model]])
	
	##############
	# inferences
	target_rf = ss_obs 

#      RANDOM FOREST	
        library(abcrf)
        params_model_rf = as.matrix(params_sim[[model]]) # in case of debug
        stats_model_rf = ss_sim[[model]] # in case of debug

        rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
			qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
			} 
        posterior_list = list()
        for(i in 1:nparams){
        	parameter = params_model_rf[,i]
			param_name = colnames(params_sim[[model]])[i]
        	data = data.frame(parameter, stats_model_rf)
        	mod = regAbcrf(parameter~., data, ntree=1000)
			estimate = predict(mod, target_rf, data)
			posterior_list[[param_name]] = estimate$expectation 

		}
		posterior_dataset = do.call(cbind,posterior_list)
		posterior_dataset = cbind('dataset'=1:nrow(posterior_dataset),posterior_dataset)
		write.table(posterior_dataset,file=paste(output_dir,'/posterior_',
					sub_dir_model,'.txt', sep=''), row.names=F, col.names=T, sep='\t', quote=F)

	
	

}

get_posterior(nameA=nameA, nameB=nameB, nSubdir=nSubdir, sub_dir_sim=sub_dir_sim, model=bestModel, sub_dir_model=bestModel, nPosterior=nPosterior, figure=F, timeStamp=timeStamp, path2observation=path2observation, useSFS=useSFS, ncores=ncores)
