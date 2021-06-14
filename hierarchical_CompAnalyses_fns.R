# Functions to perform a hierarchical analyses of ~fMRI data,
#   using either hierarchical ("treelets") PCA or hierarchical Independent Components Analysis.
#
# Procedure:
#   1. Compress data using PCA applied individually to each subject,
#       w/ compression applied to temporal dimension
#   2. Compress data using PCA applied to group data,
#       w/ subjs. concatenated in (PCA-reduced) temporal dimension
#   3. Enter group data into multi-model order ICA,
#       applied to spatial dimension, estimate group-level sources
#   4. Estimate hierarchy, as hPCA or hICA, on multi-model order sources,
#       w/ statistics used to identify redundant sources in initial levels of hierarchy
#   (4.1. Optionally identify & discard artifactual, nuisance and noise sources
#       before continuing to estimate higher levels of the hierarchy)
#   5. Back-reconstruct spatial maps & time series for all sources,
#       specific to each subject
#   6. Calculate descriptive statistics (e.g. group-mean spatial maps) 
#
#
# References:
#   Lee, A.B., Nadler, B., and Wasserman, L. (2008).  Treelets--An adaptive
#     multi-scale basis for sparse unordered data.  The Annals of Applied
#     Statistics 2(2), p. 435-471.
#   Secchi, P., Vantini, S., and Zanini, P.  (2016).  Hierarchical independent
#     component analysis: a multi-resolution non-orthogonal data-driven basis.
#     Computational Statistics and Data Analysis 95, p. 133-149.
#   Hyvarinen, A., and Oja, E.  (2000).  Independent Component Analysis: Algorithms
#     and Applications.  Neural Networks, 13(4-5), p. 411-430.
#

stopifnot(exists('PCA.Rachakonda2016'))  # required main fn. from 'PCA_fns.R'


##################################################################
detrend <- function(ts, detrend_rows=T){
  ##################################################################
  # Remove linear trend from time series'
  #    (code from stats::spec.pgram)
  #
  # Inputs:
  #   ts : time series matrix
  #   detrend_rows : detrend rows of matrix if T, otherwise detrend cols. 
  # Output:
  #   ts : as above, after subtracting global mean at each time point
  #
  
  if (detrend_rows){ ts = t(ts) }
  
  N <- nrow(ts)
  t <- 1L:N - (N + 1)/2
  sumt2 <- N * (N^2 - 1)/12
  for (i in 1L:ncol(ts)) ts[, i] <- ts[, i] - mean(ts[, i], na.rm=T) - 
    sum(ts[, i] * t, na.rm=T) * t/sumt2
  
  if (detrend_rows){ ts = t(ts) }
  return(ts)
} ##################################################################


##################################################################
subj_preprocData <- function(data, 
                             preproc.type='rm.voxel.mean', 
                             discard.n=NA,
                             var=NA, mask=NA, 
                             verbose=F){
  ##################################################################
  # Pre-processing data prior to subj.-level PCA,
  #   called as part of subj_PCA() & groupICA_BackReconstruction().
  #     Recommended for PCA prior to & part of ICA
  #
  # Inputs:
  #   data : (options)
  #           1. path to subj. file, w/ data 'var' as [voxels x time] matrix
  #           2. matrix w/ dims. [voxels x time] matrices
  #   preproc.type : pre-processing type (options, or skip if NA)
  #       rm.voxel.mean  : center all voxel time series' to mean = 0
  #       rm.time.mean   : center each timepoint (~vol.) to mean = 0
  #       variance.norm  : detrend all voxels & scale to standard deviation
  #   var  : if 'data' w/ option 1, name of saved variable w/ subj. data
  #   mask : if 'data' is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  #   discard.n : if 'data' is a Nifti file(s), 
  #                 number of initial vols. to discard to allow for magnet stabilization
  # Output:
  #   data  : [voxels x time] matrix pre-processed as above
  # Requires:
  #   load.data() : data munging fn. from PCA_fns.R
  #   detrend()   : detrending for during variance normalization
  #
  
  if (verbose){  cat('\n......pre-processing subj. data:  ')}
  
  preproc.type = as.character(preproc.type)
  discard.n = as.integer(discard.n)
  data = load.data(data, var, mask)
  
  if (!all(is.na(discard.n)) && discard.n > 0){
    if (verbose){  cat(paste0('discarding initial ', discard.n, ' volumes...'))}
    n1 = discard.n+1
    p = ncol(data)
    data = data[,n1:p]
  }
  
  if (all(is.na(preproc.type))){
    return(data)
    
  }else if (!(preproc.type %in% c('rm.voxel.mean', 'rm.time.mean', 'variance.norm'))){
    if (verbose){
      cat('\n  WARNING: could not understand preprocessing type,')
      cat(' defaulting to removing voxel means')
    }
    return(sweep(data, 1, rowMeans(data)))
    
  }else{
    if (preproc.type == 'rm.voxel.mean'){
      if (verbose){  cat('removing voxel means')  }
      return(sweep(data, 1, rowMeans(data)))
      
    }else if (preproc.type == 'rm.time.mean'){
      if (verbose){  cat('removing timepoint means')  }
      return(sweep(data, 2, colMeans(data)))
      
    }else if (preproc.type == 'variance.norm'){
      if (verbose){  cat('detrending voxels & scaling to s.d.')  }
      data = detrend(data, detrend_rows=T)
      data = t(apply(data, 1, function(yy) return(yy / sd(yy))))
      return(data)
    }
  }
} ##################################################################


##################################################################
subj_PCA <- function(data, K1=NA, 
                     preproc.type='rm.voxel.mean',
                     discard.n=NA, whitenPCs=T,
                     data_wc='*.nii', var=NA, mask=NA,
                     prefix=NA, save_path=NA,
                     return.data=F, parallel=F, verbose=F){
  ##################################################################
  # Performances subject-level dimension reduction w/ PCA
  #
  # Inputs:
  #   data : (options)
  #           1. vector of subj. files, w/ data 'var' as [voxels x time] matrix
  #           2. path to directory containing above
  #           3. list of [voxels x time] matrices
  #           4. array w/ dims. [subjs. x voxels x time]
  #   K1   : PCA subspace dimension,
  #           if NA defaults to min. dim. in data
  #   preproc.type : pre-processing type  (options, or skip if NA)
  #       rm.voxel.mean  : center all voxel time series' to mean = 0
  #       rm.time.mean   : center each timepoint (~vol.) to mean = 0
  #       variance.norm  : detrend all voxels & scale to standard deviation
  #   discard.n : if 'data' is a Nifti file(s), 
  #                 number of initial vols. to discard to allow for magnet stabilization
  #   whitenPCs : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  #   data_wc   : if 'data' w/ option 2, wildcard to filter subj. files
  #   var       : if 'data' w/ option 1, name of saved variable w/ subj. data
  #   mask      : if 'data' is a Nifti file(s), optional binary mask
  #   prefix    : identifying prefix for experiment, 
  #                 attached to saved output filename
  #   save_path   : path to dir. for saved files
  #   return.data : return list of subj. PCA spaces if T,
  #                   or path to dir. w/ saved files if F
  # Output: 
  #   saved files w/ *_pca1.RData added as suffix & vars.:
  #     pca$X  : [v by K1] matrix of subj. PCA space
  #     pca$U  : [Mp x K1] reducing matrix of eigenvectors
  #     pca$Lambda : eigenvalues
  # Requires:
  #   subj_preprocData()   : pre-processes data
  #   load.data()          : data munging from PCA_fns.R
  #   detrend()            : detrending for variance normalization
  #   PCA.Rachakonda2016() : main fn. in PCA_fns.R
  #
  
  if (verbose){  cat('\nCalculating subj.-level PCA spaces...')}
  if (parallel){  library(doParallel)  }
  
  ### Defaults ###
  subj_pca_suffix  = '_PCA1.RData'       # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData'    # filesuffix for group-level temporal PCA from group_PCA()
  group_ica_suffix = 'spatialICA.RData'  # filesuffix for group spatial ICA from group_spatialICA()
  br_suffix        = '_ICAsources.RData' # filesuffix for back-reconstructed subj. sources
  
  pca_algorithm = NA  # Initialize computation settings for PCA...
  pca_unstacked = NA  # ...as determined by 1st subj. & continued thereafter
  
  ### Inputs ###
  if (is.character(data) && (length(data)==1) && dir.exists(data)){
    stopifnot(length(list.files(data)) > 0)
    data = list.files(data, full.names=T)
    data = data[!dir.exists(data)] # exclude subdirs.
    
    if (!all(is.na(data_wc))){     # filter files
      data = data[grepl(data_wc, data)]
    }
  }
  if (is.vector(data) | is.list(data)){
    S = length(data)
  }else if (is.array(data)){
    S = dim(data)[1]
  }else{
    print('WARNING: Could not understand input data!')
    S = length(data)
  }
  
  if (!all(is.na(var)) && is.character(data) 
      && !any(grepl('.RData$', data))){
    var = NA   # only applicable for data in saved R files
  }else if (!all(is.na(mask)) && is.character(mask) 
            && !any(grepl('.nii$', data))){
    mask = NA  # only applicable for paths to nifti vols.
  }
  
  K1 = as.integer(K1)
  whitenedPCs = whitenPCs = as.logical(whitenPCs)
  prefix = as.character(prefix)
  save_path = as.character(save_path)
  if (return.data){  PCA_s = list() }
  
  
  ### Main fn. ###
  if (parallel){
    if (verbose){
      for (s in 1:S){  cat(paste0('\n...PCA compressing (in parallel):  ', data[s]))  }
    }
    
    ### Get best algorithm & run all subjs. w/ same settings ###
    recommendations = find_PCAtype_best(data[1], k=K1, 
                                        var=var, mask=mask,
                                        verbose=verbose)
    pca_algorithm = recommendations$PCAtype
    pca_unstacked = recommendations$unstacked
    
    PCA_s = foreach(s=icount(S),
                    .export=c('load.data', 'zeroMean_Yi', 'whitenedPCs',
                              'detrend', 'subj_preprocData')) %dopar% {
                                if (is.character(data) && all(file.exists(data))){
                                  ### Create save file name ###
                                  data_files = data[s]  # plural for consistency w/ other fns.
                                  save.fname = basename(data_files)
                                  if (!is.na(prefix)){
                                    save.fname = paste0(prefix, save.fname)
                                  }
                                  save.fname = sub('.nii$', '', save.fname)
                                  save.fname = sub('.RData$', '', save.fname)
                                  save.fname = paste0(save.fname, subj_pca_suffix)
                                  
                                  ### Load & pre-process data ###
                                  data.preproc = subj_preprocData(data[s], 
                                                                  preproc.type, discard.n,
                                                                  var=var, mask=mask,
                                                                  verbose=verbose)
                                  
                                }else if (is.list(data)){
                                  ### Create save file name ###
                                  data_files = NA
                                  save.fname = paste0('s',s,subj_pca_suffix)
                                  if (!is.na(prefix)){
                                    save.fname = paste0(prefix, save.fname)
                                  }
                                  
                                  ### Load & pre-process data ###
                                  data.preproc = subj_preprocData(data[[s]], 
                                                                  preproc.type, discard.n,
                                                                  var=var, mask=mask,
                                                                  verbose=verbose)
                                }
                                
                                
                                ### Main fn.:  Subj-level PCA ###
                                pca = PCA.Rachakonda2016(data.preproc, K1, 
                                                         var=var, mask=mask,
                                                         whitenPCs=whitenPCs,
                                                         algorithm=pca_algorithm,
                                                         unstacked=pca_unstacked,
                                                         verbose=verbose)
                                
                                
                                ### Saving ###
                                if (!is.na(save_path)){
                                  save.fname = file.path(save_path, save.fname)
                                }else{
                                  save.fname = file.path(dirname(data_files), save.fname)
                                }
                                if (verbose){  cat(paste0('\n......saving as: ', save.fname))}
                                save(list=c('pca', 'K1', 'preproc.type', 'discard.n', 'whitenedPCs',
                                            'data_files', 'var', 'mask', 'prefix', 's'), 
                                     file=save.fname)
                                stopifnot(file.exists(save.fname))
                                if (return.data){ return(pca) }
                              }
    
  }else{
    for (s in 1:S){
      if (verbose){  cat(paste0('\n...PCA compressing:  ', data[s]))  }
      
      if (is.character(data) && all(file.exists(data))){
        ### Create save file name ###
        data_files = data[s]
        save.fname = basename(data_files)
        if (!is.na(prefix)){
          save.fname = paste0(prefix, save.fname)
        }
        save.fname = sub('.nii$', '', save.fname)
        save.fname = sub('.RData$', '', save.fname)
        save.fname = paste0(save.fname, subj_pca_suffix)
        
        ### Load & pre-process data ###
        data.preproc = subj_preprocData(data[s], 
                                        preproc.type, discard.n,
                                        var=var, mask=mask,
                                        verbose=verbose)
        
      }else if (is.list(data)){
        ### Create save file name ###
        data_files = NA
        save.fname = paste0('s', s, subj_pca_suffix)
        if (!is.na(prefix)){
          save.fname = paste0(prefix, save.fname)
        }
        
        ### Load & pre-process data ###
        data.preproc = subj_preprocData(data[[s]], 
                                        preproc.type, discard.n,
                                        var=var, mask=mask,
                                        verbose=verbose)
      }
      
      
      ### Main fn.:  Subj-level PCA ###
      #     Automatically selects PCA algorithm, 
      #       will run unstacked in parallel, if nodes available
      
      pca = PCA.Rachakonda2016(data.preproc, K1, 
                               var=var, mask=mask,
                               whitenPCs=whitenPCs,
                               algorithm=pca_algorithm,
                               unstacked=pca_unstacked,
                               verbose=verbose)

      ### Run all subsequent subjs. w/ same settings ###
      pca_algorithm = pca$algorithm  # initialized to NA above...
      pca_unstacked = pca$unstacked  # ...then determined by 1st subj. & continued thereafter
      
      
      ### Saving ###
      if (!is.na(save_path)){
        save.fname = file.path(save_path, save.fname)
      }else{
        save.fname = file.path(dirname(data_files), save.fname)
      }
      if (verbose){  cat(paste0('\n......& saving as: ', save.fname, '\n'))}
      save(list=c('pca', 'K1', 'preproc.type', 'discard.n', 'whitenedPCs',
                  'data_files', 'var', 'mask', 'prefix', 's'), 
           file=save.fname)
      stopifnot(file.exists(save.fname))
      if (return.data){ PCA_s[[s]] = pca  }
    }
  }
  if (verbose){  cat('\n')  }
  
  ##################################################################
  if (return.data){
    return(PCA_s)     # list of subj.-specific PCAs
  }
} ##################################################################


##################################################################
get_subjPCA_dat <- function(subj_pca_file, prefix=NA, 
                            return.PCs=T, return.orig.data=F,
                            verbose=F){
  ##################################################################
  # Loads & formats subj-level PCA data,
  #   & applies PCA pre-processing to original data as output of subj_PCA()
  #
  # Requires:
  #   subj_preprocData()   : pre-processes data
  #   load.data()          : data munging fn. from PCA_fns.R
  #   detrend()            : detrending for pre-processing
  #
  
  if (verbose){cat(paste0('\n...loading subj.-level PCA from:\n      ',subj_pca_file))}
  
  Space.subj = new.env()
  load(subj_pca_file, envir=Space.subj)
  stopifnot(all(c('pca', 'K1', 
                  'data_files', 'var', 'prefix') %in% 
                  names(Space.subj)))
  
  if (all(is.na(prefix))){
    prefix = get('prefix', envir=Space.subj)
  }else{
    stopifnot(prefix == get('prefix', envir=Space.subj))
  }
  subj_data_file = get('data_files', envir=Space.subj)
  stopifnot(is.character(subj_data_file) 
            && file.exists(subj_data_file))
  var = get('var', envir=Space.subj)
  mask = get('mask', envir=Space.subj)
  preproc.type = get('preproc.type', envir=Space.subj)
  discard.n = get('discard.n', envir=Space.subj)
  subjPCA = get('pca', envir=Space.subj)
  
  subjWhitenedPCs = get('whitenedPCs', envir=Space.subj)
  if (!subjWhitenedPCs){  subjPCA$Lambda = NA  }
  
  ### Apply PCA pre-processing to data ###
  if (return.orig.data){
    data.preproc = subj_preprocData(subj_data_file, 
                                    preproc.type, discard.n,
                                    var=var, mask=mask,
                                    verbose=verbose)
    data.preproc = zeroMean_Yi(data.preproc)
  }else{  data.preproc = NA  }
  
  return.PCs = as.logical(return.PCs)
  if (!return.PCs){  subjPCA$X = NA  }
  
  
  ##################################################################
  return(list('Y'   = data.preproc,     # pre-processed subj. data used for PCA
              'X'   = subjPCA$X,        # subj. Principal Components
              'U_i' = subjPCA$U,        # subj.-level PCA eigenvectors
              'L_i' = subjPCA$Lambda,   # subj.-evel PCA eigenvalues
              'subj_pca_file' = subj_pca_file,
              'prefix'        = prefix))
} ##################################################################


##################################################################
group_PCA <- function(data, K2=NA, 
                      whitenPCs=T,
                      prefix=NA, save_path=NA,
                      return.data=F, verbose=F){
  ##################################################################
  # Performances group-level dimension reduction w/ PCA
  #
  # Inputs:
  #   data : (options)
  #           1. vector of subj. PCA files, output from subj_PCA()
  #           2. path to directory containing above
  #           3. list of [voxels x time subspace] matrices
  #           4. array w/ dims. [subjs. x voxels x time subspace]
  #   K2   : reduced group subspace dimension,
  #           if NA defaults to min. dim. in data
  #   whitenPCs   : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  #   prefix      : identifying prefix for experiment, 
  #                  attached to saved output filename
  #   save_path   : path to dir. for saved files
  # Output: 
  #   saved files w/ *_groupPCA.RData added as suffix & vars.:
  #     pca$X  : [v by K2] matrix of group PCA space
  #     pca$U  : [K2 x concatenated time] reducing matrix of eigenvectors
  #     pca$Lambda : eigenvalues
  # Requires:
  #   get_subjPCA_dat()    :  loads subj.-level PCs
  #   PCA.Rachakonda2016() :  main PCA fn. from PCA_fns.R
  #
  
  if (verbose){  cat('\nCalculating group-level PCA space...')}
  
  ### Inputs & defaults ###
  subj_pca_suffix  = '_PCA1.RData'       # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData'    # filesuffix for group-level temporal PCA from group_PCA()
  group_ica_suffix = 'spatialICA.RData'  # filesuffix for group spatial ICA from group_spatialICA()
  br_suffix        = '_ICAsources.RData' # filesuffix for back-reconstructed subj. sources  
  
  prefix = as.character(prefix)
  whitenedPCs = as.logical(whitenPCs)

  ### Prelims ###
  subj_pca_wc = paste0('*', subj_pca_suffix)
  
  if (is.character(data) && all(file.exists(data))){
    if (dir.exists(data)){
      if (is.na(save_path)){
        save_path = data
      }
      data_files = list.files(data, subj_pca_wc, full.names=T)
    }else{
      data_files = data
      if (is.na(save_path)){
        save_path = file.path(unique(dirname(data_files))[1])
      }
    }
    
    data = vector('list', length(data_files))
    if (verbose){
      cat('\n...loading subj.-level PCA data from:')
      cat(paste0('\n      ', file.path(unique(dirname(data_files))[1]), subj_pca_wc))
    }
    for (f in 1:length(data_files)){
      data[[f]] = get_subjPCA_dat(data_files[f], prefix=prefix, 
                                  return.PCs=T, verbose=F)$X
    }

  }else{
    data_files = NA
  }
  
  ### Main fn. ###
  pca = PCA.Rachakonda2016(data, K2,             # Automatically selects PCA algorithm,
                           whitenPCs=whitenPCs,  #    will run in parallel if available
                           verbose=verbose)
  
  ### Saving ###
  save.fname = group_pca_suffix
  if (!is.na(prefix)){
    save.fname = paste0(prefix, '_', save.fname)
  }
  if (!is.na(save_path)){
    save.fname = file.path(save_path, save.fname)
  }
  if (verbose){  cat(paste0('\n......saving as: ', save.fname))}
  save(list=c('pca', 'data_files', 'K2', 
              'whitenedPCs', 'prefix'), file=save.fname)
  stopifnot(file.exists(save.fname))

  if (verbose){  cat('\n\n')  }
  ##################################################################
  if (return.data){  return(pca)  } # group PCA output
} ##################################################################


##################################################################
get_groupPCA_dat <- function(group_pca_file, prefix=NA, 
                             model_orders=NA, return.PCs=T,
                             verbose=F){
  ##################################################################
  # Loads & formats group-level PCA,
  #   output by group_PCA(),
  #     prior to ICA back-reconstruction
  #
  
  if (verbose){cat(paste0('\n...loading group-level PCA data from:\n      ',
                          group_pca_file))}
  
  return.PCs = as.logical(return.PCs)
  if (all(is.na(return.PCs))){  return.PCs = F  }
  
  Space.pca = new.env()
  stopifnot(is.character(group_pca_file) 
            && file.exists(group_pca_file))
  load(group_pca_file, envir=Space.pca)
  stopifnot(all(c('pca', 'data_files', 
                  'K2', 'prefix') %in% names(Space.pca)))
  if (all(is.na(prefix))){
    prefix = get('prefix', envir=Space.pca)
  }else{
    stopifnot(prefix == get('prefix', envir=Space.pca))
  }
  subj_pca_files = get('data_files', envir=Space.pca)
  groupPCA = get('pca', envir=Space.pca)
  groupWhitenedPCs = get('whitenedPCs', envir=Space.pca)
  if (!groupWhitenedPCs){
    groupPCA$Lambda = NA
  }
  K2 = get('K2', envir=Space.pca)
  
  ### Concat subspaces for multi-model order ICA ###
  model_orders = as.integer(model_orders)
  if (!all(is.na(model_orders))){ 
    G = NULL
    L = NULL
    if (return.PCs){  X = NULL  }else{  X = NA  }
    for (k in model_orders){
      G = cbind(G, groupPCA$U[,1:k], deparse.level=0)
      L = c(L, groupPCA$Lambda[1:k])
      if (return.PCs){
        X = cbind(X, groupPCA$X[,1:k], deparse.level=0)
      }
    }
  }else{
    G = groupPCA$U
    L = groupPCA$Lambda
    if (return.PCs){  X = groupPCA$X  }else{  X = NA  }
  }
  
  ##################################################################
  return(list('X' = X,          # group-level Principal Components, in cols.
              'G' = G,          # matrix of eigenvectors for group PCA subspace(s)
              'L' = L,          # group PCA eigenvalues for subspace(s)
              'num_comps' = K2, # overall dimension of group PCA subspace
              'model_orders' = model_orders,     # dims. of concatenated subspaces
              'subj_pca_files' = subj_pca_files, # subject PCA files used in group PCA
              'prefix'         = prefix))        # file prefix for experiment
} ##################################################################



##################################################################
group_spatialICA <- function(data, K=NA, 
                             prefix=NA, save_path=NA,
                             return.data=F, verbose=F){
  ##################################################################
  # Performs spatial Independent Components Analysis,
  #   on PCA compressed group-level data
  #   consisting of temporarily-concatenated & PCA compressed subjs.
  #
  # Inputs:
  #   data : (options)
  #           1. file w/ group PCA space, output from group_PCA()
  #           2. path to directory containing above
  #           3. [voxels x time subspace] matrix of reduced dim. group data
  #   K   : Number of sources to estimate, 
  #           either integer, or vector for multi-model order ICA,
  #           if NA defaults to min. dim. in data
  #   prefix        : identifying prefix to attach to saved output
  #   save_path     : path to dir. for saved files
  # Output: 
  #   saved file w/ *_groupICA.RData added as suffix & var. for model m:
  #     ica[[m]]$X : original data, pre-processed
  #     ica[[m]]$K : pre-ICA whitening matrix        [~time x reduced dim.]
  #     ica[[m]]$W : est. unmixing matrix            [reduced dim. x sources]
  #     ica[[m]]$A : mixing matrix, inverse of W     [sources x reduced dim.]
  #     ica[[m]]$S : source matrix                   [voxels x sources]
  #     ica[[m]]$skew_sign : sign of original skew, from sources in fastICA()
  # Requires:
  #   get_groupPCA_dat()        :  load  group-level PCA info
  #   fastICA()                 :  fast ICA algorithm from fastICA package
  #   align_ICAsources_toSkew() : changes sign of ICs to give p.d.f. a pos. skew
  #
  
  if (verbose){  cat('\nCalculating temporally-concatenated group spatial ICA...')}
  library(fastICA)
  
  ### Inputs & defaults ###
  K = as.integer(K)
  prefix = as.character(prefix)
  save_path = as.character(save_path)
  
  subj_pca_suffix  = '_PCA1.RData'       # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData'    # filesuffix for group-level temporal PCA from group_PCA()
  group_ica_suffix = 'spatialICA.RData'  # filesuffix for group spatial ICA from group_spatialICA()
  br_suffix        = '_ICAsources.RData' # filesuffix for back-reconstructed subj. sources  
  
  group_pca_wc = paste0('*', group_pca_suffix)
  
  if (is.character(data) && all(file.exists(data))){
    data_files = data
    
    if ((length(data_files)==1) && dir.exists(data_files)){
      stopifnot(length(list.files(data_files, group_pca_wc)) == 1)
      data_files = list.files(data_files, group_pca_wc, full.names=T)[1]
    }else{
      stopifnot(length(data_files) == 1)
    }
    
    data = get_groupPCA_dat(data_files, prefix,
                            return.PCs=T,
                            verbose=verbose)$X
    
    if (is.na(save_path)){
      save_path = unique(dirname(data_files))[1]
    }
    
  }else{  data_files = NA  }
  
  ### Sanity checks ###
  min.K = 2
  max.K = min(dim(data))
  if (all(is.na(K))){
    cat(paste0('   ICA model order(s) not specified, defaulting to K = ', max.K))
    K = max.K
  }
  if (any(K < min.K)){
    cat(paste0('   ICA model order(s) requested below min. possible, limiting analysis to K >=', min.K))
    K = K[K >= min.K]
  }
  if (any(K > max.K)){
    cat(paste0('   ICA model order(s) requested above max. possible, limiting analysis to K <=', max.K))
    K = K[K <= max.K]
  }
  
  
  ### Main fn. ###
  ica = list()
  m = 0
  for (k in K){
    m = m + 1
    ica[[m]] = list()
    
    if (verbose){
      cat(paste0('\n...calculating ICA using FastICA algorithm of Hyvarinen et al. (2000) w/ model order  K=',k,'...\n'))
      ica[[m]] = fastICA(data[,1:k], n.comp=k, verbose=T)
    }else{
      ica[[m]] = fastICA(data[,1:k], n.comp=k, method='C')  # compiled C++ fastICA
    }
    
    if (m > 1){  ica[[m]]$X = NULL  }  # skip redundant storage for multi-model order ICA
    ica[[m]]$model_order = k
  }
  
  ### Saving ###
  save.fname = group_ica_suffix
  if (!is.na(prefix)){
    save.fname = paste0(prefix, '_', save.fname)
  }else{
    save.fhame = sub('^_', '', save.fname)
  }
  if (!is.na(save_path)){
    save.fname = file.path(save_path, save.fname)
  }
  if (verbose){  cat(paste0('\n......& saving as: ', save.fname, '\n\n'))}
  save(list=c('ica', 'data_files', 'K', 'prefix'), file=save.fname)
  stopifnot(file.exists(save.fname))

  
  ##################################################################
  if (return.data){  return(ica)  } # group spatial ICA output
} ##################################################################

##################################################################
get_groupICA_dat <- function(ica_file, prefix=NA, 
                             align.sources=F,
                             flatten.multimodel.ICA=F,
                             flatten.source.est=F,
                             return.sources=T,
                             return.input.data=F,
                             verbose=F){
  ##################################################################
  # Loads & formats group-level spatial ICA output by group_spatialICA(),
  #   prior to subj.-level ICA back-reconstruction
  #     or hierarchical analyses
  #
  # Inputs:
  #  ica_file   : path to saved ICA file
  #  prefix     : experiment prefix
  #  align.sources : align all ICA source signs to a positive 3rd moment (skew)
  #  flatten.multimodel.ICA : flatten & format data from multiple ICA models
  #  flatten.source.est     : combine source estimation steps for convenience, 
  #                             by multiplying ICA pre-whitening & unmixing matrices
  #  return.sources         : if False, skip returning est. sources to minimize mem. use
  #  return.input.data      : if False (default), do not return data input to ICA
  # Output:  list(if flatten.multimodel.ICA) or list of lists w/ elements:
  #   'X'  : data input into ICA [voxels x ~time]
  #   'K'  : pre-whitening matrix      [~time x reduced dim.]
  #   'W'  : ICA unmixing matrices     [reduced dim. x sources]
  #   'KW' : matrix product of K and W, if flatten.source.est
  #   'A'  : ICA mixing matrices       [sources x reduced dim.]
  #   'S'  : ICA source matrix         [voxels x sources]
  #   'skew_signs'     : original skew signs of ICA sources
  #   'model_orders'   : ICA model orders
  #   'group_pca_file' : path to data used for ICA
  #   'prefix'         : file prefix for experiment
  #
  
  if (verbose){cat(paste0('\n...loading ICA data from:\n      ',ica_file))}
  
  Space.ica = new.env()
  load(ica_file, envir=Space.ica)
  stopifnot(all(c('ica', 'K', 'data_files', 'prefix') %in% names(Space.ica)))
  ica = get('ica', envir=Space.ica)
  
  if (all(c('K', 'W', 'A') %in% names(ica))){
    multimodel.ica = F
    flatten.multimodel.ica = F
  }else{
    multimodel.ica = T
    stopifnot(all(c('K', 'W', 'A') %in% names(ica[[1]])))
  }
  
  if (all(is.na(prefix))){
    prefix = get('prefix', envir=Space.ica)
  }else{
    stopifnot(prefix == get('prefix', envir=Space.ica))
  }
  data_file = get('data_files', envir=Space.ica)
  
  ### Estimate & align +/- ICA sources for pos. skew ###
  if (exists('aligning_sources') && aligning_sources){
    align.sources = F   # safety switch to prevent endlessly looping fn. call
  }
  ### Estimate sources, if indicated ###
  if (align.sources){   # Align +/-ICA sources for pos. skew
    ica = align_ICAsources_toSkew(ica, 
                                  return.data=T, 
                                  return.sources=return.sources,
                                  verbose=verbose)
  }else{
    if (multimodel.ica){
      stopifnot('X' %in% names(ica[[1]])) # cannot est. sources w/o original ICA input
      for (m in 1:length(ica)){
        if (return.sources){
          ica[[m]]$S = ica[[1]]$X[,c(1:ica[[m]]$model_order)] %*% ica[[m]]$K %*% ica[[m]]$W
        }
        ica[[m]]$A = solve(ica[[m]]$W)    # ensure correct signs in inverses
        ica[[m]]$skew_sign = NA
      }
    }else{
      stopifnot('X' %in% names(ica)) # cannot est. sources w/o original ICA input
      if (return.sources){
        ica$S = ica$X %*% ica$K %*% ica$W
      }
      ica$A = solve(ica$W)           # ensure correct signs in inverses
      ica$skew_sign = NA
    }
  }
  
  
  ### Concat results for multi-model order ICA ###
  if (flatten.multimodel.ICA){
    if (verbose){  cat('\n......flattening multi-modal ICA output...')}
    X = NULL
    K = NULL
    W = NULL
    A = NULL
    S = NULL
    KW = NULL
    model_orders = c()
    skew_signs = c()
    for (m in 1:length(ica)){
      model_orders = c(model_orders, ica[[m]]$model_order)
      
      if (flatten.source.est){  
        ### Collapse source est. steps into single, smaller matrix KW ###
        K = NA
        W = NA
        A = NA
        if (is.null(KW)){
          KW = ica[[m]]$K %*% ica[[m]]$W
        }else{
          m1 = nrow(KW)
          n1 = ncol(KW)
          m2 = nrow(ica[[m]]$K)
          n2 = ncol(ica[[m]]$W)
          KW = rbind(cbind(KW, matrix(0, m1, n2), deparse.level=0),
                     cbind(matrix(0, m2, n1), ica[[m]]$K %*% ica[[m]]$W,
                           deparse.level=0), deparse.level=0)
        }
        
      }else{
        KW = NA
        if (is.null(K)){
          K = ica[[m]]$K
          W = ica[[m]]$W
          A = ica[[m]]$A
        }else{
          n1 = nrow(K)
          n2 = ncol(K)
          m1 = nrow(ica[[m]]$K)
          m2 = ncol(ica[[m]]$K)
          K = rbind(cbind(K, matrix(0, n1, m2), deparse.level=0),
                    cbind(matrix(0, m1, n2), ica[[m]]$K, 
                          deparse.level=0), deparse.level=0)

          n1 = nrow(W)
          n2 = ncol(W)
          m1 = nrow(ica[[m]]$W)
          m2 = ncol(ica[[m]]$W)
          W = rbind(cbind(W, matrix(0, n1, m2), deparse.level=0),
                    cbind(matrix(0, m1, n2), ica[[m]]$W,
                          deparse.level=0), deparse.level=0)
          stopifnot(ncol(K) == nrow(W))
          
          n1 = nrow(A)
          n2 = ncol(A)
          m1 = nrow(ica[[m]]$A)
          m2 = ncol(ica[[m]]$A)
          A = rbind(cbind(A, matrix(0, n1, m2), deparse.level=0),
                    cbind(matrix(0, m1, n2), ica[[m]]$A, 
                          deparse.level=0), deparse.level=0)
        }
      }
      
      if ('skew_sign' %in% names(ica[[m]])){
        skew_signs = c(skew_signs, ica[[m]]$skew_sign)
      }
      if (return.sources){
        if (is.null(S)){
          S = ica[[m]]$S
        }else{
          S = cbind(S, ica[[m]]$S, deparse.level=0)
        }
      }else{  S = NA  }
      if (return.input.data){
        if ('X' %in% names(ica[[m]])){
          X = cbind(X, ica[[m]]$X, deparse.level=0)
        }
      }else{  X = NA  }
    }
    
    ica = list('X' = X,        # data input into ICA
               'K' = K,        # concatenated pre-whitening matrices, prior to ICA
               'W' = W,        # semi-diagonal ICA unmixing matrices
               'KW' = KW,      # matrix product of K and W
               'A' = A,        # semi-diagonal ICA mixing matrices
               'S' = S,        # concatenated ICA sources, in cols.
               'skew_signs'   = skew_signs,   # original skew signs of ICA sources
               'model_orders' = model_orders, # ICA model orders
               'group_pca_file' = data_file,  # path to data used for ICA
               'prefix'         = prefix)     # file prefix for experiment
  }
  
  if (multimodel.ica && !flatten.multimodel.ICA){
    model_orders = c()
    for (m in 1:length(ica)){
      model_orders = c(model_orders, ncol(ica[[m]]$W))
      if (flatten.source.est){
        ica[[m]]$KW = ica[[m]]$K %*% ica[[m]]$W
      }
    }
    ica[[1]]$model_orders = model_orders
    
    if (!('group_pca_file' %in% names(ica[[1]]))){
      ica[[1]]$group_pca_file = data_file     # path to data used for ICA
    }
    if (!('prefix' %in% names(ica[[1]]))){
      ica[[1]]$prefix = prefix                # file prefix for experiment
    }
  }else if (!multimodel.ica){  # single ICA model order
    if (flatten.source.est){  ica$KW = ica$K %*% ica$W  }else{  ica$KW = NA  }
    if (!return.sources){     ica$S = NA  }
    if (!return.input.data){  ica$X = NA  }
    ica$model_orders = ncol(ica$W)
    
    if (!('group_pca_file' %in% names(ica))){
      ica$group_pca_file = data_file          # path to data used for ICA
    }
    if (!('prefix' %in% names(ica))){
      ica$prefix = prefix                     # file prefix for experiment
    }
  }
  
  ##################################################################
  return(ica)
} ##################################################################

##################################################################
align_ICAsources_toSkew <- function(data, k.indices=NA,
                                    save_file=NA,
                                    save_var=NA,
                                    return.sources=F, 
                                    return.data=T, 
                                    verbose=F){
  ##################################################################
  # Tweaks out from Independent Components Analysis (ICA)
  #   aligning all ICA source signs to a positive 3rd moment (skew),
  #   a rough work around for ICA sign indeterminacy,
  #   that ensures +/- IC wts. have some sort of meaning in the data
  #
  # Input:
  #   data : output from fastICA() package, or filepath containing output
  #   k.indices : list of integer vectors, 
  #                 limits alignment & estimation to subset of sources
  #   return.sources : if False, skip returning sources to minimize mem. use
  # Output: above list with tweaked elements
  #   X : original data, pre-processed
  #   K : pre-whitening matrix
  #   W : est. unmixing matrix, w/ pos. skew
  #   A : mixing matrix, inverse of W
  #   S : source matrix, est. by projecting data onto unmixing matrix (if indicated)
  #   skew_sign : sign of original skew, from sources in fastICA()
  #
  
  library(MASS)
  
  if (verbose){ cat('\n...aligning +/- ICA sources to skew positively...') }
  
  ### Inputs & defaults ###
  save_file = as.character(save_file)
  save_var = as.character(save_var)
  
  subj_pca_suffix  = '_PCA1.RData'       # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData'    # filesuffix for group-level temporal PCA from group_PCA()
  group_ica_suffix = 'spatialICA.RData'  # filesuffix for group spatial ICA from group_spatialICA()
  br_suffix        = '_ICAsources.RData' # filesuffix for back-reconstructed subj. sources
  
  group_ica_wc = paste0('*', group_ica_suffix)
  
  ##################
  skew <- function(dat){
    dat = as.vector(dat)
    m = mean(dat)
    sigma = sd(dat)
    return( sum((dat - m)^3) / (sigma^3) )
  } ################
  
  singleton.list = F  # assumes multimodel order ICA, list of fastICA output
  if (is.character(data) && all(file.exists(data))){
    data_files = data
    
    aligning_sources = T  # safety to avoid endless looping fn. calls
    data = get_groupICA_dat(data_files,
                            align.sources=F,
                            flatten.multimodel.ICA=F,
                            return.sources=return.sources,
                            return.input.data=F,
                            verbose=verbose)
    if (!all(c('X', 'K', 'W') %in% names(data))){
      stopifnot(all(c('X', 'K', 'W') %in% names(data[[1]])))
    }else{
      singleton.list = T
      data = list(data)
    }
    
  }else{ # assume list formatted as FastICA() outputs, & format output similarily
    data_files = NA
    save_file = NA
    stopifnot(is.list(data))
    if (!all(c('X', 'K', 'W') %in% names(data))){
      stopifnot(all(c('X', 'K', 'W') %in% names(data[[1]])))
    }else{
      singleton.list = T
      data = list(data)
    }
  }
  
  if (!all(is.na(k.indices))){
    if (singleton.list){  k.indices = list(k.indices)  }
    stopifnot(is.list(k.indices))
    stopifnot(length(data) == length(k.indices))
  }
  
  for (d in 1:length(data)){  # iterate over levels for hierarchical ICA
    ica = data[[d]]
    
    if ('X' %in% names(ica)){
      X = ica$X
    }else if ('X' %in% names(data[[1]])){
      X = data[[1]]$X
      if ('model_order' %in% names(ica)){
        X = X[, 1:ica$model_order]  # limit group PCA subspace to specified dim. for multi-model order ICA
      }
    }
    if (ncol(X) != nrow(ica$K)){
      X = X[, 1:nrow(ica$K)]
    }
    
    if (all(is.na(k.indices))){
      k.inds = c(1:ncol(ica$W))  # number of sources estimated
    }else{
      if (is.logical(k.indices[[d]])){
        k.inds = c(1:length(k.indices[[d]]))[k.indices[[d]]]
      }else{
        k.inds = sort(as.integer(k.indices[[d]]))  # preserve ordering of sources in mixing matrix & spatial maps
      }
    }
    
    ica$S = X %*% ica$K %*% ica$W[,k.inds, drop=F]  # ICA equation, w/ data pre-whitening matrix K 
    sk = rep(1, ncol(ica$W))
    sk[k.inds] = apply(ica$S, 2, skew)
    
    if (any(sk < -.Machine$double.eps)){
      sk = sign(sk)
      
      ### Flip Signs:   W is the [~time x sources] unmixing matrix   ###
      ica$W = ica$W %*% diag(sk)       # right-multiplication flips signs of sources
      if (return.sources){
        ica$S = X %*% ica$K %*% ica$W[,k.inds, drop=F]  # re-unmix sources from data
      }else{  ica$S = NULL  }
      ica$A = ginv(ica$W)              # re-calc. mixing matrix w/ sign change(s)
      ica$skew_sign = sk               # save original direction of skew signs
    }else{
      ica$A = ginv(ica$W)              # necessary to correct signs in inverses
      if (return.sources){
        ica$S = X %*% ica$K %*% ica$W[,k.inds, drop=F]  # unmix sources from data
      }else{  ica$S = NULL  }
      ica$skew_sign = abs(sk)
    }
    ica$k.indices = k.inds            # save indices of aligned & est. sources
    
    data[[d]] = ica
  }
  if (singleton.list){  data = data[[1]]  }  # return list to original format
  
  ### Save output ###
  if (!all(is.na(save_file))){
    if (is.logical(save_file) && save_file && !all(is.na(data_files))){
      save_file = data_files
    }
    stopifnot(file.exists(save_file))
    if (all(is.na(save_var))){
      save_var = 'ica'  # default name for fastICA output
    }
    tmp = new.env()
    load(save_file, envir=tmp)
    assign(save_var, data, envir=tmp)
    save(list=ls(envir=tmp),
         file=save_file,
         envir=tmp)
  }
  
  #########################################
  if (return.data){  return(data)  }
} #########################################



##################################################################
groupICA_hierarchicalAnalysis <- function (data, 
                                           similarity.method='dcor',
                                           hPCA_algorithm=T,
                                           h.levels=NA,
                                           test.levels=F,
                                           prefix=NA, 
                                           save_path=NA,
                                           return.data=F, 
                                           parallel=F,
                                           verbose=F){
  ##################################################################
  # Perform hierarchical analysis, on group-level spatial ICA data,
  #   Either hierarchical PCA ("treelets PCA"), or hierarchical ICA.
  #
  # Inputs:
  #   data : (options)
  #           1. file w/ group ICA space, output from group_spatialICA()
  #           2. path to directory containing above
  #           3. [voxels x time subspace] matrix of reduced dim. group data
  #   similarity.method  : similarity measure type {'dcor', 'bcdcor', 'cor', etc.}
  #   hPCA_algorithm : if T apply hierarchical PCA algorithm ("treelet PCA", recommended),
  #                       else apply hierarchical ICA algorithm
  #   h.levels    : Number of levels of the hierarchy to estimate (from bottom).
  #                   If NA, defaults to number of sources or dim. of data.
  #                   If 'test' or 'stats', levels determined by stats. testing.
  #                   If 'resume' or 'continue', additional levels added to prev. hierarchy 
  #   test.levels : Statistically test similarities at each level of the hierarchy,
  #                   applying Eigenvalue Ratio Test
  #                   to identify potential mergers of duplicate vars. at low levels,
  #                   & appropriate test for sim. method at higher levels
  #   prefix        : identifying prefix to attach to saved output
  #   save_path     : path to dir. for saved files
  # Output: saved file w/ *_hierarchy.RData added as suffix,
  #   formatted as list w/ elements for each level:
  #     $X : original data, pre-processed
  #     $K : pre-ICA whitening matrix        [~time x reduced dim.]
  #     $W : est. unmixing matrix            [reduced dim. x sources]
  #     $A : mixing matrix, inverse of W     [sources x reduced dim.]
  # Requires:
  #   get_groupICA_dat()        : load & format ICA info or multi-order ICA
  #   similarity.measure()      : calculates similarity & finds dependencies
  #   similarity_largeData()    : re-implementation of similarity_hica() from fastHICA
  #   hierarchicalICA_basis()   : re-implementation of basis_hica() from fastHICA
  #   hierarchicalPCA_basis()   : for treelet PCA levels
  #   eigenValue_ratio_test()   : tests smallest eigenvalue vs. sum of all eigenvals
  #   fastICA()                 : fast ICA algorithm from fastICA package
  #   find_initial_n.s._levels() : finds initial sequence of n.s. levels
  #  
  
  if (verbose){  cat('\nEstimating hierarchical ICA / PCA...')}
  library(fastICA)
  library(MASS)
  
  ### Inputs & defaults ###
  save.all.levels = F     # if False, skip saving data from initial, redundant non-sig. levels
  similarity.current = NA
  
  similarity.method = as.character(similarity.method)
  hPCA_algorithm = as.logical(hPCA_algorithm)
  test.levels = as.logical(test.levels)
  if (!all(is.na(h.levels)) && is.character(h.levels)){
    stopifnot(grepl('stat', h.levels) || grepl('test', h.levels) ||     # keywords to set levels based on stats. testing
              grepl('continue', h.levels) || grepl('resume', h.levels)) # keywords to resume constructing saved hierarchy
  }else{
    h.levels = as.integer(h.levels)
  }
  prefix = as.character(prefix)
  save_path = as.character(save_path)
  
  subj_pca_suffix  = '_PCA1.RData'       # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData'    # filesuffix for group-level temporal PCA from group_PCA()
  group_ica_suffix = 'spatialICA.RData'  # filesuffix for group spatial ICA from group_spatialICA()
  hier_suffix = 'hierarchy.RData'   # filesuffix for hierarchical ICA/PCA from groupICA_hierarchicalAnalysis()
  br_suffix   = '_ICAsources.RData' # filesuffix for back-reconstructed subj. sources  
  
  group_ica_wc = group_ica_suffix
  if (!all(is.na(prefix))){  group_ica_wc = paste0(prefix, '*_', group_ica_wc)  }
  hier_wc = hier_suffix
  if (!all(is.na(prefix))){  hier_wc = paste0(prefix, '*_', hier_wc)  }
  
  ### Load data ###
  data_files = NA
  if (is.character(data) && all(file.exists(data))){
    data_files = data
    
    if (is.na(save_path)){
      if ((length(data_files)==1) && dir.exists(data_files)){
        save_path = data_files
      }else{
        save_path = unique(dirname(data_files))[1]
      }
    }
    
    ### Load partially complete hierarchy & resume ###
    if (grepl('continue', h.levels) || grepl('resume', h.levels)){
      if ((length(data_files) == 1) && dir.exists(data_files) &&   
          (length(list.files(data_files, hier_wc)) == 1)){
        hier_file = list.files(data_files, hier_wc, full.names=T)[1]
      }
      
      hier = get_hier_dat(hier_file, prefix=prefix, 
                          align.sources=F,
                          return.input=T, 
                          return.sources=F, 
                          return.similarity=T,
                          verbose=verbose)
      hier_ns = get_hier_n.s._dat(hier_file, prefix=prefix,
                                  discard.initial.levels=T,
                                  align.sources=F,
                                  return.input=T, 
                                  return.sources=F,  
                                  return.similarity=T,
                                  verbose=verbose)
      
      if (length(hier) > 0){
        h.max = hier[[length(hier)]]$h.level
      }else{  h.max = -Inf  }
      if (length(hier_ns) > 0){
        h_ns.max = hier_ns[[length(hier_ns)]]$h.level
      }else{  h_ns.max = -Inf  }
      stopifnot(is.finite(h.max) || is.finite(h_ns.max))
      if (h.max > h_ns.max){
        stopifnot(all(c('X', 'K', 'W', 'group_ica_file',  'similarity.current')
                  %in% names(hier[[1]])))
        l = length(hier)
        data = hier[[1]]$X %*% hier[[l]]$K %*% hier[[l]]$W   # reconstruct data from last level of hierarchy
        data_files = hier[[1]]$group_ica_file
        similarity.current = hier[[1]]$similarity.current    # last state of similarity matrix
        hier[[1]][c('group_ica_file',  'similarity.current', 'h.levels.criteria', 'prefix')] = NULL
      }else if (h.max < h_ns.max){
        stopifnot(all(c('X', 'K', 'W', 'group_ica_file',  'similarity.current') 
                  %in% names(hier_ns[[1]])))
        l = length(hier_ns)
        data = hier_ns[[1]]$X %*% hier_ns[[l]]$K %*% hier_ns[[l]]$W  # reconstruct data from last level of hierarchy
        data_files = hier_ns[[1]]$group_ica_file
        similarity.current = hier_ns[[1]]$similarity.current         # last state of similarity matrix
        hier_ns[[1]][c('group_ica_file',  'similarity.current', 'h.levels.criteria', 'prefix')] = NULL
      }else{
        stopifnot(any(is.finite(h.max), is.finite(h_ns.max)))
      }
      h.levels = NA  # algorithm starts where above left off, appends to results
      
    }else{
      ### Load non-hierarchical ICA data ###
      if ((length(data_files) == 1) && dir.exists(data_files) &&   
          (length(list.files(data_files, group_ica_wc)) == 1)){
        group_ica_file = list.files(data_files, group_ica_wc, full.names=T)[1]
      }
      data = get_groupICA_dat(group_ica_file, prefix=prefix, 
                              align.sources=F,
                              flatten.multimodel.ICA=T,
                              flatten.source.est=T,
                              return.sources=T,
                              return.input.data=F,
                              verbose=verbose)
      data_files = group_ica_file
      
      if ('model_orders' %in% names(data)){
        pre_hier_model_orders = data$model_orders
      }else if ('model_orders' %in% names(data[[1]])){
        pre_hier_model_orders = data[[1]]$model_orders
      }
      if ('KW' %in% names(data)){
        pre_hier_KW = data$KW  # matrix product of ICA pre-whitening & unmixing matrices, for sources used for hPCA/hICA
      }else if ('KW' %in% names(data[[1]])){
        pre_hier_KW = data[[1]]$KW
      }else{  pre_hier_KW = NA  }
      if ('S' %in% names(data)){
        data = data$S   # apply hPCA/hICA to all ICA sources, including multi-modal order ICA sources
      }
    }
  }
  
  
  ### Hierarchy params. & data dims. ###
  L.max = dim(data)[2] - 1    # max. levels is the number of ICA sources - 1
  if (is.na(h.levels) || (is.numeric(h.levels) && ((h.levels > L.max)  || (h.levels <= 0)))){
    h.levels = L.max
  }
  V = dim(data)[1]            # number of voxels
  
  ### Stat testing params. ###
  if (test.levels){ # exclude levels that merge duplicates
    delta.eigs = 0.2      # ratio of 2nd eigenvalue to sum of eigs.
    pvalue.eigs = 0.001   # p-value for above, MCP merged into below
    pvalue.sim  = 0.001   # p-value for merging levels, uncorrected
    if (is.numeric(h.levels)){  # p-value for merging levels, FWE-corrected
      pvalue.sim = pvalue.sim / h.levels 
    }else{
      pvalue.sim = pvalue.sim / L.max
    }
  }else{            # skip stats
    delta.eigs = NA
    pvalue.eigs = NA
    pvalue.sim = NA
  }
  
  
  ### Main fn. ###
  if (hPCA_algorithm){
    if (verbose){
      cat('\n...calculating hierarchical ("treelets") PCA algorithm of Lee et al. (2008)...\n')
    }
    hier_algorithm = 'hPCA'
    
    basis_L = hierarchicalPCA_basis(data, maxlev=h.levels,
                                    sim.method=similarity.method,
                                    sim.prev=similarity.current,
                                    delta=delta.eigs, pval=pvalue.eigs,
                                    parallel=parallel, verbose=verbose)
    
  }else{
    if (verbose){
      cat('\n...calculating hierarchical ICA using the algorithm of Secchi et al. (2016)...\n')
    }
    hier_algorithm = 'hICA'
    
    basis_L = hierarchicalICA_basis(data, maxlev=h.levels,
                                    sim.method=similarity.method,
                                    sim.prev=similarity.current,
                                    delta=delta.eigs, pval=pvalue.eigs,
                                    parallel=parallel, verbose=verbose)
  }
  stopifnot(length(basis_L$basis) > 0)
  
  ### Sort, assemble & format levels ###
  if (verbose){  cat('\n\nSummary of each level of hierarchy:')}
  if (!exists('hier') || all(is.na(hier))){  hier = list()  }
  if (!exists('hier_ns') || all(is.na(hier_ns))){  hier_ns = list()  }
  hl = 0       # list index for main formatted hierarchy
  hl.ns = 0    # list index for formatted n.s. levels
  hl.prev = 0  # levels of hierarchy from previous application of fn.
  if ((length(hier) > 0) && ('h.levels.max' %in% names(hier[[1]]))){
    hl.prev = max(hl.prev, hier[[1]]$h.levels.max, na.rm=T)
  }
  if ((length(hier_ns) > 0) && ('h.levels.max' %in% names(hier_ns[[1]]))){
    hl.prev = max(hl.prev, hier_ns[[1]]$h.levels.max, na.rm=T)
  }
  if (any(is.na(basis_L$similarity$aggregation[,1]))){
    ### Find & remove levels w/ manually excluded vars. from aggregation table ###
    l.k.rm = which(is.na(basis_L$similarity$aggregation[,1]))
    basis_L$similarity$aggregation = basis_L$similarity$aggregation[-l.k.rm,, drop=F]
  }
  
  for (l in 1:length(basis_L$basis)){
    if (test.levels){
      skip.level = F
      if (basis_L$eigenvalues[[l]]$lambda_ratio.test){  
        skip.level = T        # skip levels merging redundant sources
        test = basis_L$eigenvalues[[l]]  # format output ~ t.test()
        test$statistic = test$lambda_ratio.statistic
        names(test$statistic) = 'Lambda Ratio Test'
        test$p.value = test$lambda_ratio.pvalue
        test$parameter = test$lambda_ratio.delta
        names(test$parameter) = 'delta'
        test$similarity = basis_L$eigenvalues[[l]]$similarity

      }else{
        if (similarity.method %in% c('calc_corrs', 'cor', 'pearson')){
          r = basis_L$eigenvalues[[l]]$similarity
          t = r * sqrt(V-2) / sqrt(1 - r^2)  # Anderson. (2003). Wiley, New Jersey. p. 121. 
          ttest = list('statistic' = t,
                        'parameter' = V-2,
                        'p.value'   = dt(t, V-2),
                        'estimate'  = r,
                        'similarity'= r)
          names(ttest$statistic) = 't'
          names(ttest$parameter) = 'df'
          names(ttest$estimate)  = 'cor'
          
        }else if (similarity.method == 'bcdor'){
          cat('\nWarning:  insufficient observations (n=2 comps.) to test bias-corrected distance corr.!')
          ttest = list('statistic' = NA, 'p.value' = -Inf)
        }else{
          cat(paste0('\nWarning:  statistical inference for  ',similarity.method,'  not implemented!'))
          ttest = list('statistic' = NA, 'p.value' = -Inf)
        }
        test = ttest
        if (ttest$p.value >= pvalue.sim){
          skip.level = T    # skip levels merging fully independent or uncorrelated sources
        }
      }
      if (skip.level){
        ### Non-sig. level ###
        if (verbose){
          testname = names(test$statistic)
          if (testname == 't'){  testname = paste0('T(',test$parameter,')')  }
          sim = signif(test$similarity,3)
          pval = signif(test$p.value,3)
          if (testname == 'Lambda Ratio Test'){
            cat(paste0('\nlevel  ',l+hl.prev,':   similarity=',sim,'  indicates redundant sources   by ',testname,', p < ', pval))
          }else{
            cat(paste0('\nlevel  ',l+hl.prev,':   similarity=',sim,'  n.s. by ',testname,', p > ', pval))
          }
        }
        
        hl.ns = length(hier_ns) + 1
        aggregation_level = unname(basis_L$similarity$aggregation[l+hl.prev,])

        hier_ns[[hl.ns]] = list('h.level' = l + hl.prev,               # saved hierarchy level
                                'K' = basis_L$prewhitening[[l]],       # pre-whitening, if applied prior to ICA/PCA
                                'A' = t(basis_L$basis[[l]]),           # mixing matrix (transposed basis matrix)
                                'W' = t(ginv(basis_L$basis[[l]])),     # unmixing matrix, inverse of level's basis matrix
                                'k1'         = aggregation_level[1],   # index of new comp, incorporated into subsequent levels of hieararchy 
                                'k2'         = aggregation_level[2],   # index of new 2nd comp, stored & constant therafter
                                'similarity' = aggregation_level[3],   # sim. between comps at prev. level, entered into level's ICA/PCA
                                'variances'  = aggregation_level[4:5], # variances of new comps
                                'active_inds'= basis_L$active_inds[[l]], # indices of active vars., incorporated in higher levels
                                'algorithm'  = hier_algorithm)         # hPCA or hICA algorithm  used to construct level
        if ('similarity' %in% names(test)){  test$similarity = NULL  }
        hier_ns[[hl.ns]] = append(hier_ns[[hl.ns]], test)
        next
        
      }else{
        if (verbose && !all(is.na(test$statistic))){
          testname = names(test$statistic)
          if (testname == 't'){  testname = paste0('T(',test$parameter,')')  }
          est = signif(test$estimate,3)
          pval = signif(test$p.value,3)
          cat(paste0('\nlevel  ',l+hl.prev,':   similarity=',est,'         *sig. by ',testname,', p = ', pval))
        }
      }
    }
    
    ### Assemble hierarchy ###
    hl = length(hier) + 1
    aggregation_level = basis_L$similarity$aggregation[l+hl.prev,]
    
    hier[[hl]] = list('h.level' = l + hl.prev,               # saved hierarchy level
                      'K' = basis_L$prewhitening[[l]],       # pre-whitening, if applied prior to ICA/PCA
                      'A' = t(basis_L$basis[[l]]),           # mixing matrix (transposed basis matrix)
                      'W' = t(ginv(basis_L$basis[[l]])),     # unmixing matrix, inverse of level's basis matrix
                      'k1'         = aggregation_level[1],   # index of new comp, incorporated into subsequent levels of hieararchy 
                      'k2'         = aggregation_level[2],   # index of new 2nd comp, stored & constant therafter
                      'similarity' = aggregation_level[3],   # sim. between comps at prev. level, entered into level's ICA/PCA
                      'variances'  = aggregation_level[4:5], # variances of new comps
                      'active_inds'= basis_L$active_inds[[l]], # indices of active vars., incorporated in higher levels
                      'algorithm'  = hier_algorithm)         # hPCA or hICA algorithm used to construct level
    if (length(aggregation_level) > 5){
      hier[[hl]]$sim_w_PC1          = aggregation_level[6:7] # similarities of new comps w/ leading PC1
    }
  }
  
  ### Remove irrelevant data from levels before saving ###
  if (!save.all.levels && (length(hier_ns) > 0)){
    end.hl.seq = find_initial_n.s._levels(hier_ns)
    
    if (end.hl.seq > 1){
      for (l in 1:(end.hl.seq-1)){
        hier_ns[[l]]$K = NULL
        hier_ns[[l]]$A = NULL
        hier_ns[[l]]$W = NULL
      }
    }
  }
  if ((length(hier) > 0) && ('X' %in% names(hier[[1]]))){
    data = hier[[1]]$X   # save original data, if stopped & restarted
    hier[[1]]$X = NULL
  }
  if ((length(hier_ns) > 0) && ('X' %in% names(hier_ns[[1]]))){
    data = hier_ns[[1]]$X   # save original data, if stopped & restarted
    hier_ns[[1]]$X = NULL
  }
  if ((length(hier) > 0) && ('similarity.current' %in% names(hier[[1]]))){
    hier[[1]]$similarity.current = NULL
  }
  if ((length(hier_ns) > 0) && ('similarity.current' %in% names(hier_ns[[1]]))){
    hier_ns[[1]]$similarity.current = NULL
  }
  
  ### Get info to construct additional levels ###
  h.levels.criteria  = h.levels             # input used to determine number of levels of hierarchy
  h.levels.max.hier = 0                     # highest level in main hierarchy
  if (length(hier) > 0){
    h.levels.max.hier = hier[[length(hier)]]$h.level
  }
  h.levels.max.hier_ns = 0                  # highest non-sig. level in hierarchy
  if (length(hier_ns) > 0){
    h.levels.max.hier_ns  = hier_ns[[length(hier_ns)]]$h.level
  }
  h.levels.max = max(h.levels.max.hier, h.levels.max.hier_ns)  # highest current level in hierarchy
  h.levels.new       = length(basis_L$basis)    # number of levels from recent run
  similarity.current = basis_L$similarity   # sim. matrix as list, used to construct additional hierarchy levels
  
  if (!exists('pre_hier_KW') || all(is.na(pre_hier_KW))){  # product of ICA pre-whitening matrix K & unmixing matrix W...
    if ('KW' %in% names(data)){ # ...applied to prior to hPCA/hICA, used in back-reconstruction
      pre_hier_KW = data$KW   
    }else if ('KW' %in% names(data[[1]])){
      pre_hier_KW = data[[1]]$KW
    }else if ((length(hier) > 0) && ('pre_hier_KW' %in% names(hier[[1]]))){
      pre_hier_KW = hier[[1]]$pre_hier_KW
    }else if ((length(hier_ns) > 0) && ('pre_hier_KW' %in% names(hier_ns[[1]]))){
      pre_hier_KW = hier_ns[[1]]$pre_hier_KW
    }else{
      pre_hier_KW = NA
    }
  }
  if (!exists('pre_hier_model_orders') || all(is.na(pre_hier_model_orders))){ # pre-hPCA/hICA source ICA model orders
    if ('model_orders' %in% names(data)){
      pre_hier_model_orders = data$model_orders
    }else if (is.list(data) && (length(data) > 0) && 'model_orders' %in% names(data[[1]])){
      pre_hier_model_orders = data[[1]]$model_orders
    }else if ((length(hier) > 0) && ('pre_hier_model_orders' %in% names(hier[[1]]))){
      pre_hier_model_orders = hier[[1]]$pre_hier_model_orders
    }else if ((length(hier_ns) > 0) && ('pre_hier_model_orders' %in% names(hier_ns[[1]]))){
      pre_hier_model_orders = hier_ns[[1]]$pre_hier_model_orders
    }else{
      pre_hier_model_orders = NA
    }
  }
  
  
  ### Saving ###
  save.fname = hier_suffix
  if (!is.na(prefix)){
    save.fname = paste0(prefix, '_', save.fname)
  }else{
    save.fhame = sub('^_', '', save.fname)
  }
  if (!is.na(save_path)){
    save.fname = file.path(save_path, save.fname)
  }
  if (verbose){  cat(paste0('\n\n......& saving as: ', save.fname, '\n\n'))}
  save(list=c('hier', 'hier_ns', 'h.levels.criteria', 'similarity.current',
              'similarity.method', 'hier_algorithm',
              'h.levels.max', 'h.levels.new', 'h.levels.max.hier', 'h.levels.max.hier_ns',
              'pre_hier_KW', 'pre_hier_model_orders',
              'data', 'data_files', 'prefix'), 
       file=save.fname)
  stopifnot(file.exists(save.fname))
  
  ##################################################################
  if (return.data){  return(hier)  }
} ##################################################################


##################################################################
get_hier_dat <- function(hier_file, prefix=NA,
                         align.sources=F,
                         hier_levels=NA,
                         return.input=T,
                         return.sources=T,
                         return.similarity=F,
                         verbose=F){
  ##################################################################
  # Loads group-level hierarchical ICA/PCA,
  #   output by groupICA_hierarchicalAnalysis(),
  #     prior to subj.-level ICA back-reconstruction
  #
  # Requires:  align_ICAsources_toSkew()
  #
  
  if (verbose){cat(paste0('\n...loading hierarchical analysis data from:\n      ', hier_file))}
  
  stopifnot(file.exists(hier_file))
  
  Space.hier = new.env()
  load(hier_file, envir=Space.hier)
  stopifnot(all(c('hier', 'h.levels.criteria', 'data_files', 'prefix') %in% names(Space.hier)))
  hier = get('hier', envir=Space.hier)
  if (length(hier) == 0){  return(hier)  }

  stopifnot(all(c('K','W','A','k1','k2') %in% names(hier[[1]])))
  if (all(is.na(prefix))){
    prefix = get('prefix', envir=Space.hier)
  }else{
    stopifnot(prefix == get('prefix', envir=Space.hier))
  }
  
  if (!all(is.na(hier_levels))){
    h.levels.saved = unlist(lapply(hier, function(yy) return(yy$h.level)))
    
    if (is.character(hier_levels) && ((hier_levels == 'final') || (hier_levels == 'last'))){
      hier_levels = length(hier)  # return all sources at highest/final level of hierarchy
    }else if (is.character(hier_levels) && (hier_levels == 'all')){
      hier_levels = 1:length(hier) # return all sources at each/every level of heirarchy
    }else if (is.numeric(hier_levels)){
      hier_levels = which(h.levels.saved %in% hier_levels)  # return sources at select levels
    }
    hier = hier[hier_levels]  # discard other levels
  }else{
    hier_levels = F
  }
  
  ### Add meta- info to 1st level ###
  hier[[1]]$h.levels.criteria = get('h.levels.criteria', envir=Space.hier)
  hier[[1]]$group_ica_file = get('data_files', envir=Space.hier)
  hier[[1]]$prefix = prefix
  hier[[1]]$h.levels.max = get('h.levels.max', envir=Space.hier)
  hier[[1]]$pre_hier_KW = get('pre_hier_KW', envir=Space.hier)
  hier[[1]]$pre_hier_model_orders = get('pre_hier_model_orders', envir=Space.hier)
  if ('k.indices.rm' %in% names(Space.hier)){
    hier[[1]]$k.indices.rm = get('k.indices.rm', envir=Space.hier)
  }
  
  if (return.input || return.sources){
    hier[[1]]$X = get('data', envir=Space.hier)  # get input to fn.
  }
  
  ### Estimate sources, if indicated ###
  if (return.sources){
    k.indices = list()
    for (l in 1:length(hier)){
      if (any(as.logical(hier_levels))){
        k.indices[[l]] = c(1:ncol(hier[[l]]$W))
      }else{
        k.indices[[l]] = sort(c(hier[[l]]$k1, hier[[l]]$k2)) # preserve original ordering of sources in mixing matrix & spatial maps
      }
    }
    
    if (align.sources){
      hier = align_ICAsources_toSkew(hier, 
                                     k.indices=k.indices,
                                     return.sources=return.sources,
                                     return.data=T,
                                     verbose=verbose)
    }else{
      for (l in 1:length(hier)){
        if (!('S' %in% names(hier[[l]])) || all(is.na(hier[[l]]$S))){
          stopifnot('X' %in% names(hier[[1]])) # cannot est. sources w/o original ICA input
          hier[[l]]$S = hier[[1]]$X %*% hier[[l]]$K %*% hier[[l]]$W[,k.indices[[l]], drop=F]
        }
      }
    }
  }
  
  ### Return similarity matrix, if indicated ###
  if (return.similarity){
    hier[[1]]$similarity.current = get('similarity.current', envir=Space.hier)
  }
  
  ##################################################################
  return(hier)      # list w/ hierarchical ICA/PCA sources & basis'
} ##################################################################

##################################################################
get_hier_n.s._dat <- function(hier_file, prefix=NA,
                              discard.initial.levels=T,
                              align.sources=F,
                              hier_levels=NA,
                              return.input=T,
                              return.sources=T,
                              return.similarity=F,
                              verbose=F){
  ##################################################################
  # Loads non-significant group-level hierarchical ICA/PCA sources,
  #   output by groupICA_hierarchicalAnalysis(),
  #     prior to subj.-level ICA back-reconstruction.
  # If "discard.initial.levels", discards all but last in initial n.s. levels,
  #   composed of mergers of redundant sources,
  #     returning last in level in initial n.s. sequence & subsequent n.s. levels
  #
  # Requires:    find_initial_n.s._levels()
  #              align_ICAsources_toSkew()
  #
  
  if (verbose){cat(paste0('\n...loading initial n.s. hierarchical data from:\n      ', hier_file))}
  
  stopifnot(file.exists(hier_file))
  
  Space.hier = new.env()
  load(hier_file, envir=Space.hier)
  stopifnot(all(c('hier_ns', 'h.levels.criteria', 'data_files', 'prefix') %in% names(Space.hier)))
  hier_ns = get('hier_ns', envir=Space.hier)
  if (length(hier_ns) == 0){  return(hier_ns)  }
  
  stopifnot(all(c('k1','k2','active_inds') %in% names(get('hier_ns', envir=Space.hier)[[1]])))
  if (all(is.na(prefix))){
    prefix = get('prefix', envir=Space.hier)
  }else{
    stopifnot(prefix == get('prefix', envir=Space.hier))
  }
  
  if (discard.initial.levels){ # discard initial, redundant levels
    end.hl.seq = find_initial_n.s._levels(hier_ns)
    hL_ns = length(hier_ns)
    if (all(is.na(hier_levels))){
      hier_ns = hier_ns[end.hl.seq:hL_ns]
    }
  }
  
  if (!all(is.na(hier_levels))){
    h.levels.saved = unlist(lapply(hier_ns, function(yy) return(yy$h.level)))
    
    if (is.character(hier_levels) && ((hier_levels == 'final') || (hier_levels == 'last'))){
      hier_levels = length(hier_ns)  # return all sources at highest/final level of hierarchy
    }else if (is.character(hier_levels) && (hier_levels == 'all')){
      hier_levels = 1:length(hier_ns) # return all sources at each/every level of heirarchy
    }else if (is.numeric(hier_levels)){
      hier_levels = which(h.levels.saved %in% hier_levels)  # return sources at select levels
    }
    if (discard.initial.levels){   # ensure final n.s. level is kept
      hier_levels = unique(c(end.hl.seq, hier_levels))
    }
    hier_ns = hier_ns[hier_levels]  # discard other levels
  }else{
    hier_levels = F
  }
  

  ### Add meta- info to 1st level ###
  hier_ns[[1]]$h.levels.criteria = get('h.levels.criteria', envir=Space.hier)
  hier_ns[[1]]$group_ica_file = get('data_files', envir=Space.hier)
  hier_ns[[1]]$prefix = prefix
  hier_ns[[1]]$h.levels.max = get('h.levels.max', envir=Space.hier)
  hier_ns[[1]]$pre_hier_KW = get('pre_hier_KW', envir=Space.hier)
  hier_ns[[1]]$pre_hier_model_orders = get('pre_hier_model_orders', envir=Space.hier)
  if ('k.indices.rm' %in% names(Space.hier)){
    hier_ns[[1]]$k.indices.rm = get('k.indices.rm', envir=Space.hier)
  }
  
  if (return.input || return.sources){
    hier_ns[[1]]$X = get('data', envir=Space.hier)  # get input to fn.
  }
  
  ### Estimate sources, if indicated ###
  if (return.sources){
    k.indices = list()
    for (l in 1:length(hier_ns)){
      if (any(as.logical(hier_levels))){
        k.indices[[l]] = c(1:ncol(hier_ns[[l]]$W))
      }else{
        k.indices[[l]] = sort(hier_ns[[l]]$active_inds)  # preserve original ordering of sources in mixing matrix & spatial maps
      }
    }
    
    if (align.sources){
      hier_ns = align_ICAsources_toSkew(hier_ns,
                                        k.indices=k.indices,
                                        return.sources=return.sources,
                                        return.data=T,
                                        verbose=verbose)
    }else{
      for (l in 1:length(hier_ns)){
        if (!('S' %in% names(hier_ns[[l]])) || all(is.na(hier_ns[[l]]$S))){
          stopifnot('X' %in% names(hier_ns[[1]])) # cannot est. sources w/o original ICA input
          hier_ns[[l]]$S = hier_ns[[1]]$X %*% hier_ns[[l]]$K %*% hier_ns[[l]]$W[,k.indices[[l]], drop=F]
        }
      }
    }
  }
  ### Return similarity matrix, if indicated ###
  if (return.similarity){
    hier_ns[[1]]$similarity.current = get('similarity.current', envir=Space.hier)
  }
  
  ##################################################################
  return(hier_ns)   # list w/ initial hierarchical ICA/PCA sources & basis'
} ##################################################################

##################################################################
find_initial_n.s._levels <- function(hier_ns){
  ##################################################################
  # Finds last level in initial sequence of non-sig. hierarchy levels,
  #   to skip saving/returning n.s. & irrelevant info from early levels
  #
  
  stopifnot(is.list(hier_ns))
  if (length(hier_ns) == 1){  return( 1 )}
  stopifnot('h.level' %in% names(hier_ns[[1]]))
  
  ns.h.levels = unlist(lapply(hier_ns, function(yy) return(yy$h.level)))
  if (length(ns.h.levels) == 0){
    end.hl.seq = NA
    return(  hier_ns  )  # nothing else to do
  }else if (length(ns.h.levels) == 1){
    end.hl.seq = 1
  }else{ # find end of sequence of initial n.s. levels
    end.hl.seq = length(ns.h.levels)
    for (l in 2:length(ns.h.levels)){
      if ((ns.h.levels[l] - ns.h.levels[l-1]) == 1){
        next
      }else{
        end.hl.seq = l - 1
        break
      }
    }
  }
  ##################################################################
  return(end.hl.seq)  # index of last level in n.s. sequence
} ##################################################################

##################################################################
eigenValue_ratio_test <- function(X, d=0.1, m=1, pval=0.001, 
                                  return.pca=T, return.PCs=F, 
                                  verbose=F){
  ##################################################################
  # Eigenvalue ratio test,
  #   to test for sig. of sum of smallest m+1 eigenvalues:
  #
  # H_0: (l_m+1 + l_m+2 + ... l_p) / (l_1 + ... + l_p) > d
  #      sum of smallest eigenvalues / all eigenvalues > d
  # H_A:      l_numerator / l_denominator              < d
  #   first m principal components represent all measurements
  #
  #   Anderson, T.W. (2003). An introduction to multivariate 
  #     satistical analysis (3rd ed). Wiley, New Jersey. p. 480.
  #
  
  X = as.matrix(X)
  d = as.numeric(d)
  m = as.integer(m)
  pval = as.numeric(pval)
  
  pca = prcomp(X, 
               retx=return.PCs, center=T, scale.=T)
  n = nrow(X) - 1
  p = length(pca$sdev)
  l_num = pca$sdev[(m+1):p]^2  # smallest eigenvalues
  l_denom = pca$sdev^2
  trS = sum(l_denom)  # trace of sample cov. matrix
  
  l_ratio = sum(l_num) / sum(l_denom)
  cr = 2 * (d)^2 + 2 * ((1 - d) / (trS))^2 * sum(l_num^2) # ref. eq. (7)
  z = qnorm(1 - pval)            # upper critical point of N(0,1), w/ sig. level pval
  critical_value = -z * sqrt(cr) # lower critical point, scaled to asymptotic s.d.
  stat = sqrt(n) * (l_ratio - d) # test statistic as observed vs. hypothesized ratio
  test = stat < critical_value   # rejection region is less than lower critical point
  
  CI.upper = l_ratio + 
    z * sqrt(2*sum(l_num)^2 * sum(l_denom^2) + 2*sum(l_denom)^2 * sum(l_num^2)) / 
    (sqrt(n) * sum(l_denom)^2)
  
  if (verbose){
    cat(paste0('\nEigenvalue Ratio Test:'))
    cat(paste0('\n  Sum of smallest eigenvalues = ', signif(sum(l_num))))
    cat(paste0('\n  Sum of all eigenvalues = ', signif(sum(l_denom))))
    cat(paste0('\n  Ratio of above sums = ', signif(l_ratio)))
    cat(paste0('\n  C.I. upper bound for ratio of eigenvalue sums = ', 
               signif(CI.upper)))
    cat(paste0('\n  Confidence level for upper bound = ', 1 - pval))
    cat(paste0('\n  Test statistic: sqrt(n)*(lambda_ratio - d) = ', signif(stat)))
    cat(paste0('\n  Critical value: test stat.  < ', signif(critical_value)))
    cat(paste0('\n  Null hypothesis rejected at p<',pval,': ', test))
  }
  
  pca$lambda = l_denom        # eigenvalues of covariance matrix
  pca$lambda_ratio = l_ratio  # ratio of sum of smallest eigenvalues to sum of all eigenvalues
  pca$lambda_ratio.delta = d  # input cutoff for ratio
  pca$lambda_ratio.m     = m  # number of largest eigenvalues to keep, if H_0 is rejected
  pca$lambda_ratio.pvalue = pval         # sig. level of test
  pca$lambda_ratio.statistic = stat      # test statistic
  pca$lambda_ratio.test   = test         # True if null hypothesis is rejected
  pca$lambda_ratio.CI.upper = CI.upper   # upper bound of confidence interval for ratio
  pca$lambda_radio.confidence = 1 - pval # confidence level for above upper bound
  
  ##################################################################
  if (return.pca){  return(pca)  }else{  return(test)  }
} ##################################################################


#################################################################
calc_corrs <- function(X, Y=NA, verbose=F){
  #################################################################
  # Wrapper fn. for C++ code to calculate correlations between large number of variables:
  #   corr_mat = corr_C(X,Y)
  #     where X,Y are [voxels by time]
  # Requires: corr_C(), cppFunction below
  #
  # NOTE ON USAGE: corr_C must be compiled w/n foreach loop,
  #     individually on each node, in order to be available during parallization
  #
  
  if (verbose){  cat('\nCalculating correlations...')  }
  stopifnot(exists('corr_C'))

  X = as.matrix(X)
  if (all(is.na(Y))){
    scale_Y = F
    Y = X
  }else{
    scale_Y = T
    Y = as.matrix(Y)
  }
  if (dim(X)[2] == 1){ X = t(X) } #use row-vectors for time series, for seed-based FC
  if (dim(Y)[2] == 1){ Y = t(Y) }
  stopifnot(dim(X)[2] == dim(t(Y))[1])
  
  ##################################################################
  return(corr_C(X,Y,scale_Y))
} ##################################################################

### Required code for fn., needs to be compiled ###
cpp_fn_string = 
            'NumericMatrix corr_C(NumericMatrix Xr, NumericMatrix Yr, bool scale_center_y = true) {
            int m = Xr.nrow(),
            j = Xr.ncol();
            int n = Yr.nrow(),
            k = Yr.ncol();
            
            arma::mat X = arma::mat(Xr.begin(), m, j, false);
            arma::mat Y = arma::mat(Yr.begin(), n, k, false);
            
            if (size(X) == size(Y)){
            if (approx_equal(X, Y, "absdiff", std::numeric_limits<double>::epsilon())){
            scale_center_y = false;
            }
            }
            X.each_col() -= mean(X,1);
            arma::colvec Xn = sqrt(sum(square(X),1));
            X.each_col() /= Xn;
            
            if (scale_center_y){
            Y.each_col() -= mean(Y,1);
            arma::colvec Yn = sqrt(sum(square(Y),1));
            Y.each_col() /= Yn;
            }else{
            Y = X;
            }
            
            arma::mat R = X * Y.t();
            R.replace(arma::datum::nan, 0);
            return wrap(R);
            
            } '
library(Rcpp)
cppFunction(depends='RcppArmadillo', cpp_fn_string)
#################################################################



##################################################################
similarity.measure <- function(x1=NA, x2=NA, method=F, 
                               return.exports=F, compile=T){
  ##################################################################
  # Meta-fn. to measure similarity using variety of methods, options:
  #   return.exports: finds hidden dependent fn. calls for parallelization,
  #   compile:  compiles C++ code, if needed
  
  library(energy)
  library(Rcpp)
  
  stopifnot((!all(is.na(x1)) && !all(is.na(x2))) || 
              (is.character(method) && return.exports))
  
  x1 = as.matrix(x1)
  x2 = as.matrix(x2)
  
  if ((method == 'dcor') && (max(nrow(x1),nrow(x2),na.rm=T) < 512)){
    if (return.exports){  return('dcor')  }
    stat <- dcor(x1, x2)  # standard dist. corr. from energy package
  }else if (method == 'dcor'){
    if (return.exports){  return('dcor2d')  }
    stat <- dcor2d(x1, x2, type='V')  # fast dist. corr., a V-statistic
  }else if ((method == 'bcdcor') && (max(nrow(x1),nrow(x2),na.rm=T) < 512)){
    if (return.exports){  return('bcdcor')  }
    stat <- bcdcor(x1, x2)  # bias-corrected dist. corr.
  }else if (method == 'bcdcor'){
    if (return.exports){  return('dcor2d')  }
    stat <- dcor2d(x1, x2, type='U')  # fast bias-corrected dist. corr., a U-statistic
  }else if (method %in% c("pearson", "kendall", "spearman")){
    #         NOTE: Kendall's tau slow & inefficient for large data
    if (return.exports){  return(c('abs', 'cor'))  }
    stat <- abs(cor(x1, x2, method=method))
  }else if ((method == 'calc_corrs') && exists('calc_corrs')){
    if (return.exports){  return(c('abs', 'calc_corrs', 'corr_C'))  }
    stat <- abs(calc_corrs(t(x1), t(x2)))  # compiled, very fast C++ corrs.
  }else{
    if (return.exports){  return(c('abs', 'cor'))  }
    stat <- abs(cor(x1, x2))
  }
  return(stat)
} ##################################################################


##################################################################
similarity_largeData <- function (X, 
                                  dim.subset=dim(X)[1], 
                                  method='dcor', 
                                  parallel=F, verbose=F){
  ##################################################################
  # Re-implementation of similarity_hica() from R library fastHICA,
  #
  #     https://cran.r-project.org/web/packages/fastHICA/index.html
  #
  #   Modified to use fast & memory-efficient distance corr. calculation.
  #
  # Inputs:
  #   X          : data matrix      [obs. x vars.]
  #   dim.subset : (optional) indices for subset of rows of X
  #   method     : similarity method {'dcor', 'bcdcor', 'cor'}
  # Output:
  #   -long-form similarity matrix & used subset of data
  # Requires:
  #   similarity.measure() : calculates similarity & finds dependencies
  #
  
  if (verbose){
    cat(paste0('\n...calculating similarities between  ',dim(X)[2],' vars.  using method: ',method,'\n\n'))
  }
  
  if (dim.subset < 1) {
    stop("dim.subset must be an integer positive value")
  }
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (dim.subset >= n){
    sam <- 1:n
  }else{
    sam <- sample(1:n, dim.subset)
  }
  
  if (parallel){
    sim.fns <- similarity.measure(method=method, 
                                  return.exports=T) # get hidden fn. dependencies
    
    stat_value <- foreach(i=1:(p - 1),
                          .combine=rbind) %:% 
      foreach(j=(i + 1):p, 
              .export=sim.fns,
              .combine=rbind) %dopar% {
                stat  <- similarity.measure(X[sam, i], X[sam, j], method)
                return(c(stat, i, j))
              }
    stat_value <- unname(stat_value)
  }else{
    stat_value <- matrix(0, nrow = ((p + 1) - 1) * ((p) - 1)/2, 
                         ncol = 3)
    count <- 1
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        stat_value[count, 1] <- similarity.measure(X[sam, i], X[sam, j], method)
        stat_value[count, 2] <- i
        stat_value[count, 3] <- j
        count = count + 1
      }
    }
  }
  return(list(similarity_matrix = stat_value, 
              subset = sam, method = method))
} ##################################################################


##################################################################
hierarchicalICA_basis <- function (X,
                                   maxlev=dim(X)[2]-1,
                                   dim.subset=dim(X)[1],
                                   sim.method='dcor',
                                   sim.prev=NA,
                                   selectIC_w_PC1=F,
                                   delta=0.1, pval=0.001, 
                                   parallel=F, verbose=F) {
  ##################################################################
  # Re-implementation of basis_hica() from R library fastHICA,
  #
  #     https://cran.r-project.org/web/packages/fastHICA/index.html
  #
  #   Modified to align & store ICA pre-whitening matrix for back-reconstruction,
  #     fast & memory-efficient distance corr. calculation,
  #     & stats. on eigenvalues to discard levels w/ mergers of redundant sources
  #
  # Inputs:
  #   X          : data matrix      [obs. x vars.]
  #   maxlev     : max. level of the hierarchy, total number of mergers of vars.
  #                 If 'stats' or 'test', max. level set by stats. testing of initial levels
  #                   & current similarity matrix is returned with results
  #   dim.subset : (optional) indices for subset of rows of X
  #   sim.method : similarity measure type {'dcor', 'bcdcor', 'cor', etc.}
  #   sim.prev   : previous similarity matrix, used to construct additional levels of hierarchy,
  #                 output from similarity_largeData()
  #   selectIC_w_PC1 : if T, select IC to keep based on similarity w/ leading PC,
  #                       hierarchy levels then created roughly based on shared variance
  #                       between select ICs at each level
  #   delta      : at each level, test ratio of 2nd eigenvalue to all eigs,
  #                 & flag levels w/ only one eig. is sig. to indicate merger of redundant comps.
  #                   If NA, skip stats. testing.
  #   pval       : p-value for above test, skip stats. testing if NA
  # Output:  list w/ elements:
  #   X            : input data
  #   basis        : matrix w/ vector basis
  #   similarity   : similarity matrix, used to contruct additional levels
  #   aggregation  : matrix with merging info, structured as:
  #                   levels in rows
  #                   ind of 1st comp., used in higher levels, in 1st col.
  #                   ind of 2nd comp., stored & unchanged, in 2nd col.
  #                   similarity value (dist. corr.) in 3rd col.
  #                   variance of 1st new comp. in 4th col.
  #                   variance of 2nd new comp. in 5th col.
  #   prewhitening : applied ICA pre-whitening matrix, or identity 
  #   eigenvalues  : list w/ eigenvalue & stats. test info by hierarchy level
  #   active_inds  : list w/ indices of vars. considered at higher levels of hierarchy
  # Requires:
  #   similarity.measure()    : calculates similarity & finds dependencies
  #   similarity_largeData()  : calculates similarity matrix
  #   eigenValue_ratio_test() : tests relevance of PCs
  #
  
  library(fastICA)
  library(energy)
  
  Xin <- X
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (is.character(maxlev)){
    stopifnot(grepl('stat', maxlev) || grepl('test', maxlev))
    maxlev = p - 1
    maxlev.stop = T  # stop after establishing initial mergers between redundant vars.
  }else {  maxlev.stop = F  }
  if (maxlev < 1 || maxlev > (p - 1)) {
    stop(paste("the maximum level must to be between 1 and ", p - 1, sep = ""))
  }
  if (all(is.na(dim.subset))){
    dim.subset = dim(X)[1]
  }else{
    dim.subset = as.integer(dim.subset)
  }
  
  Ktot <- diag(rep(1, p))
  basis <- vector("list", maxlev)
  prewhit <- vector("list", maxlev)
  eigs <- vector("list", maxlev)
  active_vars <- vector("list", maxlev)
  
  if (!all(is.na(sim.prev))){
    sim_hica <- sim.prev
  }else{
    sim_hica <- similarity_largeData(X, dim.subset, 
                                     method=sim.method,
                                     parallel=parallel, verbose=verbose)
  }
  stat_value <- sim_hica$similarity_matrix
  sam <- sim_hica$subset
  sim.method <- sim_hica$method
  if ('aggregation' %in% names(sim_hica)){
    ind = sim_hica$aggregation
    if (ncol(ind) < 7){ # for prev. hPCA
      d = 7 - ncol(ind)
      ind = cbind(ind, matrix(NA,nrow(ind),d))
    }
  }else{
    ind = NULL
  }
  
  for (h in 1:maxlev) {
    if (all(is.na(stat_value)) || is.null(stat_value) || nrow(stat_value) == 0){
      if (verbose){cat('\n......hierarchy fully constructed, all possible levels exhausted by current & prev. run(s)\n')}
      stat_value = NA
      basis[h:maxlev] = NULL
      prewhit[h:maxlev] = NULL
      eigs[h:maxlev] = NULL
      active_vars[h:maxlev] = NULL
      h = h - 1
      break
    }
    if (verbose){cat(paste0("\nPerforming step ", h, " out of ", maxlev,
                            ", max. similarity = ",signif(max(stat_value[,1]),3)))}
    index1 <- unname(stat_value[which.max(stat_value[, 1]), 2])
    index2 <- unname(stat_value[which.max(stat_value[, 1]), 3])

    if (!is.na(delta) && is.numeric(delta) &&
        !is.na(pval) && is.numeric(pval)){
      
      pca <- eigenValue_ratio_test(X[,c(index1, index2)],
                                   d=delta, m=1, pval=pval, 
                                   return.pca=T, verbose=verbose)
      if (verbose){  cat('\n')  }
      
    }else{
      pca <- prcomp(X[,c(index1, index2)], retx=F, center=T, scale.=T)
      pca$lambda_ratio.test <- F
    }
    pc1 <- pca$x[,1]
    
    if (maxlev.stop && !pca$lambda_ratio.test){
      if (h==1){
        cat('\n\nWARNING: Hierarchical ICA algorithm stopped at 1st level, due to stats. testing criteria!')
        cat('\n  No initial sequence of levels with single eigenvalue found,')
        cat(' recommend raising stats. params., or prevent stopping of algorithm by setting "maxlev" or "h.level" input to NA')
      }else{
        if (verbose){
          cat(paste0('\n...stopping hierarchy construction at level  ',h,'  due to statistical stopping criteria.'))
          cat(paste0('\n      To resume calculations, re-start fn. with input "maxlev" or "h.level" = "resume"'))
        }
      }
      basis[h:maxlev] <- NULL
      prewhit[h:maxlev] <- NULL
      eigs[h:maxlev] <- NULL
      active_vars[h:maxlev] <- NULL
      h <- h - 1
      break     # reached end of initial mergers
    }else if (pca$lambda_ratio.test){  # if H0 is rejected, smaller eigenvalue is irrelevant, merge w/ leading PC
      A <- pca$rotation # eigenvectors of cov. matrix
      K <- diag(c(1,1)) # no pre-whitening 
      pca$lambda_inds <- c(index1, index2) # record indices used in PCA
      
      var1 <- pca$sdev[1]^2
      var2 <- pca$sdev[2]^2
      
    }else{ # merge w/ ICA
      p_sdev <- pca$sdev
      
      if (verbose){
        ica <- fastICA(cbind(X[, index1], X[, index2], deparse.level=0), 2,
                       w.init = cbind(c(p_sdev[1], 0), c(0, p_sdev[2])), verbose=T)
      }else{
        ica <- fastICA(cbind(X[, index1], X[, index2], deparse.level=0), 2,
                       w.init = cbind(c(p_sdev[1], 0), c(0, p_sdev[2])), method='C')
      }
      A <- t(ica$A) # ICA mixing matrix, w/ basis vectors in cols. 2/2 transpose
      K <- ica$K    # ICA pre-whitening matrix, needed for back-reconstruction
      
      var1 <- sum(A[, 1]^2)
      var2 <- sum(A[, 2]^2)
    }
    
    A <- cbind(A[, 1]/sqrt(sum(A[, 1]^2)), A[, 2]/sqrt(sum(A[,2]^2)), deparse.level=0)
    Aparz <- diag(rep(1, p))
    Aparz[c(index1, index2), c(index1, index2)] <- A
    Kparz <- diag(rep(1, p))
    Kparz[c(index1, index2), c(index1, index2)] <- K
    
    prewhit[[h]] <- Ktot %*% Kparz
    basis[[h]] <- Aparz
    eigs[[h]] <- pca[grepl('^lambda*', names(pca))]
    eigs[[h]]$similarity = max(stat_value[,1])
    
    Ktot <- Ktot %*% Kparz %*% solve(t(Aparz))  # store cumulative whitening matrices & inverse basis', for next round
    
    X[, index1] <- X[, index1] - mean(X[, index1])
    X[, index2] <- X[, index2] - mean(X[, index2])
    X <- X %*% Kparz %*% solve(t(Aparz))
    
    sim1_PC1 <- similarity.measure(pc1, X[, index1], sim.method)
    sim2_PC1 <- similarity.measure(pc1, X[, index2], sim.method)
    if (selectIC_w_PC1){  # Select variable to merge into higher levels based on similarity w/ leading PC
      val1 <- sim1_PC1
      val2 <- sim2_PC1
    }else{                # Select variable to merge into higher levels based on variance of new ICs  
      val1 <- var1
      val2 <- var2
    }
    
    if (val1 < val2) {
      ind <- rbind(ind, 
                   c(index2, index1, max(stat_value[,1]), var2, var1, sim2_PC1, sim1_PC1), deparse.level=0)
      index_select <- index2
      
    }else{
      ind <- rbind(ind, 
                   c(index1, index2, max(stat_value[,1]), var1, var2, sim1_PC1, sim2_PC1), deparse.level=0)
      index_select <- index1
      
    }
    
    del <- NULL
    for (cc in 1:dim(stat_value)[1]) {
      if (stat_value[cc, 2] == index1 || stat_value[cc, 2] == index2 ||
          stat_value[cc, 3] == index1 || stat_value[cc, 3] == index2) {
        del <- c(del, cc)
      }
    }
    stat_value <- stat_value[-del, ]
    
    
    # for (f in 1:dim(X)[2]) {      # skip parallelization, more robust & almost as fast
    #   if (f != index_select && sum(f == ind[, 2]) == 0) {
    #     aaa <- similarity.measure(X[sam, f], X[sam, index_select], sim.method)
    #     stat_value <- rbind(stat_value, 
    #                         c(aaa, f, index_select))
    #   }
    # }
    # active_vars[[h]] <- c(unique(as.integer(stat_value[,2])), index_select)
    
    if (parallel){               #  w/ parallelization, may not be worth the effort to speed up
      sim.fns <- similarity.measure(method=sim.method, return.exports=T)  # get dependencies for fn.

      stat_value_updates <- foreach(f=icount(dim(X)[2]),
                                    .export=sim.fns,
                                    .combine=rbind) %dopar% {
                                      if (f != index_select && sum(f == ind[, 2]) == 0) {
                                        aaa <- similarity.measure(X[sam, f], X[sam, index_select], sim.method)
                                        return(c(aaa, f, index_select))
                                      }else{
                                        return(NULL)
                                      }
                                    }
      stat_value <- rbind(stat_value,
                          stat_value_updates, deparse.level=0)
    }else{
      for (f in 1:dim(X)[2]) {
        if (f != index_select && sum(f == ind[, 2]) == 0) {
          aaa <- similarity.measure(X[sam, f], X[sam, index_select], sim.method)
          stat_value <- rbind(stat_value,
                              c(aaa, f, index_select), deparse.level=0)
        }
      }
    }
    active_vars[[h]] <- c(unique(as.integer(stat_value[,2])), index_select)
    
  }
  if ((h == p - 1) || all(is.na(stat_value))){  stat_value = NA  }  # no other hierarchy levels possible
  
  ##################################################################
  return(list(X = Xin, 
              basis = basis, 
              similarity = list('similarity_matrix'=stat_value,
                                'subset'=sam,
                                'method'=sim.method,
                                'aggregation'=ind),
              active_inds = active_vars,
              prewhitening = prewhit, 
              eigenvalues = eigs,
              level.final = h))
} ##################################################################

##################################################################
hierarchicalPCA_basis <- function (X, 
                                   maxlev=dim(X)[2]-1, 
                                   dim.subset=dim(X)[1], 
                                   sim.method='dcor',
                                   sim.prev=NA,
                                   delta=0.1, pval=0.001, 
                                   parallel=F, verbose=F) {
  ##################################################################
  # Modification of basis_hica() from R library fastHICA,
  #
  #     https://cran.r-project.org/web/packages/fastHICA/index.html
  #
  #  Implements the hierarchical Treelets PCA algorithm of Lee et al. (2008),
  #     more flexible similarity measure options,
  #     & applys stats. on eigenvalues to flag levels merging redundant vars.
  #
  # Inputs:
  #   X          : data matrix      [obs. x vars.]
  #   maxlev     : max. level of the hierarchy, total number of mergers of vars.
  #                 If 'stats' or 'test', max. level set by stats. testing of initial levels
  #                   & current similarity matrix is returned with results  
  #   dim.subset : (optional) indices for subset of rows of X
  #   sim.method : similarity measure type {'dcor', 'bcdcor', 'cor', etc.}
  #   sim.prev   : previous similarity matrix, used to construct additional levels of hierarchy,
  #                 output from similarity_largeData()
  #   delta      : at each level, test ratio of 2nd eigenvalue to all eigs,
  #                 & flag levels w/ only one eig. is sig. to indicate merger of redundant comps.
  #                   If NA, skip stats. testing.
  #   pval       : p-value for above test, skip stats. testing if NA
  # Output:  list w/ elements:
  #   X            : input data
  #   basis        : matrix w/ vector basis
  #   similarity   : similarity matrix, used to contruct additional levels
  #   aggregation  : matrix with merging info, structured as:
  #                   levels in rows
  #                   ind of 1st comp., used in higher levels, in 1st col.
  #                   ind of 2nd comp., stored & unchanged, in 2nd col.
  #                   similarity value (dist. corr.) in 3rd col.
  #                   variance of 1st new comp. in 4th col.
  #                   variance of 2nd new comp. in 5th col.
  #   prewhitening : applied ICA pre-whitening matrix, or identity 
  #   eigenvalues  : list w/ eigenvalue & stats. test info by hierarchy level
  #   active_inds  : list w/ indices of vars. considered at higher levels of hierarchy
  # Requires:
  #   similarity.measure()    : calculates similarity & finds dependencies
  #   similarity_largeData()  : calculates similarity matrix
  #   eigenValue_ratio_test() : tests relevance of PCs
  #
  
  library(energy)
  
  Xin <- X
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (is.character(maxlev)){
    stopifnot(grepl('stat', maxlev) || grepl('test', maxlev))
    maxlev <- p - 1
    maxlev.stop <- T  # stop after establishing initial mergers between redundant vars.
  }else {  maxlev.stop <- F  }
  if (maxlev < 1 || maxlev > (p - 1)) {
    stop(paste("the maximum level must to be between 1 and ", p - 1, sep=""))
  }
  if (all(is.na(dim.subset))){
    dim.subset <- dim(X)[1]
  }else{
    dim.subset <- as.integer(dim.subset)
  }
  
  Ktot <- diag(rep(1, p))
  basis <- vector("list", maxlev)
  prewhit <- vector("list", maxlev)
  eigs <- vector("list", maxlev)
  active_vars <- vector("list", maxlev)
  
  if (!all(is.na(sim.prev))){
    sim_hpca <- sim.prev
  }else{
    sim_hpca <- similarity_largeData(X, dim.subset, 
                                     method=sim.method,
                                     parallel=parallel, verbose=verbose)
  }
  stat_value <- sim_hpca$similarity_matrix
  sam <- sim_hpca$subset
  sim.method <- sim_hpca$method
  if ('aggregation' %in% names(sim_hpca)){
    ind = sim_hpca$aggregation
  }else{
    ind = NULL
  }
  
  for (h in 1:maxlev) {
    if (all(is.na(stat_value)) || is.null(stat_value) || nrow(stat_value)==0){
      if (verbose){cat('\n......hierarchy fully constructed, all possible levels exhausted by current & prev. run(s)\n')}
      stat_value <- NA
      basis[h:maxlev] <- NULL
      prewhit[h:maxlev] <- NULL
      eigs[h:maxlev] <- NULL
      active_vars[h:maxlev] <- NULL
      h <- h - 1
      break
    }
    if (verbose){cat(paste0("\nPerforming step ", h, " out of ", maxlev,
                            ", max. similarity = ",signif(max(stat_value[,1]),3)))}
    index1 <- unname(stat_value[which.max(stat_value[, 1]), 2])
    index2 <- unname(stat_value[which.max(stat_value[, 1]), 3])

    if (!is.na(delta) && is.numeric(delta) &&
        !is.na(pval) && is.numeric(pval)){
      
      pca <- eigenValue_ratio_test(X[, c(index1, index2)], 
                                   d=delta, m=1, pval=pval, 
                                   return.pca=T, verbose=verbose)
      if (verbose){  cat('\n')  }
      
    }else{
      pca <- prcomp(X[, c(index1, index2)], retx=F, center=T, scale.=T)
      pca$lambda_ratio.test <- F
    }
    if (maxlev.stop && !pca$lambda_ratio.test){
      if (h==1){
        cat('\n\nWARNING: Hierarchical PCA algorithm stopped at 1st level, due to stats. testing criteria!')
        cat('\n  No initial sequence of levels with single eigenvalue found,')
        cat(' recommend raising stats. params., or prevent stopping of algorithm by setting "maxlev" or "h.level" input to NA')
      }else{
        if (verbose){
          cat(paste0('\n...stopping hierarchy construction at level  ',h-1,'  due to statistical stopping criteria.'))
          cat(paste0('\n      To resume calculations, re-start fn. with input "maxlev" or "h.level" = "resume"'))
        }
      }
      basis[h:maxlev] <- NULL
      prewhit[h:maxlev] <- NULL
      eigs[h:maxlev] <- NULL
      active_vars[h:maxlev] <- NULL
      h <- h - 1
      break     # reached end of initial mergers
    }
    
    A <- pca$rotation  # eigenvectors of cov. matrix
    if (sign(cor(X[,index1], X[,index2])) < 0){  # check alignment of merged vars. in PC1 & PC2
      if (prod(sign(A[,1])) > 0){
        A[1,] <- A[1,] * -1   # correct alignment of index1 relative to index2
      }
    }else{
      if (prod(sign(A[,1])) < 0){
        A[1,] <- A[1,] * -1
      }
    }
    if (!isSymmetric(A)){  # ensure equivalent indexing of rows vs. cols.
      A[,2] <- A[,2] * -1
    }
    
    pca$lambda_inds <- c(index1, index2) # record indices used in PCA
    
    var1 <- pca$sdev[1]^2
    var2 <- pca$sdev[2]^2
    
    Ajacobi <- diag(rep(1, p))
    Ajacobi[c(index1, index2), c(index1, index2)] <- A

    prewhit[[h]] <- Ktot  # store prior PCA whitenings & cumulative rotations from prev. levels
    basis[[h]] <- Ajacobi # basis matrix, consisting of eigenvectors in cols. for merged vars.
    eigs[[h]] <- pca[grepl('^lambda*', names(pca))]  # eigenvalue ratio testing, if any
    eigs[[h]]$similarity <- max(stat_value[,1])      # similarity measure used to select merged vars.
    
    Ktot <- Ktot %*% Ajacobi  # store cumulative basis' projections for next round

    X[, index1] <- X[, index1] - mean(X[, index1])
    X[, index2] <- X[, index2] - mean(X[, index2])
    X <- X %*% Ajacobi

    ind <- rbind(ind, 
                 c(index1, index2, max(stat_value[,1], na.rm=T), var1, var2), deparse.level=0)
    index_select <- index1  # index1 always contains PC1 after above sorting, regardless of ranking & order of inds

    del <- NULL
    for (cc in 1:dim(stat_value)[1]) {
      if (stat_value[cc, 2] == index1 || stat_value[cc, 2] == index2 ||
          stat_value[cc, 3] == index1 || stat_value[cc, 3] == index2) {
        del <- c(del, cc)
      }
    }
    stat_value <- stat_value[-del, ]
    
    
    for (f in 1:dim(X)[2]) {   # no parallelization, simpler, more robust & usually faster
      if (f != index_select && sum(f == ind[, 2]) == 0) {
        aaa <- similarity.measure(X[sam, f], X[sam, index_select], sim.method)
        stat_value <- rbind(stat_value,
                            c(aaa, f, index_select), deparse.level=0)
      }
    }
    active_vars[[h]] <- c(unique(as.integer(stat_value[,2])), index_select)
    
  #   if (parallel){               # w/ parallelization, may not be worth the effort to speed up
  #     sim.fns <- similarity.measure(method=sim.method, return.exports=T)
  # 
  #     stat_value_updates <- foreach(f=icount(dim(X)[2]),
  #                                   .export=sim.fns,
  #                                   .combine=rbind) %dopar% {
  #       if (f != index_select && sum(f == ind[, 2]) == 0) {
  #         aaa <- similarity.measure(X[sam, f], X[sam, index_select], sim.method)
  #         return(c(aaa, f, index_select))
  #       }else{
  #         return(NULL)
  #       }
  #     }
  #     stat_value <- rbind(stat_value,
  #                         stat_value_updates, deparse.level=0)
  #   }else{
  #     for (f in 1:dim(X)[2]) {
  #       if (f != index_select && sum(f == ind[, 2]) == 0) {
  #         aaa <- similarity.measure(X[sam, f], X[sam, index_select], sim.method)
  #         stat_value <- rbind(stat_value,
  #                             c(aaa, f, index_select), deparse.level=0)
  #       }
  #     }
  #   }
  #   active_vars[[h]] <- c(unique(as.integer(stat_value[,2])), index_select)
    
  }

  if ((h == p - 1) || all(is.na(stat_value))){  stat_value <- NA  }  # no other hierarchy levels possible
  
  ##################################################################
  return(list(X = Xin, 
              basis = basis, 
              similarity = list('similarity_matrix'=stat_value,
                                'subset'=sam,
                                'method'=sim.method,
                                'aggregation'=ind),
              active_inds = active_vars,
              prewhitening = prewhit, 
              eigenvalues = eigs,
              level.final = h))
} ##################################################################


##################################################################
edit_hier_dat <- function(hier_file, 
                          prefix=NA,
                          k.indices.rm=NA,
                          verbose=F){
  ##################################################################
  # Edits paused/interrupted group-level hierarchical ICA/PCA,
  #   as output by groupICA_hierarchicalAnalysis(),
  #     before resuming calculations by same fn.
  # 
  # Inputs:
  #   hier_file    : file path to paused groupICA_hierarchicalAnalysis() output
  #   prefix       : prefix for experiment
  #   k.indices.rm : indices of vars. to remove from active vars.,
  #                   effectively excluding from incorporation 
  #                   into higher levels of hierarchy
  #
  
  stopifnot(!all(is.na(hier_file)))
  if (dir.exists(hier_file)){
    hier_suffix = 'hierarchy.RData'   # filesuffix for hierarchical ICA/PCA from groupICA_hierarchicalAnalysis()
    hier_wc = paste0('*', hier_suffix)
    if (!all(is.na(prefix))){  hier_wc = paste0(prefix,'.', hier_wc)  }
    hier_file = list.files(hier_file, hier_wc, full.names=T)[1]
  }
  stopifnot(file.exists(hier_file))
  
  Space.hier = new.env()
  load(hier_file, envir=Space.hier)
  stopifnot(all(c('hier', 'h.levels.criteria', 
                  'data_files', 'prefix', 
                  'similarity.current') %in% names(Space.hier)))
  if(!all(is.na(prefix))){
    stopifnot(prefix == get('prefix', envir=Space.hier))
  }
  
  if (!all(is.na(k.indices.rm))){
    if (verbose){
      cat(paste0('\nExcluding vars. from similarity matrix, preventing inclusion into further levels of hierarchy:'))
      cat(paste0('\n   k = ',paste(k.indices.rm, collapse =',')))
      cat('\n\n')
    }
    
    k.indices.rm = as.integer(k.indices.rm)
    similarity.current = get('similarity.current', envir=Space.hier)
    
    stat_value = similarity.current$similarity_matrix
    ind = similarity.current$aggregation
    
    if (all(is.na(stat_value))|| dim(stat_value)[1] == 0){
      h.levels.criteria = get('h.levels.criteria', envir=Space.hier)
      if (nrow(similarity.current$aggregation) == h.levels.criteria){
        cat(paste('\n...Warning:  cannot exclude indices from subsequent levels of hierarchy, hierarchy is fully constructed!'))
        return(NULL)
      }else{
        cat(paste('\n...Warning: cannot exclude indices from hierarchy, similarity_matrix=',as.character(stat_value),'\n'))
        stopifnot(!all(is.na(stat_value)))
        stopifnot(dim(stat_value)[1] > 0)
      }
    }
    
    ### Remove vars. from similarity matrix list of pairs ###
    del <- NULL
    for (k in k.indices.rm){
      for (cc in 1:dim(stat_value)[1]) {
        if (stat_value[cc, 2] == k || stat_value[cc, 3] == k) {
          del <- c(del, cc)
        }
      }
    }
    if (!is.null(del)){
      stat_value <- stat_value[-del, ]
    }
    stopifnot(!any(stat_value[,2] %in% k.indices.rm))
    stopifnot(!any(stat_value[,3] %in% k.indices.rm))
    
    ### Ensure removed vars. satisfy exclusion condition ###
    K = length(k.indices.rm)
    ind_exclusions = cbind(rep(NA,K), k.indices.rm, rep(0,K), rep(0,K), rep(0,K), 
                           deparse.level=0)
    if (ncol(ind) > ncol(ind_exclusions)){
      d = ncol(ind) - ncol(ind_exclusions)
      ind_exclusions = cbind(ind_exclusions, matrix(NA,K,d),
                             deparse.level=0)
    }
    ind = rbind(ind,
                ind_exclusions, deparse.level=0)
    stopifnot(all(k.indices.rm %in% ind[,2]))
    
    similarity.current$similarity_matrix = stat_value
    similarity.current$aggregation = ind
    assign('similarity.current', similarity.current, envir=Space.hier)
    
    ### Remove vars. from active inds ###
    h.levels.max = get('h.levels.max', envir=Space.hier)
    if (h.levels.max > 0){
      if (get('h.levels.max.hier', envir=Space.hier) == h.levels.max){
        hier = get('hier', envir=Space.hier)
        active_inds = hier[[length(hier)]]$active_inds
        active_inds = active_inds[-which(active_inds %in% k.indices.rm)]
        hier[[length(hier)]]$active_inds = active_inds
        assign('hier', hier, envir=Space.hier)
        
      }else if (get('h.levels.max.hier_ns', envir=Space.hier) == h.levels.max){
        hier_ns = get('hier_ns', envir=Space.hier)
        active_inds = hier_ns[[length(hier_ns)]]$active_inds
        active_inds = active_inds[-which(active_inds %in% k.indices.rm)]
        hier_ns[[length(hier_ns)]]$active_inds = active_inds
        assign('hier_ns', hier_ns, envir=Space.hier)
      
      }
    }
    assign('k.indices.rm', k.indices.rm, envir=Space.hier)
  }
  
  save(file=hier_file, list=ls(envir=Space.hier), envir=Space.hier)
} ##################################################################




##################################################################
groupICA_BackReconstruction <- function(groupICA_path, prefix=NA, 
                                        align.sources=T,
                                        initial_ICA_sources=F,
                                        initial_ns_sources=T,
                                        hier_sources=T,
                                        hier_2nd_sources=F,
                                        hier_levels=NA,
                                        parallel=F, verbose=F){
  ##################################################################
  # Subj.-specific source back-reconstruction for all subjs. in group
  #   Estimates subj-level spatial maps & time series,
  #     using method GICA3 as detailed in:
  #
  # Erhardt et al. (2011). Comparison of multi-subject ICA methods
  #   for analysis of fMRI data. Hum Brain Mapp 32(12), p. 2075-2095.
  #
  # Input:
  #   groupICA_path : Path to dir. w/ saved outputs of 
  #                     subj_PCA(), group_PCA(), group_spatialICA()
  #   prefix        : Identifying prefix to attach to saved output
  #   align.sources : align all ICA source signs to a positive 3rd moment (skew)
  #   initial_ICA_sources      : Reconstruct ICs from original ICA (pre-hier. ICA/PCA)
  #   initial_ns_sources       : Reconstruct non-sig. sources from initial levels of hierarchy,
  #                               summarizing redundnant sources from early mergers of leading PCs
  #   hier_sources             : Reconstruct comps. from hierarchical ICA/PCA,
  #                                 reconstructs single primary comp. at each level
  #   hier_2nd_sources         : Reconstruct secondary comp. from each level,
  #                               not incorporated into subsequent levels.
  #                                 For hPCA, 2nd PC represents difference in var. between ICs
  #                                 For hICA, 2nd IC is an alternative comp. w/ smaller var.
  #   hier_levels              : Reconstruct all comps from select level of hierarchy. Options:
  #                                  "final" : reconstructs all comps from highest level only
  #                                  "all"   : reconstructs all comps from all levels
  #                                  vector of integers : reconstructs specific levels
  # Output: 
  #   saved files w/ *_ICAsources.RData added as suffix & vars.:
  #     S_i  : estimated subj.-specific source spatial maps as [voxels x sources] matrix
  #     R_i  : estimated subj.-specific source time series  as [time x sources] matrix
  #     ICA_model_order       : vector of ICA model order by component (NA for hPCA/hICA comps.)
  #     initial_ICA_inds      : indices of original ICA sources w/n above matrices
  #     hier_inds             : indices of hierarchical ICA/PCA sources w/n above matrices
  #     Hierarchy_level       : vector of levels of the hierarchy for b.r. sources (0=initial ICA data)
  #     hier_k1_inds          : index of primary new comp., created at each level
  #     hier_k2_inds          : index of secondary new comp., created at each level & stored thereafter
  #     hier_ns_inds          : indices of initial/redundant/n.s. sources, summarizing early mergers
  # Requires:
  #   load.data()          : data munging from PCA_fns.R
  #   zeroMean_Yi()        :  "      "      "     "
  #   detrend()            : detrending for pre-processing
  #   subj_preprocData()   : pre-processes data used for subj.-level PCA
  #`` get_groupICA_dat()   : load & format ICA info or multi-order ICA
  #   get_groupPCA_dat()   :   "   &   "   group-level PCA info
  #   get_subjPCA_dat()    :   "   &   "   subj.-level PCA info
  #   subj_BackReconstructICA()     : reconstructs original ICA sources
  #   subj_BackReconstructHier_ns() :      "       sources from non-sig. levels of hPCA/hICA
  #   subj_BackReconstructHier()    :      "       sources from sig. levels of hPCA/hICA
  #   subj_BackReconstructSources() :      "       subj.-specific IC spatial maps & t.s.'
  #
  
  if (verbose){ cat('\nStarting subj.-level source back-reconstruction...')}
  
  ### Inputs ###
  groupICA_path = as.character(groupICA_path) 
  prefix = as.character(prefix) 
  align.sources = as.logical(align.sources)
  initial_ICA_sources = as.logical(initial_ICA_sources)
  hier_sources = as.logical(hier_sources)
  initial_ns_sources = as.logical(initial_ns_sources)
  hier_2nd_sources = as.logical(hier_2nd_sources)
  if (all(is.na(hier_levels)) || all(hier_levels == F)){
    hier_sources_levels = F
    if (!all(is.na(hier_levels))){  hier_levels = NA  }
  }else if (is.character(hier_levels) && (hier_levels %in% c('final', 'last', 'all'))){
    hier_sources_levels = T
  }else if (is.numeric(hier_levels)){
    hier_sources_levels = T
  }else{
    hier_sources_levels = any(as.logical(as.integer(hier_levels)))
    if (is.na(hier_sources_levels)){  hier_sources_levels = F  }
  }
  parallel = as.logical(parallel)
  verbose = as.logical(verbose)
  
  hier_sources_sig = hier_sources_levels || hier_sources || hier_2nd_sources
  hier_sources_nonsig = hier_sources_levels || initial_ns_sources
  hier_sources_any = hier_sources_sig || hier_sources_nonsig
  
  if (verbose){
    if (initial_ICA_sources){ cat('\n...for original/initial ICA sources...')  }
    if (hier_sources_any){ cat('\n...for hierarchical sources...') }
    if (initial_ns_sources){
      cat('\n......including all sources in initial non-sig. levels of hierarchy, after merging redundant vars....')
    }
    if (hier_sources_levels){
      if (is.character(hier_levels) && ((hier_levels == 'final') || (hier_levels == 'last'))){
        cat('\n......including all sources at last (highest) level of hierarchy...')
      }else if (is.character(hier_levels) && (hier_levels == 'all')){
        cat('\n......including redundant hierarchical sources at every level...')
        cat('\n.........(Warning: computationally intensive!')
      }else if ((length(hier_levels) > 0) && is.numeric(hier_levels)){
        cat('\n......including all sources at select level(s):  ')
        cat(paste(hier_levels, collapse = ','))
      }
    }
    if (hier_sources && hier_2nd_sources){
      cat('\n......including new sources (i.e., sum & diff. of variances for hPCA/both ICs for hICA) at each level...')
    }else if (hier_sources){
      cat('\n......including new primary source (i.e., sum of variances for hPCA/main IC for hICA) at each level')
    }
  }
  
  ### Defaults ###
  subj_pca_suffix  = '_PCA1.RData'       # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData'    # filesuffix for group-level temporal PCA from group_PCA()
  group_ica_suffix = 'spatialICA.RData'  # filesuffix for group spatial ICA from group_spatialICA()
  hier_suffix = 'hierarchy.RData' # filesuffix for hierarchical ICA/PCA from groupICA_hierarchicalAnalysis()
  br_suffix   = '_ICAsources.RData' # filesuffix for back-reconstructed subj. sources  
  
  group_ica_wc = group_ica_suffix
  if (!all(is.na(prefix))){  group_ica_wc = paste0(prefix, '.*_', group_ica_wc)  }
  hier_wc = hier_suffix
  if (!all(is.na(prefix))){  hier_wc = paste0(prefix, '.*_', hier_wc)  }
  
  save_path = groupICA_path         # save to same dir.
  
  ### Sanity checks ###
  stopifnot(is.character(groupICA_path) 
            && dir.exists(groupICA_path) 
            && (length(list.files(groupICA_path, group_ica_wc)) > 0))
  if (length(list.files(groupICA_path, group_ica_wc)) > 1){
    cat(paste0('WARNING: multiple ICA runs detected in ', groupICA_path,
               ',\n only the 1st will be processed'))
  }
  if (hier_sources_any){
    stopifnot(length(list.files(groupICA_path, hier_wc)) > 0)
    if (length(list.files(groupICA_path, hier_wc)) > 1){
      cat(paste0('WARNING: multiple hierarchical ICA runs detected in ', groupICA_path,
                 ',\n only the 1st will be processed'))
    }
  }
  
  ### Prelims ###
  if (hier_sources_any){
    hier_file  = list.files(groupICA_path, hier_wc, full.names=T)[1]
  }else{
    hier_file = NA
    group_ica_file = list.files(groupICA_path, group_ica_wc, full.names=T)[1]
  }
  
  if (hier_sources_sig){
    hier = get_hier_dat(hier_file, prefix=prefix, 
                        align.sources=align.sources, 
                        hier_levels=hier_levels,
                        return.input=F,
                        return.sources=align.sources,
                        verbose=verbose)
    if (length(hier) > 0){
      group_ica_file = hier[[1]]$group_ica_file
      h.levels = unlist(lapply(hier, function(yy) return(yy$h.level)))
    }else if (!all(is.na(hier_levels))){  # continue & check for requested levels in n.s. hier levels
      rm('hier')
      hier_sources_sig = F
      hier_sources_nonsig = T
    }else{  stopifnot(length(hier) > 0)  }
  }
  
  if (hier_sources_nonsig){
    hier_ns = get_hier_n.s._dat(hier_file, prefix=prefix, 
                                discard.initial.levels=T,
                                align.sources=align.sources,
                                return.input=F,
                                return.sources=align.sources,
                                verbose=verbose)
    if (length(hier_ns) > 0){
      if (all(is.na(hier_levels)) || !is.character(hier_levels) || !as.logical(hier_levels)){
        hier_ns = hier_ns[1]   # discard higher n.s. levels, if not requested for b.r.
      }
      group_ica_file = hier_ns[[1]]$group_ica_file
      h_ns.levels = unlist(lapply(hier_ns, function(yy) return(yy$h.level)))
    }else if (!all(is.na(hier_levels))){
      rm('hier_ns')
      hier_sources_nonsig = F
    }else{  stopifnot(length(hier_ns) > 0)  }
  }
  
  
  ica = get_groupICA_dat(group_ica_file, prefix,
                         align.sources=align.sources,
                         flatten.multimodel.ICA=F,
                         return.sources=initial_ICA_sources,
                         return.input.data=F,
                         verbose=verbose)
  
  if (all(c('model_orders', 'group_pca_file') %in% names(ica))){
    stopifnot(all(c('K','W','A') %in% names(ica)))
    ica = list(ica)  # format single-order ICA as multi-model order ICA output
  }else{
    stopifnot(all(c('K','W','A') %in% names(ica[[1]])))
  }
  stopifnot('group_pca_file' %in% names(ica[[1]]))
  group_pca_file = ica[[1]]$group_pca_file
  
  if ('model_orders' %in% names(ica[[1]])){
    if (exists('hier') && ('pre_hier_model_orders' %in% names(hier[[1]]))){
      stopifnot(all(hier[[1]]$pre_hier_model_orders == ica[[1]]$model_orders))
    }
    if (exists('hier_ns') && ('pre_hier_model_orders' %in% names(hier_ns[[1]]))){
      stopifnot(all(hier_ns[[1]]$pre_hier_model_orders == ica[[1]]$model_orders))
    }
    ica_model_orders = ica[[1]]$model_orders
  }else if ('model_order' %in% names(ica[[1]])){
    ica_model_orders = ica[[1]]$model_order
  }else{
    ica_model_orders = ncol(ica[[1]]$W)
  }
  
  
  PCA_g = get_groupPCA_dat(group_pca_file, prefix, 
                           model_orders=ica_model_orders, 
                           return.PCs=T,
                           verbose=verbose)
  G    = PCA_g$G
  L_gr = PCA_g$L
  subj_pca_files = PCA_g$subj_pca_files
  stopifnot(is.character(subj_pca_files) 
            && all(file.exists(subj_pca_files)))
  S = length(subj_pca_files)
  
  
  ### Hierarchy level selection ###
  if (!all(is.na(hier_levels)) && hier_sources_any){
    if (!exists('h.levels')){  h.levels = NA  }
    if (!exists('h_ns.levels')){  h_ns.levels = NA  }
    stopifnot(!all(is.na(c(h.levels, h_ns.levels))))
    
    if (is.character(hier_levels) && ((hier_levels == 'final') || (hier_levels == 'last'))){
      hier_levels = max(c(h.levels, h_ns.levels), na.rm=T)  # return all sources at highest/final level of hierarchy
    }else if (is.character(hier_levels) && (hier_levels == 'all')){
      hier_levels = sort(h.levels, h_ns.levels)             # return all sources at each/every level of hierarchy
    }else{
      hier_levels = as.integer(hier_levels)                 # return all sources at select level(s) of hierarchy
    }
    
    if (verbose){
      if (!all(hier_levels %in% c(h_ns.levels, h.levels))){
        cat('\n      WARNING: info for requested hierarchy level(s) not saved, back-reconstruction not possible for these levels:  ')
        cat(paste(hier_levels[!(hier_levels %in% c(h_ns.levels, h.levels))], collapse= ','))
      }
      if (any(hier_levels %in% c(h_ns.levels, h.levels)) && !all(hier_levels %in% h.levels)){
        cat('\n      WARNING: requested hierarchy level(s) are non-sig., prior to back-reconstruction:  ')
        cat(paste(hier_levels[!(hier_levels %in% h.levels)], collapse= ','))
      }
    }
    
    hier_levels = sort(hier_levels[(hier_levels %in% c(h_ns.levels, h.levels))])
    stopifnot(length(hier_levels) > 0)
    
  }else{
    hier_levels = F
  }

  
  ### Calibrate components ###
  if (initial_ICA_sources && align.sources){
    M = length(ica)
    V = nrow(PCA_g$X)
    IC_corrSigns = vector('list', M)
    
    for (m in 1:M){
      K = ica_model_orders[m]
      groupSources = matrix(NA, V, K)
      groupSources = PCA_g$X[,1:K] %*% ica[[m]]$K %*% ica[[m]]$W    # unmix sources from group PCs
      IC_corrSigns[[m]] = rep(1,K)
      
      for (k in 1:K){  # check alignment of unmixed sources w/ saved sources
        IC_corrSigns[[m]][k] = sign(cor(groupSources[,k], ica[[m]]$S[,k])) # correct sign indeterminacy in PCA/ICA, alignment to skew, etc.
      }
      if (verbose && any(IC_corrSigns[[m]] < 0)){
        k.indices = c(1:K)
        cat(paste0('\n...ICA dim. K=',K,',  switching +/- orientation of sources to match un-mixed sources from group PCA space:'))
        cat(paste0('\n      ', paste(k.indices[which(IC_corrSigns[[m]] < 0)], collapse=', ')))
      }
    }
  }else{  IC_corrSigns = NA  }
  
  
  if (hier_sources_nonsig && align.sources){
    stopifnot(length(hier_ns) > 0)
    hL = length(hier_ns)
    V = nrow(PCA_g$X)
    
    hIC_ns_corrSigns = vector('list', hL)
    for (l in 1:hL){
      h = hier[[l]]$h.level
      
      if (any(as.logical(hier_levels))){
        k.inds = c(1:ncol(hier_ns[[l]]$S))
        k.sorting = k.inds
      }else if (initial_ns_sources && (l==1)){
        k.inds = sort(hier_ns[[l]]$active_inds) # only reconstruct active comps, ignoring redundant comps.
        k.sorting = sort.int(hier_ns[[l]]$active_inds, index.return=T)$ix
      }
      if (ncol(hier_ns[[l]]$S) != length(k.inds)){
        stopifnot(ncol(hier_ns[[l]]$S) > length(k.inds))
        stopifnot(ncol(hier_ns[[l]]$S) >= max(k.inds))
        hier_ns[[l]]$S = hier_ns[[l]]$S[,k.inds, drop=F]
      }
      
      hL_K = length(k.inds)
      groupSources = matrix(NA, V, hL_K)
      hIC_ns_corrSigns[[l]] = rep(1,hL_K)  # default to not switching comp. t.s. signs
      stopifnot(ncol(groupSources) == ncol(hier_ns[[l]]$S))
      
      if (('pre_hier_KW' %in% names(hier_ns[[1]])) &&
          (!all(is.na(hier_ns[[1]]$pre_hier_KW)))){
        groupSources = PCA_g$X %*% hier_ns[[1]]$pre_hier_KW %*% hier_ns[[l]]$K %*% hier_ns[[l]]$W[,k.inds, drop=F]  # unmix sources from group PCs
      }else{
        groupSources = PCA_g$X %*% hier_ns[[l]]$K %*% hier_ns[[l]]$W[,k.inds, drop=F]  # unmix sources from group PCs
      }
      
      for (k in 1:hL_K){
        hIC_ns_corrSigns[[l]][k] =  sign(cor(groupSources[,k], hier_ns[[l]]$S[,k]))
      }
      hIC_ns_corrSigns[[l]] = hIC_ns_corrSigns[[l]][k.sorting]  # order to match output of main b.r. fn.
    }
  }else{  hIC_ns_corrSigns = NA  }
  
  
  if (hier_sources_sig && align.sources){
    hL = length(hier)
    V = nrow(PCA_g$X)
    
    hIC_corrSigns = vector('list', hL)
    for (l in 1:hL){
      h = hier[[l]]$h.level
      
      if (any(as.logical(hier_levels))){
        k.inds = c(1:ncol(hier[[l]]$S))
        k.sorting = k.inds
      }else if (hier_sources && hier_2nd_sources){
        k.inds = sort(c(hier[[l]]$k1, hier[[l]]$k2))
        k.sorting = sort.int(c(hier[[l]]$k1, hier[[l]]$k2), index.return=T)$ix
      }else if (hier_sources){
        k.inds = c(hier[[l]]$k1)
        k.sorting = which(sort.int(c(hier[[l]]$k1, hier[[l]]$k2), index.return=T)$ix == 1)  # find k1 ordering in new sources
        hier[[l]]$S = hier[[l]]$S[,k.sorting, drop=F]  # isolate estimated k1 in sources
        k.sorting = 1  # identity index
      }else if (hier_2nd_sources){
        k.inds = c(hier[[l]]$k2)
        k.sorting = which(sort.int(c(hier[[l]]$k1, hier[[l]]$k2), index.return=T)$ix == 2)  # find k2 ordering in new sources
        hier[[l]]$S = hier[[l]]$S[,k.sorting, drop=F]  # isolate estimated k2 in sources
        k.sorting = 1  # identity index
      }
      if (ncol(hier[[l]]$S) != length(k.inds)){
        stopifnot(ncol(hier[[l]]$S) > length(k.inds))
        stopifnot(ncol(hier[[l]]$S) >= max(k.inds))
        hier[[l]]$S = hier[[l]]$S[,k.inds, drop=F]
      }
      
      hL_K = length(k.inds)
      groupSources = matrix(NA, V, hL_K)
      hIC_corrSigns[[l]] = rep(1,hL_K)  # default to not switching IC t.s. signs
      stopifnot(ncol(groupSources) == ncol(hier[[l]]$S))
      
      if (('pre_hier_KW' %in% names(hier[[1]])) &&
          (!all(is.na(hier[[1]]$pre_hier_KW)))){
        groupSources = PCA_g$X %*% hier[[1]]$pre_hier_KW %*% hier[[l]]$K %*% hier[[l]]$W[,k.inds, drop=F]    # unmix sources from group PCs
      }else{
        groupSources = PCA_g$X %*% hier[[l]]$K %*% hier[[l]]$W[,k.inds, drop=F]    # unmix sources from group PCs
      }

      for (k in 1:hL_K){
        hIC_corrSigns[[l]][k] =  sign(cor(groupSources[,k], hier[[l]]$S[,k]))
      }
      hIC_corrSigns[[l]] = hIC_corrSigns[[l]][k.sorting]  # order to match output of main b.r. fn.
    }
  }else{  hIC_corrSigns = NA  }
  
  
  ### Group ICA & PCA Back-reconstruction ###
  vrbs1 = vrbs2 = vrbs3 = verbose
  if (parallel){
    if (verbose){
      for (s in 1:S){  cat(paste0('\n......reconstructing sources for (in parallel):  ', basename(subj_pca_files[s])))  }
    }
    foreach(s=icount(S), 
            .export=c('load.data', 'zeroMean_Yi',
                      'detrend', 'subj_preprocData')) %dopar% {
                        if (verbose){  cat(paste0('\n...Back-reconstructing sources for:  ',
                                                  basename(subj_pca_files[s])))}
                        
                        ### Get subject's PCA ###
                        PCA_s = get_subjPCA_dat(subj_pca_files[s], prefix,
                                                return.PCs=F, return.orig.data=T)
                        
                        ### Back-reconstruction(s) ###
                        br = list('S_i' = NULL, 'R_i'= NULL)
                        if (initial_ICA_sources){
                          br = subj_BackReconstructICA(ica, 
                                                       s, PCA_s, G, L_gr,
                                                       align.sources, IC_corrSigns,
                                                       br, verbose=vrbs1 && vrbs2 && vrbs3)
                          vrbs1 = F
                        }
                        if (hier_sources_nonsig){
                          br = subj_BackReconstructHier_ns(hier_ns, 
                                                           initial_ns_sources, hier_levels,
                                                           s, PCA_s, G, L_gr,
                                                           align.sources, hIC_ns_corrSigns,
                                                           br, verbose=vrbs1 && vrbs2 && vrbs3)
                          vrbs2 = F
                        }
                        if (hier_sources_sig){
                          br = subj_BackReconstructHier(hier,
                                                        hier_sources, hier_2nd_sources,
                                                        hier_levels,
                                                        s, PCA_s, G, L_gr,
                                                        align.sources, hIC_corrSigns,
                                                        br, verbose=vrbs1 && vrbs2 && vrbs3)
                          vrbs3 = F
                        }
                        
                        ### Saving ###
                        save.fname = basename(subj_pca_files[s])
                        if (!is.na(prefix) && !grepl(prefix, save.fname)){
                          save.fname = paste0(prefix, save.fname)
                        }
                        save.fname = sub(subj_pca_suffix, br_suffix, save.fname)
                        save.fname = file.path(save_path, save.fname)
                        if (verbose){  cat(paste0('\n......saving as: ', save.fname))}
                        save(file=save.fname, list=c('br'))
                  }
  }else{
    for (s in 1:S){
      if (verbose){  cat(paste0('\n...Back-reconstructing sources for:  ',
                                basename(subj_pca_files[s])))}
      
      ### Get subject's PCA ###
      PCA_s = get_subjPCA_dat(subj_pca_files[s], prefix,
                              return.PCs=F, return.orig.data=T)
      
      ### Back-reconstruction(s) ###
      br = list('S_i' = NULL, 'R_i'= NULL)
      if (initial_ICA_sources){
        br = subj_BackReconstructICA(ica, 
                                     s, PCA_s, G, L_gr,
                                     align.sources, IC_corrSigns,
                                     br, verbose=vrbs1 && vrbs2 && vrbs3)
        vrbs1 = F
      }
      if (hier_sources_nonsig){
        br = subj_BackReconstructHier_ns(hier_ns, 
                                         initial_ns_sources, hier_levels,
                                         s, PCA_s, G, L_gr,
                                         align.sources, hIC_ns_corrSigns,
                                         br, verbose=vrbs1 && vrbs2 && vrbs3)
        vrbs2 = F
      }
      if (hier_sources_sig){
        br = subj_BackReconstructHier(hier,
                                      hier_sources, hier_2nd_sources,
                                      hier_levels,
                                      s, PCA_s, G, L_gr,
                                      align.sources, hIC_corrSigns,
                                      br, verbose=vrbs1 && vrbs2 && vrbs3)
        vrbs3 = F
      }
      
      ### Saving ###
      save.fname = basename(subj_pca_files[s])
      if (!is.na(prefix) && !grepl(prefix, save.fname)){
        save.fname = paste0(prefix, save.fname)
      }
      save.fname = sub(subj_pca_suffix, br_suffix, save.fname)
      save.fname = file.path(save_path, save.fname)
      if (verbose){  cat(paste0('\n......saving as: ', save.fname))}
      save(file=save.fname, list=c('br'))
    }
  }
  if (verbose){  cat(paste0('\n...saved in dir.:  ', save_path, '\n\n'))  }
} ##################################################################

##################################################################
subj_BackReconstructICA <- function(ica, 
                                    s, PCA_s, G, L_gr=NA, 
                                    align.sources=T, IC_corrSigns=NA,
                                    br=list(), verbose=F){
  ##################################################################
  # Subj.-specific source back-reconstruction for ICA
  #     using 'subj_BackReconstructSources()'
  #       called as part of groupICA_BackReconstruction().
  #
  # Input:
  #   ica   : ICA output, from get_groupICA_dat()
  #   s     : subject index
  #   PCA_s : subj.-specific PCA data, from get_subjPCA_dat()
  #   G     : group PCA eigenvector matrix, from get_groupPCA_dat()
  #   L_gr  : group PCA eigenvalues, from get_groupPCA_dat
  #   align.sources : align sources to skew
  #   IC_corrSigns  : option list of vectors w/ IC signs for above
  #   br    : previous back-reconstructed components for subj.
  # Output:     list w/ elements
  #   S_i              : matrix of source spatial maps [voxels x comps.]
  #   R_i              : matrix of source time series  [time x comps.]
  #   ICA_model_order  : vector of ICA model orders
  #   initial_ICA_inds : vector of indices for ICA components in S_i & R_i
  #   orig_var_inds    : vector of original indices of comps. in S_i & R_i
  #   Hierarchy_level  : hierarchical ICA/PCA level for return comps.
  #   hier_algorithm   : hierarchical ICA/PCA algorithm used for level
  #   hier_ns_inds     : indices of n.s. levels of the hierarchy
  #   hier_k1_inds     : indices of higher var. comps. at each level of hierarchy
  #   hier_k2_inds     : indices of lower var. comps. at each level of hierarchy 
  #   hier_ns_inds     : indices for all non-sig. comps. from hICA / hPCA
  #   hier_inds        : indices for all sig. comps. from hICA / hPCA
  # Requires:
  #   subj_BackReconstructSources()
  #   zeroMean_Yi()    : datamunging fn. from PCA_fns.R
  #
  
  ### Inputs / defaults ###
  br = as.list(br)
  if ('model_orders' %in% names(ica[[1]])){
    ica_model_orders = ica[[1]]$model_orders
  }else if ('model_order' %in% names(ica[[1]])){
    ica_model_orders = ica[[1]]$model_order
  }else{        # assumes output from single-model order ICA
    ica_model_orders = ncol(ica[[1]]$W)
  }

  if (!('S_i' %in% names(br)) || is.null(br$S_i)){
    initial_ICA_i0 = 1 # start index of initial ICA sources in output
  }else{
    initial_ICA_i0 = ncol(br$S_i) + 1
    if (!('ICA_model_order' %in% names(br))){  br$ICA_model_order = rep(NA, ncol(br$S_i))  }
    if (!('orig_var_inds' %in% names(br))){    br$orig_var_inds = rep(NA, ncol(br$S_i))  }
  }
  

  for (m in 1:length(ica)){
    K = ica_model_orders[m]
    if (m > 1){  
      g1 = sum(ica_model_orders[1:(m-1)]) + 1
    }else{  g1 = 1  }
    g2 = g1 + ica_model_orders[m] - 1
    
    br.origICA = subj_BackReconstructSources(ica[[m]], 
                                             G[,g1:g2], PCA_s$U_i, PCA_s$Y, 
                                             L_gr[g1:g2], PCA_s$L_i, 
                                             pre_ica=NA, s=s, 
                                             verbose=verbose)
    
    if (align.sources && any(IC_corrSigns[[m]] < 0)){  # Switch sign, after back-reconstruction
      br.origICA$R_i = br.origICA$R_i %*% diag(IC_corrSigns[[m]])
    }
    br$S_i = cbind(br$S_i, br.origICA$S_i, deparse.level=0)
    br$R_i = cbind(br$R_i, br.origICA$R_i, deparse.level=0)
    br$ICA_model_order = c(br$ICA_model_order, rep(K, K))
  }
  
  br$initial_ICA_inds = c(br$initial_ICA_inds, c(initial_ICA_i0:ncol(br$S_i)))
  br$Hierarchy_level = c(br$Hierarchy_level, rep(0, ncol(br$S_i))) # initial sources at level 0 of hierarchy
  br$orig_var_inds   = c(br$orig_var_inds, c(initial_ICA_i0:ncol(br$S_i)))
  
  ##################################################################
  return(br)
} ##################################################################

##################################################################
subj_BackReconstructHier_ns <- function(hier_ns, 
                                        initial_ns_sources=T,
                                        hier_levels=F, 
                                        s, PCA_s, G, L_gr=NA, 
                                        align.sources=T, hIC_ns_corrSigns=NA,
                                        br=list(), verbose=F){
  ##################################################################
  # Subj.-specific source back-reconstruction
  #   for non-sig. levels of group hierarchical ICA / hierarchical PCA, 
  #     using 'subj_BackReconstructSources()'
  #       called as part of groupICA_BackReconstruction().
  #
  # Input:
  #   hier_ns     : list w/ hierarchy info by level, from get_hier_n.s._dat()
  #   initial_ns_sources : reconstruct sources in last of initial n.s. hierarchy levels,
  #                         (~all sources used for higher levels in hierarchy)  
  #   hier_levels : vector of specific levels of hierarchy to back-reconstruct
  #   s     : subject index
  #   PCA_s : subj.-specific PCA data, from get_subjPCA_dat()
  #   G     : group PCA eigenvector matrix, from get_groupPCA_dat()
  #   L_gr  : group PCA eigenvalues, from get_groupPCA_dat
  #   align.sources : align sources to skew
  #   hIC_ns_corrSigns  : option list of vectors w/ IC signs for above
  #   br    : previous back-reconstructed components for subj.
  # Output:     list w/ elements
  #   S_i              : matrix of source spatial maps [voxels x comps.]
  #   R_i              : matrix of source time series  [time x comps.]
  #   ICA_model_order  : vector of ICA model orders
  #   initial_ICA_inds : vector of indices for ICA components in S_i & R_i
  #   orig_var_inds    : vector of original indices of comps. in S_i & R_i
  #   Hierarchy_level  : hierarchical ICA/PCA level for return comps.
  #   hier_algorithm   : hierarchical ICA/PCA algorithm used for level
  #   hier_ns_inds     : indices of n.s. levels of the hierarchy
  #   hier_k1_inds     : indices of higher var. comps. for each level of hierarchy
  #   hier_k2_inds     : indices of lower var. comps. for each level of hierarchy 
  #   hier_ns_inds     : indices for all non-sig. comps. from hICA / hPCA
  #   hier_inds        : indices for all sig. comps. from hICA / hPCA
  # Requires:
  #   subj_BackReconstructSources()
  #   zeroMean_Yi()    : datamunging fn. from PCA_fns.R
  #
  
  ### Inputs / defaults ###
  br = as.list(br)
  
  if (!('S_i' %in% names(br)) || is.null(br$S_i)){
    hier_ns_i0 = 1  # start index of hierarchical sources in output
  }else{
    hier_ns_i0 = ncol(br$S_i) + 1  # start index of hPCA/hICA sources in output
    if (!('Hierarchy_level' %in% names(br))){  br$Hierarchy_level = rep(NA, ncol(br$S_i))  }
    if (!('hier_algorithm' %in% names(br))){   br$hier_algorithm = rep(NA, ncol(br$S_i))  }
    if (!('orig_var_inds' %in% names(br))){    br$orig_var_inds = rep(NA, ncol(br$S_i))  }
    if (!('hier_ns_inds' %in% names(br))){     br$hier_ns_inds = rep(NA, ncol(br$S_i))  }
  }
  

  for (l in 1:length(hier_ns)){ # index hierarchical sources starting from bottom of hierarchy
    h = hier_ns[[l]]$h.level
    hL_K = ncol(hier_ns[[l]]$W)
    k.inds = NULL
    if (any(as.logical(hier_levels)) && (h %in% hier_levels)){
      k.inds = c(1:hL_K)
    }else if (initial_ns_sources && (l==1)){
      k.inds = hier_ns[[l]]$active_inds
    }
    if (is.null(k.inds)){
      cat(paste0('\nWARNING:  could not determine var. indices for n.s. level l = ',l))
      next
    }
    
    br.hier_ns = subj_BackReconstructSources(hier_ns[[l]], 
                                             G, PCA_s$U_i, PCA_s$Y, 
                                             L_gr, PCA_s$L_i, 
                                             pre_ica=hier_ns[[1]]$pre_hier_KW, 
                                             s=s, k.indices=k.inds,
                                             verbose=verbose)
    
    if (align.sources && any(hIC_ns_corrSigns[[l]] < 0, na.rm=T)){  # switch sign, after back-reconstruction
      br.hier_ns$R_i = br.hier_ns$R_i %*% diag(hIC_ns_corrSigns[[l]])
    }
    br$S_i = cbind(br$S_i, br.hier_ns$S_i, deparse.level=0)
    br$R_i = cbind(br$R_i, br.hier_ns$R_i, deparse.level=0)
    br$Hierarchy_level = c(br$Hierarchy_level, rep(h, length(k.inds)))
    br$hier_algorithm  = c(br$hier_algorithm, rep(hier_ns[[l]]$algorithm, length(k.inds)))
    br$orig_var_inds   = c(br$orig_var_inds,  k.inds)
    
  }
  br$hier_ns_inds = c(br$hier_ns_inds, c(hier_ns_i0:ncol(br$S_i)))
  
  ##################################################################
  return(br)
} ##################################################################


##################################################################
subj_BackReconstructHier <- function(hier, 
                                     hier_sources=T, hier_2nd_sources=T,
                                     hier_levels=F, 
                                     s, PCA_s, G, L_gr=NA, 
                                     align.sources=T, hIC_corrSigns=NA,
                                     br=list(), verbose=F){
  ##################################################################
  # Subj.-specific source back-reconstruction
  #   for all sig. levels of group hierarchical ICA / hierarchical PCA, 
  #     using 'subj_BackReconstructSources()'
  #       called as part of groupICA_BackReconstruction().
  #
  # Input:
  #   hier     : list w/ hierarchy info by level, from get_hier_dat()
  #   hier_sources      : reconstruct primary comp. from each level
  #   hier_2nd_sources  : reconstruct secondary comp. from each level
  #   hier_levels       : vector of specific levels of hierarchy to back-reconstruct
  #   s     : subject index
  #   PCA_s : subj.-specific PCA data, from get_subjPCA_dat()
  #   G     : group PCA eigenvector matrix, from get_groupPCA_dat()
  #   L_gr  : group PCA eigenvalues, from get_groupPCA_dat
  #   align.sources : align sources to skew
  #   nIC_ns_corrSigns  : option list of vectors w/ IC signs for above
  #   br    : previous back-reconstructed components for subj.
  # Output:     list w/ elements
  #   S_i              : matrix of source spatial maps [voxels x comps.]
  #   R_i              : matrix of source time series  [time x comps.]
  #   ICA_model_order  : vector of ICA model orders
  #   initial_ICA_inds : vector of indices for ICA components in S_i & R_i
  #   orig_var_inds    : vector of original indices of comps. in S_i & R_i
  #   Hierarchy_level  : hierarchical ICA/PCA level for return comps.
  #   hier_algorithm   : hierarchical ICA/PCA algorithm used for level
  #   hier_ns_inds     : indices of n.s. levels of the hierarchy
  #   hier_k1_inds     : indices of higher var. comps. for each level of hierarchy
  #   hier_ns_inds     : indices for all non-sig. comps. from hICA / hPCA
  #   hier_inds        : indices for all sig. comps. from hICA / hPCA
  # Requires:
  #   subj_BackReconstructSources()
  #   zeroMean_Yi()    : datamunging fn. from PCA_fns.R
  #
  
  ### Inputs / defaults ###
  br = as.list(br)
  
  if (!('S_i' %in% names(br)) || is.null(br$S_i)){
    hier_i0 = 1  # start index of hierarchical ICA/PCA sources in output
  }else{
    hier_i0 = ncol(br$S_i) + 1  # start index of hierarchical ICA/PCA sources in output
    if (!('Hierarchy_level' %in% names(br))){  br$Hierarchy_level = rep(NA, ncol(br$S_i))  }
    if (!('hier_algorithm' %in% names(br))){   br$hier_algorithm = rep(NA, ncol(br$S_i))  }
    if (!('orig_var_inds' %in% names(br))){    br$orig_var_inds = rep(NA, ncol(br$S_i))  }
  }

  for (l in 1:length(hier)){
    h = hier[[l]]$h.level
    hL_K = ncol(hier[[l]]$W)
    k.inds = NULL
    if (any(as.logical(hier_levels)) && (h %in% hier_levels)){
      k.inds = c(1:hL_K)
      if (!('S_i' %in% names(br)) || is.null(br$S_i)){
        br$hier_k1_inds = c(br$hier_k1_inds, hier[[l]]$k1)
        br$hier_k2_inds = c(br$hier_k2_inds, hier[[l]]$k2)
      }else{
        br$hier_k1_inds = c(br$hier_k1_inds, ncol(br$S_i) + hier[[l]]$k1)
        br$hier_k2_inds = c(br$hier_k2_inds, ncol(br$S_i) + hier[[l]]$k2)
      }
    }else if (hier_sources || hier_2nd_sources){
      if (hier_sources){      k.inds = c(k.inds, hier[[l]]$k1)  }
      if (hier_2nd_sources){  k.inds = c(k.inds, hier[[l]]$k2)  }
      if (hier_sources && !hier_2nd_sources){
        if (!('S_i' %in% names(br)) || is.null(br$S_i)){
          br$hier_k1_inds = c(br$hier_k1_inds, 1)
        }else{
          br$hier_k1_inds = c(br$hier_k1_inds, ncol(br$S_i) + 1)
        }
      }else if (hier_2nd_sources && !hier_sources){
        if (!('S_i' %in% names(br)) || is.null(br$S_i)){
          br$hier_k2_inds = c(br$hier_k2_inds, 1)
        }else{
          br$hier_k2_inds = c(br$hier_k2_inds, ncol(br$S_i) + 1)
        }
      }else if (hier_sources && hier_2nd_sources){
        if (!('S_i' %in% names(br)) || is.null(br$S_i)){
          br$hier_k1_inds = c(br$hier_k1_inds, 1)
          br$hier_k2_inds = c(br$hier_k2_inds, 2)
        }else{
          br$hier_k1_inds = c(br$hier_k1_inds, ncol(br$S_i) + 1)
          br$hier_k2_inds = c(br$hier_k2_inds, ncol(br$S_i) + 2)
        }
      }
    }
    if (is.null(k.inds)){
      cat(paste0('\nWARNING:  could not determine var. indices for level l = ',l))
      next
    }
    
    br.hier = subj_BackReconstructSources(hier[[l]], 
                                          G, PCA_s$U_i, PCA_s$Y, 
                                          L_gr, PCA_s$L_i, 
                                          pre_ica=hier[[1]]$pre_hier_KW, 
                                          s=s, k.indices=k.inds,
                                          verbose=verbose)
    
    if (align.sources && any(hIC_corrSigns[[l]] < 0, na.rm=T)){  # switch sign, after back-reconstruction
      br.hier$R_i = br.hier$R_i %*% diag(hIC_corrSigns[[l]])
    }
    br$S_i = cbind(br$S_i, br.hier$S_i, deparse.level=0)
    br$R_i = cbind(br$R_i, br.hier$R_i, deparse.level=0)
    br$Hierarchy_level = c(br$Hierarchy_level, rep(h, length(k.inds)))
    br$hier_algorithm  = c(br$hier_algorithm, rep(hier[[l]]$algorithm, length(k.inds)))
    br$orig_var_inds   = c(br$orig_var_inds, k.inds)
  }
  br$hier_inds = c(br$hier_inds, c(hier_i0:ncol(br$S_i)))
  
  ##################################################################
  return(br)
} ##################################################################

##################################################################
subj_BackReconstructSources <- function(ica, 
                                        pca_group, pca_subj, data_subj, 
                                        lambda_group=NA, lambda_subj=NA,
                                        pre_ica=NA,
                                        s=NA, k.indices=NA,
                                        verbose=T){
  ##################################################################
  # Subj.-specific source back-reconstruction
  #   Estimates subj-level spatial maps & time series,
  #     using method GICA3 as detailed in:
  #
  # Erhardt et al. (2011). Comparison of multi-subject ICA methods
  #   for analysis of fMRI data. Hum Brain Mapp 32(12), p. 2075-2095.
  #
  # Input:
  #   ica : list w/ group ICA output, as formatted by fastICA(), w/ elements
  #     K : pre-whitening matrix for ICA
  #     W : est. unmixing matrix, w/ pos. skew
  #     A : mixing matrix, inverse of W
  #     S : source matrix, est. by projecting data onto unmixing matrix
  #   pca_group : group-level eigenvectors      as [(subjs*K1) x K2] matrix
  #   pca_subj  : subj.-level eigenvector basis as [time x K1] matrix
  #   data_subj : subj. data                    as [voxels x time] matrix
  #   lambda_group : group-level eigenvalues
  #   lambda_subj  : subject-level eigenvalues
  #   pre_ica      : matrix of pre-hierarchical ICA processing
  #                   (e.g. pre-whitening & unmixing if ICA sources were used for hPCA/hICA)
  #   s         :  (optional) subj. numeric index
  #   k.indices : (optional) vector of integers to limit b.r.
  #                 to subset of sources
  # Output: list w/ elements
  #   S_itilde  : estimated subj.-specific source spatial maps as [voxels x K] matrix
  #   R_itilde  : estimated subj.-specific source time series  as [time x K] matrix
  # Requires:
  #   zeroMean_Yi() : datamunging fn. from PCA_fns.R
  #
  
  library(MASS) # ginv() fn. for generalized inverse (pseudo-inverse) of a matrix
  
  ### Input handling & Defaults ###
  s = as.integer(s)
  L_gr = as.numeric(lambda_group)
  L_i = as.numeric(lambda_subj)
  
  stopifnot(is.list(ica) && all(c('W', 'A', 'K') %in% names(ica)))
  
  if (all(is.na(k.indices))){
    K = nrow(ica$A)    # number of sources estimated
    k.indices = c(1:K)
  }else{
    if (is.logical(k.indices)){
      K = sum(k.indices, na.rm=T)
    }else{
      k.indices = as.integer(k.indices)
      K = length(k.indices)
    }
  }
  
  stopifnot(is.matrix(pca_subj) && is.numeric(pca_subj))
  U_i = pca_subj
  K1 = ncol(U_i)
  if (all(is.na(L_i))){   L_i  = rep(1,K1)  }
  
  Y_i = zeroMean_Yi(data_subj)
  
  stopifnot(is.matrix(pca_group) && is.numeric(pca_group))
  G = pca_group
  if (!is.na(s)){ # find subj.specific column space w/n group PCA space
    i.start = K1 * (s-1) + 1
    i.end =   K1 * s
    G_i = G[i.start:i.end,] # subj.-specific slice of group reduction matrix
  }else{
    G_i = G
  }
  K2 = ncol(G)
  if (all(is.na(L_gr))){  L_gr = rep(1,K2)  }
  
  ### Additional modifications of subj. data by multiplying, pre-ICA/hPCA/hICA ###
  pre_ica = as.matrix(pre_ica)
  if (all(is.na(pre_ica))){  pre_ica = diag(rep(1,K2))  }
  
  ### Sanity checks ###
  stopifnot(is.matrix(data_subj) && is.numeric(data_subj))
  stopifnot(ncol(ica$K)   == nrow(ica$W))
  stopifnot(ncol(pre_ica) == nrow(ica$K))
  stopifnot(ncol(G_i) == nrow(pre_ica))
  stopifnot(ncol(U_i) == nrow(G_i))
  stopifnot(ncol(Y_i) == nrow(U_i))
  
  
  ### Back-reconstruction ###
  if (verbose){
    cat('\n......projecting ICA unmixing matrix onto subj. data to obtain subj.-specific spatial maps')
  }
  S_itilde = Y_i %*% U_i %*% diag(1 / sqrt(L_i)) %*%
    G_i %*% diag(1 / sqrt(L_gr)) %*%
    pre_ica %*% ica$K %*% ica$W[,k.indices, drop=F]  # Erhardt et al. (2011) eq. 12

  if (verbose){
    cat('\n......projecting ICA mixing matrix onto PCA space to obtain subj.-specific time series')
  }
  R_itilde = U_i %*% diag(1 / sqrt(L_i)) %*%
    ginv(t(G_i %*% diag(1 / sqrt(L_gr)))) %*%
    pre_ica %*% ica$K %*% t(ica$A[k.indices,, drop=F])  # Erhardt et al. (2011) eq. 14
  
  
  ##################################################################
  return(list('S_i' = S_itilde,  # estimated subj.-specific source spatial maps [voxels x sources] matrix
              'R_i' = R_itilde)) # estimated subj.-specific source time series  [time x sources] matrix
} ##################################################################


##################################################################
groupICA_ScaleSources <- function(groupICA_path, prefix=NA, 
                                  scale='z-scores',
                                  save.as.nii=T,
                                  parallel=F, verbose=F){
  ##################################################################
  # Calibrates & scales subj.-specific sources signals
  #   following group ICA & back-reconstruction
  #
  # Input:
  #   groupICA_path : path to dir. w/ saved outputs of 
  #                     subj_PCA(), group_PCA(), group_spatialICA()
  #   prefix : identifying prefix to attach to saved output
  #   scale  : scaling for subj. source signals, currently only z-scores
  #   save.as.nii : save output as nifti volume
  # Output: if save.as.nii
  #   -saved 4D nifti vol. w/ scaled & calibrated spatial maps for each subj.
  #   -saved matrix as 2D nifti vol. w/ time series for each subj.
  # Requies:
  #   unflatten_img() to load ICA spatial maps into 4D nifti vol.
  #

  library(RNifti)
  
  if (verbose){ cat('\nScaling subj.-level ICA sources...')}
  
  ### Inputs & defaults ###
  subj_pca_suffix  = '_PCA1.RData'       # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData'    # filesuffix for group-level temporal PCA from group_PCA()
  group_ica_suffix = 'spatialICA.RData'  # filesuffix for group spatial ICA from group_spatialICA()
  br_suffix        = '_ICAsources.RData' # filesuffix for back-reconstructed subj. sources
  spatial_maps_suffix = '_ICAspatialMaps.nii' # filesuffix for 4D nii vol. w/ back-reconstructed spatial maps
  time_series_suffix  = '_ICAtimeseries.nii'  # filesuffix for 2D nii vol. w/ back-reconstructed t.s.
  
  groupICA_path = as.character(groupICA_path)
  prefix = as.character(prefix)
  scale = as.character(scale)
  if (length(scale) > 1){  scale = scale[1]  }
  if (!(scale %in% c('z-scores'))){  scale = 'z-scores'  }
  save.as.nii = as.logical(save.as.nii)
  
  ### Prelims ###
  stopifnot(scale %in% c('z-scores'))
  group_ica_wc = paste0('*', group_ica_suffix)
  br_wc        = paste0('*', br_suffix)
  if (!all(is.na(prefix))){
    group_ica_wc = paste0(prefix,'.', group_ica_wc)
    br_wc = paste0(prefix,'.', br_wc)
  }

  ### Get file names & paths ###
  group_ica_file = list.files(groupICA_path, group_ica_wc, full.names=T)[1]
  br_files = list.files(groupICA_path, br_wc, full.names=T)
  S = length(br_files)
  
  
  ### Load ICA, group PCA & subj. PCA data ###
  Space.ica = new.env()
  load(group_ica_file, envir=Space.ica)
  stopifnot(all(c('ica', 'K', 'data_files', 'prefix') %in% names(Space.ica)))
  if (all(is.na(prefix))){
    prefix = get('prefix', envir=Space.ica)
  }else{
    stopifnot(prefix == get('prefix', envir=Space.ica))
  }
  group_pca_files = get('data_files', envir=Space.ica)
  
  
  Space.pca = new.env()
  stopifnot(is.character(group_pca_files) 
            && file.exists(group_pca_files))
  load(group_pca_files, envir=Space.pca)
  stopifnot(all(c('pca', 'data_files', 
                  'K2', 'prefix') %in% names(Space.pca)))
  stopifnot(prefix == get('prefix', envir=Space.pca))
  subj_pca_files = get('data_files', envir=Space.pca)
  
  
  Space.subj_pca = new.env()
  stopifnot(is.character(subj_pca_files) 
            && file.exists(subj_pca_files))
  load(subj_pca_files[1], envir=Space.subj_pca)
  stopifnot(all(c('mask', 'prefix') %in% names(Space.subj_pca)))
  stopifnot(prefix == get('prefix', envir=Space.subj_pca))
  mask_file = get('mask', envir=Space.subj_pca)
  
  if (!save.as.nii){
    spatial_maps_suffix = sub('.nii$', '.RData', spatial_maps_suffix)
    time_series_suffix = sub('.nii$', '.RData', time_series_suffix)
  }
  
  
  ### Scale back-reconstructed components ###
  if (parallel){
    foreach(s=icount(S), .export=c('unflatten_img')) %dopar% {
      if (verbose){  cat(paste0('\n...scaling subj.:  ', br_files[s]))}
      load(br_files[s])
      if (all(is.na(br$S_i))){
        if (verbose){  cat(paste0('\n......Warning:  no back-reconstructed to scale for ', basename(br_files[s])))  }
        next
      }
      
      if (scale == 'z-scores'){
        if (verbose){  cat('\n...subtracting mean & scaling to s.d. for all sources  ~  z-scores')}
        br$S_i = sweep(br$S_i, 2, colMeans(br$S_i))
        br$S_i = apply(br$S_i, 2, function(yy) return(yy / sd(yy)))
        br$R_i = sweep(br$R_i, 2, colMeans(br$R_i))
        br$R_i = apply(br$R_i, 2, function(yy) return(yy / sd(yy)))
      }
      if (save.as.nii){
        spatialmaps.nii = unflatten_img(br$S_i, mask_file)
        timeseries.nii = asNifti(br$R_i)
      }

      save.fname1 = sub(br_suffix, spatial_maps_suffix, br_files[s])
      if (!is.na(prefix) && !grepl(prefix, save.fname1)){
        save.fname = paste0(prefix, save.fname1)
      }
      save.fname2 = sub(br_suffix, time_series_suffix, br_files[s])
      if (!is.na(prefix) && !grepl(prefix, save.fname2)){
        save.fname = paste0(prefix, save.fname2)
      }
      if (verbose){  cat(paste0('\n......& saving as:  ', save.fname1))  }
      if (save.as.nii){
        if (verbose){  cat(paste0('\n......& saving as:  ', save.fname2))  }
        writeNifti(spatialmaps.nii, save.fname1)
        writeNifti(timeseries.nii, save.fname2)
      }else{
        save(list='br', file=save.fname1)
      }
    }
  }else{
    for (s in 1:S){
      if (verbose){  cat(paste0('\n...scaling subj.:  ', basename(br_files[s])))}
      load(br_files[s])
      if (all(is.na(br$S_i))){
        if (verbose){  cat('\n......Warning:  no back-reconstructed to scale for this subject')}
        next
      }
      
      if (scale == 'z-scores'){
        if (verbose){  cat('\n......subtracting mean & scaling to s.d. for all sources  ~  z-scores')}
        br$S_i = sweep(br$S_i, 2, colMeans(br$S_i))
        br$S_i = apply(br$S_i, 2, function(yy) return(yy / sd(yy)))
        br$R_i = sweep(br$R_i, 2, colMeans(br$R_i))
        br$R_i = apply(br$R_i, 2, function(yy) return(yy / sd(yy)))
      }
      if (save.as.nii){
        spatialmaps.nii = unflatten_img(br$S_i, mask_file)
        timeseries.nii = asNifti(br$R_i)
      }

      save.fname1 = sub(br_suffix, spatial_maps_suffix, br_files[s])
      if (!is.na(prefix) && !grepl(prefix, save.fname1)){
        save.fname = paste0(prefix, save.fname1)
      }
      save.fname2 = sub(br_suffix, time_series_suffix, br_files[s])
      if (!is.na(prefix) && !grepl(prefix, save.fname2)){
        save.fname = paste0(prefix, save.fname2)
      }
      if (verbose){  cat(paste0('\n.........& saving as:  ', basename(save.fname1)))  }
      if (save.as.nii){
        if (verbose){  cat(paste0('\n.........& saving as:  ', basename(save.fname2)))  }
        writeNifti(spatialmaps.nii, save.fname1)
        writeNifti(timeseries.nii, save.fname2)
      }else{
        save(list='br', file=save.fname1)
      }
    }
  }
  if (verbose){  cat(paste0('\n...saved in dir.:  ', unique(dirname(br_files)), '\n\n'))  }
} ##################################################################


##################################################################
groupICA_summaryStats <- function(groupICA_path, prefix=NA, 
                                  calc.var=F,
                                  create.csv=T,
                                  verbose=F){
  ##################################################################
  # Calculates summary statistics
  #   & saves as Nifti volume
  #   for group-level spatial ICA output
  #
  # Input:
  #   groupICA_path : path to dir. w/ saved outputs of 
  #                     subj_PCA(), group_PCA(), group_spatialICA()
  #   prefix   : identifying prefix to attach to saved output
  #   calc.var : calculate standard deviation for each component
  #   create.csv : creates .csv table w/ component info in fields:
  #
  # Output: 
  #   -saved 4D nifti vol. w/ mean of subj. spatial maps for each IC
  #   -if indicated, saved 4D nifti vol. w/ standard deviation
  #     of subj. spatial maps for each IC
  #   -saved matrix as 2D nifti vol. w/ mean of subj. time series
  #   -csv table w/ columns:
  #     Component               : filename & 4d-nifti index for component
  #     ICA Model Order         : ICA model order, if NA indicates hierarchy components
  #     Hierarchy Level         : Level of component in hierarchy, NA for non-hPCA/hICA comps.
  #     Hierarchy Variable Type : Options:
  #                               'initial sig.' : comp. is from lowest significant level of hierarchy
  #                               'sum' : comp. represents merged shared variance in hPCA (~PC1)
  #                               'diff': comp. represents difference of merged vars. in hPCA (~PC2)
  #                               'IC1' : higher-variance of new vars. in hICA, merged in higher levels
  #                               'IC2' : lower-variance of new vars. in hICA, stored
  #     Algorithm               : Hierarchical PCA or ICA algorithm
  #     Variable Index          : Index for comp. w/n ordered set of original variables,
  #                                 used for locating component within data
  #
  
  library(RNifti)
  
  if (verbose){ cat('\nCalculating group-level summary statistics (mean, variance) of ICA sources...')}
  
  ### Inputs & defaults ###
  subj_pca_suffix  = '_PCA1.RData'       # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData'    # filesuffix for group-level temporal PCA from group_PCA()
  group_ica_suffix = 'spatialICA.RData'  # filesuffix for group spatial ICA from group_spatialICA()
  br_suffix        = '_ICAsources.RData' # filesuffix for back-reconstructed subj. sources
  spatial_maps_suffix = '_ICAspatialMaps.nii' # filesuffix for 4D nii vol. w/ back-reconstructed spatial maps
  time_series_suffix  = '_ICAtimeseries.nii'  # filesuffix for 4D nii vol. w/ back-reconstructed t.s.
  mean_spatial_map_suffix = 'spatialICA_SpatialMapsMeans.nii' # filesuffix for mean of spatial maps
  var_spatial_map_suffix  = 'spatialICA_SpatialMapsVars.nii'  # filesuffix for variance of spatial maps
  mean_time_series_suffix = 'spatialICA_TimeSeriesMeans.nii'  # filesuffix for mean of time series

  prefix = as.character(prefix)
  create.csv = as.logical(create.csv)
  save_path = groupICA_path         # save to same dir.
  
  ### Prelims ###
  br_wc           = paste0('*', br_suffix)
  spatial_maps_wc = paste0('*', spatial_maps_suffix)
  time_series_wc  = paste0('*', time_series_suffix)
  if (!all(is.na(prefix))){
    br_wc           = paste0(prefix,'.', br_wc)
    spatial_maps_wc = paste0(prefix,'.', spatial_maps_wc)
    time_series_wc  = paste0(prefix,'.', time_series_wc)
  }
  
  spatial_map_files = list.files(groupICA_path, spatial_maps_wc, full.names=T)
  time_series_files = list.files(groupICA_path, time_series_wc, full.names=T)
  S = length(spatial_map_files)
  if (S==0 && verbose){  cat(paste0('\n...Warning:  Could not find any spatial map files in:  ',groupICA_path,'\n'))  }
  stopifnot(S > 0)
  
  ###################
  get_Time_minMax <- function(ts_files){
    # quick fn. to find limiting number of timepoints
    #   shared among subjs.
    S = length(ts_files)
    T.max = Inf
    for (s in 1:S){
      T.max = min(T.max, niftiHeader(ts_files[s])$dim[2])
    }
    return(T.max)
  } #################
  
  T.limit = get_Time_minMax(time_series_files)
  ts.nii.hdr = niftiHeader(time_series_files[1])
  ts.nii.hdr$dim[2] = T.limit
  
  
  ### Main fn. ###
  sm.mean.nii = readNifti(spatial_map_files[1])
  ts.mean.nii = readNifti(time_series_files[1])[1:T.limit,]
  for (s in 2:S){
    sm.nii = readNifti(spatial_map_files[s])
    stopifnot(all(dim(sm.nii) == dim(sm.mean.nii)))
    sm.mean.nii = sm.mean.nii + sm.nii
    
    ts.nii = readNifti(time_series_files[s])[1:T.limit,]
    stopifnot(all(dim(ts.nii) == dim(ts.mean.nii)))
    ts.mean.nii = ts.mean.nii + ts.nii
    
  }
  sm.mean.nii = sm.mean.nii / S
  ts.mean.nii = ts.mean.nii / S
  
  if (calc.var){
    sm.var.nii = sm.mean.nii * 0
    for (s in 1:S){
      sm.nii = readNifti(spatial_map_files[s])
      stopifnot(all(dim(sm.nii) == dim(sm.var.nii)))
      sm.var.nii = sm.var.nii + (sm.nii - sm.mean.nii)^2
    }
    sm.var.nii = sm.var.nii / (S-1)
  }
  ts.mean.nii = asNifti(ts.mean.nii, reference=ts.nii.hdr)
  
  
  ### Saving ###
  save.fname1 = mean_spatial_map_suffix
  save.fname2 = mean_time_series_suffix
  save.fname3 = var_spatial_map_suffix
  if (!all(is.na(prefix))){
    save.fname1 = paste0(prefix, save.fname1)
    save.fname2 = paste0(prefix, save.fname2)
    save.fname3 = paste0(prefix, save.fname3)
  }
  save.fname1 = file.path(save_path, save.fname1)
  save.fname2 = file.path(save_path, save.fname2)
  save.fname3 = file.path(save_path, save.fname3)

  if (verbose){  cat(paste0('\n...saving vol. of IC mean spatial maps as:  ', basename(save.fname1)))  }
  writeNifti(sm.mean.nii, save.fname1)
  if (verbose){  cat(paste0('\n...saving ~2D vol. of IC mean time series as:  ', basename(save.fname2)))  }
  writeNifti(ts.mean.nii, save.fname2)
  if (calc.var){
    if (verbose){  cat(paste0('\n...saving vol. of IC spatial map variances as:  ', basename(save.fname3)))  }
    writeNifti(sm.var.nii, save.fname3)
  }
  
  ### Create table ###
  if (create.csv){
    stopifnot(file.exists(save.fname1))
    
    Space.br = new.env()
    br_file = list.files(groupICA_path, br_wc, full.names=T)[1]
    load(br_file, envir=Space.br)
    stopifnot(all(c('br') %in% names(Space.br)))
    br = get('br', envir=Space.br)
    orig_var_inds = br$orig_var_inds
    ICA_model_order = br$ICA_model_order
    Hierarchy_level = br$Hierarchy_level
    hier_algorithm = br$hier_algorithm
    hier_ns_inds  = br$hier_ns_inds
    hier_k1_inds = br$hier_k1_inds
    hier_k2_inds = br$hier_k2_inds
    rm('br')

    IC_inds = paste0(basename(save.fname1))
    IC_inds = sub('\\.nii', '', IC_inds)
    K = dim(sm.mean.nii)[4]
    IC_inds = paste0(IC_inds,',',as.integer(1:K))
    
    ICA_model_order = get('ICA_model_order', envir=Space.br)
    if (all(is.na(ICA_model_order))){  ICA_model_order = rep(NA, K)  }
    stopifnot(length(IC_inds) == length(ICA_model_order))
    
    hier_algorithm = get('hier_algorithm', envir=Space.br)
    if (all(is.na(hier_algorithm))){  hier_algorithm = rep(NA, K)  }
    stopifnot(length(IC_inds) == length(hier_algorithm))
    
    Hierarchy_level = get('Hierarchy_level', envir=Space.br)
    if (all(is.na(Hierarchy_level))){  Hierarchy_level = rep(NA, K)  }
    stopifnot(length(IC_inds) == length(Hierarchy_level))
    
    varType = rep(NA,K)
    hier_k1_inds = get('hier_k1_inds', envir=Space.br)
    hier_k2_inds = get('hier_k2_inds', envir=Space.br)
    hier_ns_inds = get('hier_ns_inds', envir=Space.br)
    init0 = pc1 = pc2 = ic1 = ic2 = F
    hpca = which(hier_algorithm == 'hPCA')
    if (length(hpca) > 0){
      if (!all(is.na(hier_k1_inds))){  pc1 = hier_k1_inds[hier_k1_inds %in% hpca]  }
      if (!all(is.na(hier_k2_inds))){  pc2 = hier_k2_inds[hier_k2_inds %in% hpca]  }
    }
    hica = which(hier_algorithm == 'hICA')
    if (length(hica) > 0){
      if (!all(is.na(hier_k1_inds))){  ic1 = hier_k1_inds[hier_k1_inds %in% hica]  }
      if (!all(is.na(hier_k2_inds))){  ic2 = hier_k2_inds[hier_k2_inds %in% hica]  }
    }
    if (!all(is.na(hier_ns_inds))){
      init0 = hier_ns_inds[hier_ns_inds %in% c(hpca, hica)]
      varType[init0] = 'initial sig.'  # indicates vars. from lowest significant level of hierarchy
    }
    if (any(as.logical(pc1))){
      varType[pc1] = 'sum'  # shared variance of merged vars. at each level of hierarchy (~PC1 for hPCA)
    }
    if (any(as.logical(pc2))){
      varType[pc2] = 'diff' # differences between merged vars. at each levelof hierarchy (~PC2 for hPCA)
    }
    if (any(as.logical(ic1))){
      varType[ic1] = 'IC1'  # higher-var. new comp. at each level of hierarchy for hICA
    }
    if (any(as.logical(ic2))){
      varType[ic2] = 'IC2'  # lower-var. new comp. at each level of hierarchy for hICA
    }
    
    orig_var_inds = get('orig_var_inds', envir=Space.br)  # original/replacement index for var. 
    stopifnot(length(IC_inds) == length(orig_var_inds))
    
    header = c('Component', 'ICA Model Order', 'Hierarchy Level', 'Hierarchy Variable Type', 'Algorithm', 'Variable Index')
    csv.body = rbind(header,
                     cbind(IC_inds, ICA_model_order, Hierarchy_level, varType, hier_algorithm, orig_var_inds,
                           deparse.level=0), deparse.level=0)
    csv.fname = sub('\\.nii', '.csv', save.fname1)
    write.table(csv.body, file=csv.fname, 
                sep=",", qmethod="double",
                row.names=F, col.names=F)
    if (verbose){  cat(paste0('\n...saving table of spatial map info as:  ', basename(csv.fname)))  }
  }
  if (verbose){  cat(paste0('\n......saved in dir.:  ', save_path, '\n\n'))  }
  
} ##################################################################


