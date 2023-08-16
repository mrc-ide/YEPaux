# Functions relating to processing of model calculations for multiple parameter sets
#-------------------------------------------------------------------------------
#' @title get_FOI_R0_dist_data
#'
#' @description Obtain parameters of distribution of spillover FOI and/or R0 values
#'   obtained for multiple regions using get_mcmc_FOI_R0_data()
#'
#' @details Take output of get_mcmc_FOI_R0_data() (data frame of all spillover FOI and/or
#'   R0 values obtained from MCMC results for multiple regions) and obtain information
#'   on the distribution of values for each region, including mean, coefficient of variation,
#'   median and quartiles of 95% and 50% critical intervals (values below which 2.5%, 25%, 75% and
#'   97.5% of values lie).
#'
#' @param FOI_R0_values Data frame of spillover FOI and/or R0 values obtained for multiple
#'   regions using get_mcmc_FOI_R0_data()
#' '
#' @export
#'
get_FOI_R0_dist_data <- function(FOI_R0_values=list()){
  assert_that(is.data.frame(FOI_R0_values))
  assert_that(colnames(FOI_R0_values)[1]=="n_region")
  assert_that(colnames(FOI_R0_values)[2]=="region")
  assert_that(colnames(FOI_R0_values)[3]=="FOI")
  if(length(colnames(FOI_R0_values))==3){
    flag_R0=0
  } else {
    flag_R0=1
  }

  regions=unique(FOI_R0_values$region)
  n_regions=length(regions)
  assert_that(all(FOI_R0_values$region[c(1:n_regions)]==regions))
  n_entries=nrow(FOI_R0_values)/n_regions
  bl=rep(NA,n_regions)
  if(flag_R0==0){
    FOI_R0_summary=data.frame(n_region=c(1:n_regions),region=regions,
                              FOI_025=bl,FOI_25=bl,FOI_50=bl,FOI_75=bl,FOI_975=bl,FOI_mean=bl,FOI_cv=bl)
  } else {
    FOI_R0_summary=data.frame(n_region=c(1:n_regions),region=regions,
                              FOI_025=bl,FOI_25=bl,FOI_50=bl,FOI_75=bl,FOI_975=bl,FOI_mean=bl,FOI_cv=bl,
                              R0_025=bl,R0_25=bl,R0_50=bl,R0_75=bl,R0_975=bl,R0_mean=bl,R0_cv=bl)
  }
  n_025=ceiling(n_entries*0.025)
  n_25=ceiling(n_entries*0.25)
  n_75=floor(n_entries*0.75)
  n_975=floor(n_entries*0.975)

  for(i in 1:n_regions){
    lines=i+(n_regions*c(0:(n_entries-1)))
    FOI_values=sort(FOI_R0_values$FOI[lines])
    R0_values=sort(FOI_R0_values$R0[lines])
    FOI_R0_summary$FOI_025[i]=FOI_values[n_025]
    FOI_R0_summary$FOI_25[i]=FOI_values[n_25]
    FOI_R0_summary$FOI_50[i]=median(FOI_values)
    FOI_R0_summary$FOI_75[i]=FOI_values[n_75]
    FOI_R0_summary$FOI_975[i]=FOI_values[n_975]
    FOI_mean=mean(FOI_values)
    FOI_R0_summary$FOI_mean[i]=FOI_mean
    FOI_R0_summary$FOI_cv[i]=sqrt(var(FOI_values))/FOI_mean
    if(flag_R0==1){
      FOI_R0_summary$R0_025[i]=R0_values[n_025]
      FOI_R0_summary$R0_25[i]=R0_values[n_25]
      FOI_R0_summary$R0_50[i]=median(R0_values)
      FOI_R0_summary$R0_75[i]=R0_values[n_75]
      FOI_R0_summary$R0_975[i]=R0_values[n_975]
      R0_mean=mean(R0_values)
      FOI_R0_summary$R0_mean[i]=R0_mean
      FOI_R0_summary$R0_cv[i]=sqrt(var(R0_values))/R0_mean
    }
  }

  return(FOI_R0_summary)
}
#-------------------------------------------------------------------------------
#' @title mcmc_R0_value_probs
#'
#' @description Get sets of values by region of probability that R0 equals or exceeds one or more values
#'
#' @details [TBA]
#'
#' @param FOI_R0_values Data frame of spillover FOI and/or R0 values obtained for multiple
#'   regions using get_mcmc_FOI_R0_data()
#' @param R0_target_values Vector of target values of R0
#' '
#' @export
#'
mcmc_R0_value_probs <- function(FOI_R0_values=list(),R0_target_values=c(1.0)){

  regions=unique(FOI_R0_values$region)
  n_regions=length(regions)
  n_entries=length(FOI_R0_values$region)/n_regions

  n_values=length(R0_target_values)
  R0_probs=data.frame(region=regions)
  for(i in 1:n_values){
    R0_probs[,i+1]=rep(NA,n_regions)
    colnames(R0_probs)[i+1]=paste0("P(R0>",R0_target_values[i],")")
  }
  for(n_region in 1:n_regions){
    lines=n_region+((c(1:n_entries)-1)*n_regions)
    R0_values=FOI_R0_values$R0[lines]
    for(i in 1:n_values){
      R0_probs[n_region,i+1]=sum(R0_values>=R0_target_values[i])/n_entries
    }
  }

  return(R0_probs)
}
#-------------------------------------------------------------------------------
# TODO - Add function for doing global burden calculations (or similar) based on
# multiple parameter sets and finding which sets give median, quartile etc. output
#-------------------------------------------------------------------------------
#' @title mcmc_datasets_sort_by_burden
#'
#' @description [TBA]
#'
#' @details [TBA]
#'
#' @param input_frame Data frame of MCMC output data produced by get_mcmc_data(), truncate_mcmc_data() and/or
#'   combine_multichain()
#' '
#' @export
#'
mcmc_datasets_sort_by_burden <- function(input_frame = list()){


  return(NULL)
}
