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
#' @param FOI_R0_values TBA
#' @param sorting_metric_values Vector of values of a metric associated with parameter sets used to generate FOI_R0_values
#'   (e.g. total number of deaths calculated across yellow fever-endemic countries when the relevant parameters are
#'   used to calculate burden), used to sort values for each region; if set to NULL, values are simply sorted from lowest
#'   to highest
#' '
#' @export
#'
get_FOI_R0_dist_data <- function(FOI_R0_values = list(), sorting_metric_values = NULL){

  regions=FOI_R0_values$regions
  n_regions=length(regions)
  n_entries=dim(FOI_R0_values$FOI)[1]
  if(is.null(sorting_metric_values)==FALSE){
    assert_that(length(sorting_metric_values)==n_entries)
    sort_order=order(sorting_metric_values)
  }
  if(length(dim(FOI_R0_values$FOI))==3){
    FOI_R0_values=list(regions=FOI_R0_values$regions,FOI=rowMeans(FOI_R0_values$FOI,dims=2),R0=rowMeans(FOI_R0_values$R0,dims=2))
  }

  FOI_R0_summary=data.frame(array(NA,dim=c(n_regions,16)))
  colnames(FOI_R0_summary)=c("n_region","region","FOI_025","FOI_25","FOI_50","FOI_75","FOI_975","FOI_mean","FOI_cv",
                             "R0_025","R0_25","R0_50","R0_75","R0_975","R0_mean","R0_cv")
  FOI_R0_summary$n_region=c(1:n_regions)
  FOI_R0_summary$region=regions
  n_025=ceiling(n_entries*0.025)
  n_25=ceiling(n_entries*0.25)
  n_75=max(1, floor(n_entries*0.75))
  n_975=max(1, floor(n_entries*0.975))

  for(i in 1:n_regions){
    if(is.null(sorting_metric_values)==FALSE){
      FOI_values=FOI_R0_values$FOI[sort_order,i]
      R0_values=FOI_R0_values$R0[sort_order,i]
    } else {
      FOI_values=sort(FOI_R0_values$FOI[,i])
      R0_values=sort(FOI_R0_values$R0[,i])
    }
    FOI_R0_summary$FOI_025[i]=FOI_values[n_025]
    FOI_R0_summary$FOI_25[i]=FOI_values[n_25]
    FOI_R0_summary$FOI_50[i]=median(FOI_values)
    FOI_R0_summary$FOI_75[i]=FOI_values[n_75]
    FOI_R0_summary$FOI_975[i]=FOI_values[n_975]
    FOI_R0_summary$FOI_mean[i]=mean(FOI_values)
    FOI_R0_summary$FOI_cv[i]=sqrt(var(FOI_values))/FOI_R0_summary$FOI_mean[i]
    FOI_R0_summary$R0_025[i]=R0_values[n_025]
    FOI_R0_summary$R0_25[i]=R0_values[n_25]
    FOI_R0_summary$R0_50[i]=median(R0_values)
    FOI_R0_summary$R0_75[i]=R0_values[n_75]
    FOI_R0_summary$R0_975[i]=R0_values[n_975]
    FOI_R0_summary$R0_mean[i]=mean(R0_values)
    FOI_R0_summary$R0_cv[i]=sqrt(var(R0_values))/FOI_R0_summary$R0_mean[i]
  }

  return(FOI_R0_summary)
}
#-------------------------------------------------------------------------------
# TODO - Add function or expand above function to calculate metric values as well as sorting parameter sets
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
mcmc_R0_value_probs <- function(FOI_R0_values=list(), R0_target_values=c(1.0)){

  regions=unique(FOI_R0_values$region)
  n_regions=length(regions)
  n_entries=length(FOI_R0_values$region)/n_regions

  n_values=length(R0_target_values)
  R0_probs=data.frame(region=regions)
  for(i in 1:n_values){
    R0_probs[, i+1]=rep(NA, n_regions)
    colnames(R0_probs)[i+1]=paste0("P(R0>=", R0_target_values[i], ")")
  }
  for(n_region in 1:n_regions){
    lines=n_region+((c(1:n_entries)-1)*n_regions)
    R0_values=FOI_R0_values$R0[lines]
    for(i in 1:n_values){
      R0_probs[n_region, i+1]=sum(R0_values>=R0_target_values[i])/n_entries
    }
  }

  return(R0_probs)
}
