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
  assert_that(colnames(FOI_R0_values)[2]=="FOI")
  if(length(colnames(FOI_R0_values))==2){
    flag_R0=0
  } else {
    flag_R0=1
  }

  n_regions=length(unique(FOI_R0_values$n_region))
  n_entries=nrow(FOI_R0_values)/n_regions
  bl=rep(NA,n_regions)
  FOI_R0_summary=data.frame(n_region=c(1:n_regions),FOI_025=bl,FOI_25=bl,FOI_50=bl,FOI_75=bl,FOI_975=bl,FOI_mean=bl,FOI_cv=bl)
  if(flag_R0==1){
    FOI_R0_summary$R0_025=FOI_R0_summary$R0_25=FOI_R0_summary$R0_50=FOI_R0_summary$R0_75=
      FOI_R0_summary$R0_975=FOI_R0_summary$R0_mean=FOI_R0_summary$R0_cv=bl
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
