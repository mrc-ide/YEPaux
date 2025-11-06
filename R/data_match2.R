extra_param_names = c("vaccine_efficacy","p_severe_inf","p_death_severe_inf",
                      "p_rep_severe","p_rep_death")
#Functions for generating sets of modelled data to compare with observed data and displaying comparative graphs
#-------------------------------------------------------------------------------
#' @title data_match_single2
#'
#' @description Function which runs the model to create simulated data corresponding to supplied observed data for a
#'   single set of parameters with one or more repetitions
#'
#' @details This function takes in observed data and produces a corresponding simulated dataset for comparison using
#'   supplied model parameters (in the params input variable), population and vaccination settings (in the
#'   input_data input variable) and other settings including parameter type, time increment and environmental covariate
#'   values (in the consts input variable). The simulated dataset will be able to be compared directly to the
#'   observed dataset - individual data values will be produced for the same years and regions (or combinations of
#'   regions), seroprevalence data will be produced for the same age groups, etc.
#'
#' @param params Values of input parameters in order FOI/FOI coefficients, R0/R0 coefficients, [etc.]
#' @param input_data List of population and vaccination data for multiple regions in standard format [TBA]
#' @param env_covar_values TBA
#' @param template TBA
#' @param fixed_extra TBA
#' @param ... = Additional parameters/flags/etc. (n_reps, mode_start, time_inc, enviro_data_const, enviro_data_var,
#'   vaccine_efficacy, p_rep_severe, p_rep_death, m_FOI_BRA, deterministic, mode_time,
#'   mode_parallel, cluster, p_severe_inf, p_death_severe_inf)
#'
#' @export
#'
data_match_single2 <- function(params = c(), input_data = list(), env_covar_values = list(),
                               template = list(), fixed_extra = list(), ...){

  #assert_that(all(params>0), msg = "All parameter values must be positive")
  n_params=ncol(params)
  assert_that(input_data_check(input_data),
              msg = "Input data must be in standard format (see https://mrc-ide.github.io/YEP/articles/CGuideAInputs.html )")
  consts<-list(...) #TODO - account for potential missing parameters

  # Checks
  assert_that(is.logical(consts$deterministic))
  assert_that(consts$mode_start %in% c(0, 1, 2),msg = "mode_start must have value 0, 1 or 2")
  if(is.null(template$case)==FALSE){
    assert_that(all(template$case$cases==round(template$case$cases,0)),msg="Case data values must be integers")
    assert_that(all(template$case$deaths==round(template$case$deaths,0)),msg="Case data values must be integers")
  }
  regions = regions_breakdown(c(template$sero$region,template$case$region))
  n_regions=length(regions)
  n_env_vars=dim(env_covar_values)[1]
  n_coeffs=2*n_env_vars
  n_extra=n_params-n_coeffs
  assert_that(dim(env_covar_values)[2]==n_regions)
  assert_that(length(input_data$region_labels)==n_regions)

  #Get additional values - TODO: Make flexible?
  vaccine_efficacy = p_severe_inf = p_death_severe_inf = p_rep_severe = p_rep_death = 1.0

  for(var_name in extra_param_names){
    if(var_name %in% names(params)){
      i = match(var_name, names(params))
      assign(var_name, as.numeric(params[i]))
    } else {
      assert_that(var_name %in% names(fixed_extra))
      assign(var_name, fixed_extra[[var_name]])
    }
  }

  #Get FOI and R0 values
  #TODO - get coeff indices from param names?
  consts$n_r=n_regions
  consts$ref_BRA=which(substr(input_data$region_labels,1,3)=="BRA")
  log_FOI_coeffs=params[c(1:n_env_vars)+n_extra]
  log_R0_coeffs=params[c(1:n_env_vars)+n_extra+n_env_vars]
  epi_params = epi_param_calc2(pars_fixed=consts, env_covar_values, log_FOI_coeffs,
                               log_R0_coeffs,vars_extra=params)

  #Generate modelled data over all regions
  dataset <- Generate_Dataset(FOI_values = epi_params$FOI_spillover, R0_values = epi_params$R0,
                              input_data, template, vaccine_efficacy,
                              consts$time_inc, consts$mode_start, consts$start_SEIRV, consts$mode_time,
                              consts$n_reps, consts$deterministic, p_severe_inf, p_death_severe_inf,
                              p_rep_severe, p_rep_death, consts$mode_parallel, consts$cluster, output_frame = FALSE,
                              consts$seed, template$region_grouping)

  return(dataset)
}
#-------------------------------------------------------------------------------
#' @title data_match_multi2
#'
#' @description Function which runs the model to create simulated data corresponding to supplied observed data for
#'   multiple parameter sets
#'
#' @details This function runs the data_match_single2() function for multiple parameter sets. It takes in observed data
#'   and produces corresponding simulated datasets for comparison using multiple sets of supplied model parameters (in
#'   the param_sets input variable), population and vaccination settings (in the input_data input variable) and other
#'   settings including parameter type, time increment and environmental covariate values (in the consts input
#'   variable). The simulated dataset will be able to be compared directly to the observed dataset - individual data
#'   values will be produced for the same years and regions (or combinations of regions), seroprevalence data will be
#'   produced for the same age groups, etc.
#'
#' @param param_sets Data frame of values of proposed parameters, one set per row, with parameter names as headings
#' @param input_data List of population and vaccination data for multiple regions, with tables to cross-reference
#'   with observed data, added using input_data_process2
#' @param env_covar_values TBA
#' @param template TBA
#' @param fixed_extra TBA
#' @param ... = Constant additional parameters/flags/etc. (n_reps, mode_start, time_inc, enviro_data_const, enviro_data_var,
#'   vaccine_efficacy, p_rep_severe, p_rep_death, m_FOI_BRA, deterministic, mode_time,
#'   mode_parallel, cluster, p_severe_inf, p_death_severe_inf)
#'
#' @export
#'
data_match_multi2 <- function(param_sets = list(), input_data = list(), env_covar_values = list(),
                              template = list(), fixed_extra = list(), ...){

  #TODO - add assert_that functions?
  assert_that(is.data.frame(param_sets), msg = "param_sets must be a data frame")

  if(is.null(template$xref_sero)){
    template$xref_sero = template_region_xref(template$sero,input_data$region_labels)
  }
  if(is.null(template$xref_case)){
    template$xref_case = template_region_xref(template$case,input_data$region_labels)
  }
  template$region_grouping = get_region_grouping(input_data$region_labels,template,mode_grouping=2)

  n_param_sets = nrow(param_sets)
  model_data_all = list()
  cat("\nSet:\n")
  for(i in 1:n_param_sets){
    cat("\t", i)
    params = param_sets[i, ]
    model_data_all[[i]] <- data_match_single2(params, input_data, env_covar_values, template, fixed_extra, ...)
  }

  return(model_data_all)
}
