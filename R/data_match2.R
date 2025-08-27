extra_param_names <- c("vaccine_efficacy","p_severe_inf","p_death_severe_inf","p_rep_severe","p_rep_death","m_FOI_BRA")
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
#' @param template TBA
#' @param ... = Additional parameters/flags/etc. (n_reps, mode_start, time_inc, enviro_data_const, enviro_data_var,
#'   vaccine_efficacy, p_rep_severe, p_rep_death, m_FOI_BRA, deterministic, mode_time,
#'   mode_parallel, cluster, p_severe_inf, p_death_severe_inf)
#'
#' @export
#'
data_match_single2 <- function(params = c(), input_data = list(), template = list(), ...){

  #assert_that(all(params>0), msg = "All parameter values must be positive")
  n_params=ncol(params)
  assert_that(input_data_check(input_data),
              msg = "Input data must be in standard format (see https://mrc-ide.github.io/YEP/articles/CGuideAInputs.html )")
  consts<-list(...)

  # Checks
  assert_that(is.logical(consts$deterministic))
  assert_that(consts$mode_start %in% c(0, 1, 3),
              msg = "mode_start must have value 0, 1 or 3 (NB 3 should be changed to 1)")
  if(is.null(template$case)==FALSE){
    assert_that(all(template$case$cases==round(template$case$cases,0)),msg="Case data values must be integers")
    assert_that(all(template$case$deaths==round(template$case$deaths,0)),msg="Case data values must be integers")
  }
  regions = regions_breakdown(c(template$sero$region,template$case$region))
  n_regions=length(regions)
  assert_that(all(regions %in% consts$enviro_data_const$region),
              msg = "Time-invariant environmental data must be available for all regions in observed data")
  if(is.null(consts$enviro_data_var)==FALSE){
    assert_that(enviro_data_var_check(consts$enviro_data_var))
    assert_that(all(regions %in% consts$enviro_data_var$regions),
                msg = "Time-variant environmental data must be available for all regions in observed data")
  }

  #Truncate input and environmental data to only include relevant regions
  input_data = input_data_truncate(input_data,regions)
  enviro_data_const = subset(enviro_data_const, enviro_data_const$region %in% regions)
  if(is.null(enviro_data_var)==FALSE){enviro_data_var = enviro_data_var_truncate(enviro_data_var,regions)}

  #Designate constant and variable covariates
  const_covars = colnames(enviro_data_const)[c(2:ncol(enviro_data_const))]
  var_covars = enviro_data_var$env_vars
  covar_names = c(const_covars,var_covars)
  n_env_vars = length(covar_names)
  n_extra_vars = n_params - (2*n_env_vars)

  i_FOI_const = c(1:n_env_vars)[covar_names %in% const_covars] + n_extra_vars
  i_FOI_var = c(1:n_env_vars)[covar_names %in% var_covars] + n_extra_vars
  i_R0_const = i_FOI_const + n_env_vars
  i_R0_var = i_FOI_var + n_env_vars

  #frac = 1.0/consts$n_reps
  #n_params = length(params)

  #Get additional values - TODO: Make flexible?
  vaccine_efficacy = p_severe_inf = p_death_severe_inf = p_rep_severe = p_rep_death = m_FOI_BRA = 1.0
  for(var_name in extra_param_names){
    if(is.numeric(consts[[var_name]]) == FALSE){
      i = match(var_name, names(params))
      assign(var_name, as.numeric(params[i]))
    } else {
      assign(var_name, consts[[var_name]])
    }
  }

  #Get FOI and R0 values
  FOI_values = R0_values = rep(0, n_regions)
  FOI_values = epi_param_calc(coeffs_const = exp(as.numeric(params[i_FOI_const])), coeffs_var = exp(as.numeric(params[i_FOI_var])),
                              enviro_data_const = consts$enviro_data_const,enviro_data_var = consts$enviro_data_var)
  for(n_region in 1:n_regions){ #Apply Brazil FOI multiplier to relevant regions
    if(substr(input_data$region_labels[n_region],1,3) == "BRA"){FOI_values[n_region] = FOI_values[n_region]*m_FOI_BRA}
  }
  R0_values = epi_param_calc(coeffs_const = exp(as.numeric(params[i_R0_const])), coeffs_var = exp(as.numeric(params[i_R0_var])),
                             enviro_data_const = consts$enviro_data_const,enviro_data_var = consts$enviro_data_var)


  #Generate modelled data over all regions
  dataset <- Generate_Dataset(FOI_values, R0_values, input_data, template, vaccine_efficacy,
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
#' @param template TBA
#' @param ... = Constant additional parameters/flags/etc. (n_reps, mode_start, time_inc, enviro_data_const, enviro_data_var,
#'   vaccine_efficacy, p_rep_severe, p_rep_death, m_FOI_BRA, deterministic, mode_time,
#'   mode_parallel, cluster, p_severe_inf, p_death_severe_inf)
#'
#' @export
#'
data_match_multi2 <- function(param_sets = list(), input_data = list(), template = list(), ...){

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
    model_data_all[[i]] <- data_match_single2(params, input_data, template, ...)
  }

  return(model_data_all)
}
