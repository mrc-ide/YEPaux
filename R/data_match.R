#Functions for generating sets of modelled data to compare with observed data and displaying comparative graphs
#-------------------------------------------------------------------------------
#' @title data_match_single
#'
#' @description Function which runs the model to create simulated data corresponding to supplied observed data for a
#'   single set of parameters with one or more repetitions
#'
#' @details This function takes in observed data and produces a corresponding simulated dataset for comparison using
#'   supplied model parameters (in the params input variable), population and vaccination settings (in the
#'   input_data input variable) and other settings including parameter type, time increment and environmental covariate
#'   values (in the const_list input variable). The simulated dataset will be able to be compared directly to the
#'   observed dataset - individual data values will be produced for the same years and regions (or combinations of
#'   regions), seroprevalence data will be produced for the same age groups, etc.
#'
#' @param params Values of input parameters in order FOI/FOI coefficients, R0/R0 coefficients, [etc.]
#' @param input_data List of population and vaccination data for multiple regions, with tables to cross-reference
#' with observed data, added using input_data_process
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param ... = Additional parameters/flags/etc. (n_reps, mode_start, time_inc, enviro_data_const, enviro_data_var,
#'   vaccine_efficacy, p_rep_severe, p_rep_death, m_FOI_Brazil, deterministic, mode_time,
#'   mode_parallel, cluster, p_severe_inf, p_death_severe_inf)
#'
#' @export
#'
data_match_single <- function(params = c(), input_data = list(), obs_sero_data = NULL, obs_case_data = NULL, ...){

  #assert_that(all(params>0), msg = "All parameter values must be positive")
  assert_that(input_data_check(input_data),
              msg = "Input data must be in standard format (see https://mrc-ide.github.io/YEP/articles/CGuideAInputs.html )")
  const_list<-list(...)
  #TODO - Add additional assert_that checks?

  #Process input data to check that all regions with sero, case and/or outbreak data supplied are present, remove
  #regions without any supplied data, and add cross-referencing tables for use when calculating likelihood. Take
  #subset of environmental data (if used) and check that environmental data available for all regions
  input_data = input_data_process(input_data, obs_sero_data, obs_case_data)
  regions = names(table(input_data$region_labels)) #Regions in new processed input data list
  n_regions = length(regions)
  assert_that(all(regions %in% const_list$enviro_data_const$region),
              msg = "Time-invariant environmental data must be available for all regions in observed data")
  const_list$enviro_data_const = subset(const_list$enviro_data_const, const_list$enviro_data_const$region %in% regions)
  if(is.null(const_list$enviro_data_var)==FALSE){
    assert_that(enviro_data_var_check(const_list$enviro_data_var))
    assert_that(all(regions %in% enviro_data_var$regions),
                msg = "Time-variant environmental data must be available for all regions in observed data")
    enviro_data_var = enviro_data_var_truncate(enviro_data_var,regions)
  }


  frac = 1.0/const_list$n_reps
  n_params = length(params)
  extra_estimated_params = c()
  if(is.na(const_list$vaccine_efficacy) == TRUE){extra_estimated_params = append(extra_estimated_params, "vaccine_efficacy")}
  if(is.na(const_list$p_rep_severe) == TRUE){extra_estimated_params = append(extra_estimated_params, "p_rep_severe")}
  if(is.na(const_list$p_rep_death) == TRUE){extra_estimated_params = append(extra_estimated_params, "p_rep_death")}
  if(is.na(const_list$m_FOI_Brazil) == TRUE){extra_estimated_params = append(extra_estimated_params, "m_FOI_Brazil")}
  names(params) = create_param_labels(const_list$enviro_data_const, const_list$enviro_data_var, extra_estimated_params)

  checks <- mcmc_checks(params, n_regions, list(type = "zero"), const_list$enviro_data_const, const_list$enviro_data_var,
                        add_values = list(vaccine_efficacy = const_list$vaccine_efficacy, p_rep_severe = const_list$p_rep_severe,
                                          p_rep_death = const_list$p_rep_death, m_FOI_Brazil = const_list$m_FOI_Brazil),
                        extra_estimated_params)

  #Get additional values - TODO: Make flexible?
  vaccine_efficacy = p_rep_severe = p_rep_death = m_FOI_Brazil = 1.0
  for(var_name in c("vaccine_efficacy", "p_rep_severe", "p_rep_death", "m_FOI_Brazil")){
    if(is.numeric(const_list[[var_name]]) == FALSE){
      i = match(var_name, names(params))
      assign(var_name, as.numeric(params[i]))
    } else {
      assign(var_name, const_list[[var_name]])
    }
  }

  #Designate constant and variable covariates
  const_covars = colnames(const_list$enviro_data_const)[c(2:ncol(const_list$enviro_data_const))]
  var_covars = const_list$enviro_data_var$env_vars
  covar_names = c(const_covars,var_covars)
  n_env_vars = length(covar_names)
  i_FOI_const = c(1:n_env_vars)[covar_names %in% const_covars]
  i_FOI_var = c(1:n_env_vars)[covar_names %in% var_covars]
  i_R0_const = i_FOI_const + n_env_vars
  i_R0_var = i_FOI_var + n_env_vars

  #Get FOI and R0 values
  FOI_values = R0_values = rep(0, n_regions)
  FOI_values = epi_param_calc(coeffs_const = as.numeric(params[i_FOI_const]), coeffs_var = params[i_FOI_var],
                              enviro_data_const = const_list$enviro_data_const,enviro_data_var = const_list$enviro_data_var)
  R0_values = epi_param_calc(coeffs_const = as.numeric(params[i_R0_const]), coeffs_var = params[i_R0_var],
                             enviro_data_const = const_list$enviro_data_const,enviro_data_var = const_list$enviro_data_var)

  #Generate modelled data over all regions
  dataset <- Generate_Dataset(FOI_values, R0_values, input_data, obs_sero_data, obs_case_data, vaccine_efficacy,
                              const_list$time_inc, const_list$mode_start, const_list$start_SEIRV, const_list$mode_time,
                              const_list$n_reps, const_list$deterministic, const_list$p_severe_inf, const_list$p_death_severe_inf,
                              p_rep_severe, p_rep_death, const_list$mode_parallel, const_list$cluster, output_frame = FALSE)

  return(dataset)
}
#-------------------------------------------------------------------------------
#' @title data_match_multi
#'
#' @description Function which runs the model to create simulated data corresponding to supplied observed data for
#'   multiple parameter sets
#'
#' @details This function runs the data_match_single() function for multiple parameter sets. It takes in observed data
#'   and produces corresponding simulated datasets for comparison using multiple sets of supplied model parameters (in
#'   the param_sets input variable), population and vaccination settings (in the input_data input variable) and other
#'   settings including parameter type, time increment and environmental covariate values (in the const_list input
#'   variable). The simulated dataset will be able to be compared directly to the observed dataset - individual data
#'   values will be produced for the same years and regions (or combinations of regions), seroprevalence data will be
#'   produced for the same age groups, etc.
#'
#' @param param_sets Data frame of log values of proposed parameters, one set per row
#' @param input_data List of population and vaccination data for multiple regions, with tables to cross-reference
#'   with observed data, added using input_data_process2
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param ... = Constant additional parameters/flags/etc. (n_reps, mode_start, time_inc, enviro_data_const, enviro_data_var,
#'   vaccine_efficacy, p_rep_severe, p_rep_death, m_FOI_Brazil, deterministic, mode_time,
#'   mode_parallel, cluster, p_severe_inf, p_death_severe_inf)
#'
#' @export
#'
data_match_multi <- function(param_sets = list(), input_data = list(), obs_sero_data = NULL,
                             obs_case_data = NULL, ...){
  assert_that(is.data.frame(param_sets), msg = "param_sets must be a data frame")

  n_param_sets = nrow(param_sets)
  model_data_all = list()
  cat("\nSet:\n")
  for(i in 1:n_param_sets){
    cat("\t", i)
    params = param_sets[i, ]
    model_data_all[[i]] <- data_match_single(params, input_data, obs_sero_data, obs_case_data, ...)
  }

  return(model_data_all)
}
#-------------------------------------------------------------------------------
#' @title sero_match_graphs
#'
#' @description Function to create a series of graphs comparing modelled and observed serological data from results
#'   generated from data_match_multi() function
#'
#' @details Takes in simulated datasets produced by data_match_multi() function along with observed seroprevalence
#'   dataset used to generate the simulated datasets, and produces graph(s) (one graph per survey) showing the observed
#'   and simulated data on the same set(s) of axes, with observed data displayed as points and lines and simulated data
#'   as coloured bands.
#'
#' @param model_data Simulated datasets produced by data_match_multi() function
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param plot_type Form in which to plot model data: "all": bands showing 95\% and 50\% of all values;
#'   "mean": bands showing 95\% and 50\% confidence intervals for mean of all values
#' @param text_size Size of text to display on graphs
#' @param hide_observed If TRUE, indicates that supplied observed data is "dummy" data only supplied to indicate years
#'   and age ranges, and should not be plotted on the graph(s)
#'
#' @export
#'
sero_match_graphs <- function(model_data = list(), obs_sero_data = list(), plot_type = "mean", text_size = 1.0,
                              hide_observed = FALSE){

  #TODO - Add image adjustment options

  assert_that(is.list(model_data)) #TODO - improve checks on model_data
  assert_that(is.data.frame(obs_sero_data))
  assert_that(plot_type %in% c("mean", "all"), msg = "plot_type must be 'mean' or 'all'")
  assert_that(is.numeric(text_size))

  if(typeof(model_data[[1]]) == "list"){
    data_type = "multi"
    n_param_sets = length(model_data)
  } else{
    assert_that(typeof(model_data[[1]]) == "double") #TODO - improve checks on model_data
    data_type = "single"
    n_param_sets = 1
  }

  obs_sero_values = obs_sero_data$positives/obs_sero_data$samples
  n_sero_values = length(obs_sero_values)
  obs_sero_values[is.nan(obs_sero_values)] = 0.0
  obs_sero_values_low = obs_sero_values_high = rep(NA, n_sero_values)
  for(i in 1:n_sero_values){
    if(obs_sero_data$samples[i]>0){
      CI = prop.test(x = obs_sero_data$positives[i], n = obs_sero_data$samples[i])
      obs_sero_values_low[i] = CI$conf.int[1]
      obs_sero_values_high[i] = CI$conf.int[2]
    } else {
      obs_sero_values_low[i] = 0
      obs_sero_values_high[i] = 0
    }
  }

  model_sero_values = array(NA, dim = c(length(obs_sero_values), n_param_sets))
  lines_all = graph_lines = c(1:length(obs_sero_values))
  model_CI50_low = model_CI50_high = model_CI95_low = model_CI95_high = rep(0, length(obs_sero_values))

  if(data_type == "multi"){
    for(i in 1:n_param_sets){model_sero_values[, i] = model_data[[i]]$model_sero_values}
  }
  if(data_type == "single"){ model_sero_values[, 1] = model_data$model_sero_values }

  if(data_type == "multi"){
    if(plot_type == "mean"){
      for(i in lines_all){
        CI_095 = CI(model_sero_values[i, ], ci = 0.95)
        model_CI95_low[i] = CI_095[3][[1]]
        model_CI95_high[i] = CI_095[1][[1]]
        CI_050 = CI(model_sero_values[i, ], ci = 0.50)
        model_CI50_low[i] = CI_050[3][[1]]
        model_CI50_high[i] = CI_050[1][[1]]
      }
    } else {
      n_095_low = ceiling(n_param_sets*0.025)
      n_095_high = max(1, floor(n_param_sets*0.975))
      n_050_low = ceiling(n_param_sets*0.25)
      n_050_high = max(1, floor(n_param_sets*0.75))
      for(i in lines_all){
        model_values_sorted = sort(model_sero_values[i, ])
        model_CI95_low[i] = model_values_sorted[n_095_low]
        model_CI95_high[i] = model_values_sorted[n_095_high]
        model_CI50_low[i] = model_values_sorted[n_050_low]
        model_CI50_high[i] = model_values_sorted[n_050_high]
      }
    }
  }
  if(data_type == "single"){
    dS = min(model_sero_values)*0.01
    model_CI95_low = model_sero_values-dS
    model_CI95_high = model_sero_values+dS
  }

  if(is.null(obs_sero_data$country_zone) == FALSE){
    data_regions = names(table(obs_sero_data$country_zone))
    graph_titles = c()
    n_graphs = 0
    for(region in data_regions){
      lines = lines_all[obs_sero_data$country_zone == region]
      subset = obs_sero_data[lines, ]
      years = as.numeric(names(table(subset$year)))
      for(year in years){
        lines2 = lines[subset$year == year]
        subset2 = subset(subset, year == year)
        n_graphs = n_graphs+1
        graph_lines[lines2] = n_graphs
        graph_titles = append(graph_titles, paste(region, year, sep = " "))
      }
    }
  } else {
    data_regions = names(table(obs_sero_data$region))
    graph_titles = c()
    n_graphs = 0
    for(region in data_regions){
      lines = lines_all[obs_sero_data$region == region]
      subset = subset(obs_sero_data, obs_sero_data$region == region)
      years = names(table(subset$year))
      for(year in years){
        lines2 = lines[subset$year == year]
        subset2 = subset(subset, year == year)
        n_graphs = n_graphs+1
        graph_lines[lines2] = n_graphs
        graph_titles = append(graph_titles, paste(region, year, sep = " "))
      }
    }
  }

  age_values = sero_obs = sero_obs_low = sero_obs_high = sero_model_low95 = NULL
  sero_model_low50 = sero_model_high95 = sero_model_high50 = NULL
  sero_graphs = list()
  for(i in 1:n_graphs){
    lines = lines_all[graph_lines == i]

    df = data.frame(age_values = obs_sero_data$age_min[lines], sero_obs = obs_sero_values[lines],
                  samples = obs_sero_data$samples[lines], sizes = obs_sero_data$samples[lines]/max(obs_sero_data$samples[lines]),
                  sero_obs_low = obs_sero_values_low[lines], sero_obs_high = obs_sero_values_high[lines],
                  sero_model_low95 = model_CI95_low[lines], sero_model_high95 = model_CI95_high[lines],
                  sero_model_low50 = model_CI50_low[lines], sero_model_high50 = model_CI50_high[lines])
    df$samples[df$samples == 0] = 1

    samples = sizes = NULL
    sero_graphs[[i]] <- ggplot(data = df) + theme_bw()+labs(title = graph_titles[i])
    sero_graphs[[i]] <- sero_graphs[[i]]+geom_ribbon(data = df, aes(x = age_values, ymin = sero_model_low95,
                                                                 ymax = sero_model_high95), fill = "blue", alpha = 0.5)
    sero_graphs[[i]] <- sero_graphs[[i]]+geom_line(data = df, aes(x = age_values, y = sero_model_low95), col = "blue")
    sero_graphs[[i]] <- sero_graphs[[i]]+geom_line(data = df, aes(x = age_values, y = sero_model_high95), col = "blue")
    if(data_type == "multi"){
      sero_graphs[[i]] <- sero_graphs[[i]]+geom_ribbon(data = df, aes(x = age_values, ymin = sero_model_low50,
                                                                   ymax = sero_model_high50), fill = "green", alpha = 0.5)
    }

    if(hide_observed == FALSE){
      sero_graphs[[i]] <- sero_graphs[[i]]+geom_point(data = df, aes(x = age_values, y = sero_obs, size = sizes), #TODO - amend size
                                                      show.legend = FALSE)
      sero_graphs[[i]] <- sero_graphs[[i]]+geom_errorbar(data = df, aes(x = age_values, ymin = sero_obs_low,
                                                                     ymax = sero_obs_high), width = 1.0)
    }

    sero_graphs[[i]] <- sero_graphs[[i]]+scale_x_continuous(name = "Age (min)", breaks = df$age_values, labels = df$age_values)
    sero_graphs[[i]] <- sero_graphs[[i]]+scale_y_continuous(name = "Seroprevalence")
    sero_graphs[[i]] <- sero_graphs[[i]]+theme(axis.text.x = element_text(size = text_size),
                                               axis.text.y = element_text(size = text_size),
                                               title = element_text(size = text_size))
  }

  return(sero_graphs)
}

#-------------------------------------------------------------------------------
#' @title case_match_graphs
#'
#' @description Function to create a series of graphs comparing modelled and observed case data from results
#'   generated from data_match_multi() function
#'
#' @details Takes in simulated datasets produced by data_match_multi() function along with observed annual case
#'   dataset used to generate the simulated datasets, and produces graph(s) (one graph per region or group of regions
#'   for which annual case data was supplied) showing the observed and simulated data on the same set(s) of axes, with
#'   observed data displayed as points and lines and simulated data as coloured bands.
#'
#' @param model_data Simulated datasets produced by data_match_multi() function
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param input_data List of population and vaccination data for multiple regions, with tables to cross-reference
#' with observed data, added using input_data_process
#' @param plot_type Form in which to plot model data: "all": bands showing 95\% and 50\% of all values; "mean": bands
#'   showing 95\% and 50\% confidence intervals for mean of all values
#' @param text_size Size of text to display on graphs
#' @param hide_observed If TRUE, indicates that supplied observed data is "dummy" data only supplied to indicate years
#'   and age ranges, and should not be plotted on the graph(s)
#'
#' @export
#'
case_match_graphs <- function(model_data = list(), obs_case_data = list(), input_data = list(), plot_type = "mean", text_size = 1.0,
                              hide_observed = FALSE){

  assert_that(is.list(model_data)) #TODO - improve checks on model_data
  assert_that(is.data.frame(obs_case_data))
  assert_that(plot_type %in% c("mean", "all"), msg = "plot_type must be 'mean' or 'all'")
  assert_that(is.numeric(text_size))

  if(typeof(model_data[[1]]) == "list"){
    data_type = "multi"
    n_param_sets = length(model_data)
  } else
  {
    assert_that(typeof(model_data[[1]]) == "double") #TODO - improve checks on model_data
    data_type = "single"
    n_param_sets = 1
  }

  obs_case_values = obs_case_data$cases
  obs_death_values = obs_case_data$deaths
  n_case_values = length(obs_case_values)
  n_year_values = match(obs_case_data$year, input_data$years_labels)
  obs_case_values_low = obs_case_values_high = obs_death_values_low = obs_death_values_high = rep(NA, n_case_values)
  for(i in 1:n_case_values){
    regions = strsplit(obs_case_data$region[i], ",")[[1]]
    n_region_values = input_data$region_labels %in% regions
    population = round(sum(input_data$pop_data[n_region_values, n_year_values[i], ]), digits = 0)
    CI = prop.test(x = obs_case_data$cases[i], n = population)
    obs_case_values_low[i] = round(CI$conf.int[1]*population, digits = 0)
    obs_case_values_high[i] = round(CI$conf.int[2]*population, digits = 0)
    CI = prop.test(x = obs_case_data$deaths[i], n = population)
    obs_death_values_low[i] = round(CI$conf.int[1]*population, digits = 0)
    obs_death_values_high[i] = round(CI$conf.int[2]*population, digits = 0)
  }

  model_case_values = model_death_values = array(NA, dim = c(length(obs_case_values), n_param_sets))
  if(data_type == "multi"){
    for(i in 1:n_param_sets){
      model_case_values[, i] = model_data[[i]]$model_case_values
      model_death_values[, i] = model_data[[i]]$model_death_values
    }
  }
  if(data_type == "single"){
    model_case_values[, 1] = model_data$model_case_values
    model_death_values[, 1] = model_data$model_death_values
  }

  model_cases_CI95_low = model_cases_CI95_high = model_cases_CI50_low = model_cases_CI50_high = rep(0, n_case_values)
  model_deaths_CI95_low = model_deaths_CI95_high = model_deaths_CI50_low = model_deaths_CI50_high = rep(0, n_case_values)
  if(data_type == "multi"){
    if(plot_type == "mean"){
      for(i in 1:n_case_values){
        CI_095 = CI(model_case_values[i, ], ci = 0.95)
        model_cases_CI95_low[i] = CI_095[3][[1]]
        model_cases_CI95_high[i] = CI_095[1][[1]]
        CI_050 = CI(model_case_values[i, ], ci = 0.50)
        model_cases_CI50_low[i] = CI_050[3][[1]]
        model_cases_CI50_high[i] = CI_050[1][[1]]
        CI_095 = CI(model_death_values[i, ], ci = 0.95)
        model_deaths_CI95_low[i] = CI_095[3][[1]]
        model_deaths_CI95_high[i] = CI_095[1][[1]]
        CI_050 = CI(model_death_values[i, ], ci = 0.50)
        model_deaths_CI50_low[i] = CI_050[3][[1]]
        model_deaths_CI50_high[i] = CI_050[1][[1]]
      }
    } else {
      n_095_low = ceiling(n_param_sets*0.025)
      n_095_high = max(1, floor(n_param_sets*0.975))
      n_050_low = ceiling(n_param_sets*0.25)
      n_050_high = max(1, floor(n_param_sets*0.75))
      for(i in 1:n_case_values){
        model_case_values_sorted = sort(model_case_values[i, ])
        model_cases_CI95_low[i] = model_case_values_sorted[n_095_low]
        model_cases_CI95_high[i] = model_case_values_sorted[n_095_high]
        model_cases_CI50_low[i] = model_case_values_sorted[n_050_low]
        model_cases_CI50_high[i] = model_case_values_sorted[n_050_high]
        model_death_values_sorted = sort(model_death_values[i, ])
        model_deaths_CI95_low[i] = model_death_values_sorted[n_095_low]
        model_deaths_CI95_high[i] = model_death_values_sorted[n_095_high]
        model_deaths_CI50_low[i] = model_death_values_sorted[n_050_low]
        model_deaths_CI50_high[i] = model_death_values_sorted[n_050_high]
      }
    }
  }
  if(data_type == "single"){
    model_cases_CI95_low = model_case_values-1
    model_cases_CI95_low[model_cases_CI95_low<0] = 0
    model_cases_CI95_high = model_case_values+1
    model_deaths_CI95_low = model_death_values-1
    model_deaths_CI95_low[model_deaths_CI95_low<0] = 0
    model_deaths_CI95_high = model_death_values+1
  }

  #TODO - Optimize use of ggplot
  data_regions = names(table(obs_case_data$region))
  n_graphs = length(data_regions)
  cases_graphs = deaths_graphs = list()
  years = case_obs = case_obs_low = case_obs_high = case_model_low95 = case_model_low50 = case_model_high95 = case_model_high50 = NULL
  death_obs = death_obs_low = death_obs_high = death_model_low95 = death_model_low50 = death_model_high95 = death_model_high50 = NULL
  for(i in 1:n_graphs){
    region = data_regions[i]
    lines = obs_case_data$region == region

    df = data.frame(years = obs_case_data$year[lines], case_obs = obs_case_values[lines], death_obs = obs_death_values[lines],
                  case_obs_low = obs_case_values_low[lines], case_obs_high = obs_case_values_high[lines],
                  case_model_low95 = model_cases_CI95_low[lines], case_model_high95 = model_cases_CI95_high[lines],
                  case_model_low50 = model_cases_CI50_low[lines], case_model_high50 = model_cases_CI50_high[lines],
                  death_obs_low = obs_death_values_low[lines], death_obs_high = obs_death_values_high[lines],
                  death_model_low95 = model_deaths_CI95_low[lines], death_model_high95 = model_deaths_CI95_high[lines],
                  death_model_low50 = model_deaths_CI50_low[lines], death_model_high50 = model_deaths_CI50_high[lines])

    cases_graphs[[i]] <- ggplot(data = df) + theme_bw()+labs(title = region)
    cases_graphs[[i]] <- cases_graphs[[i]]+geom_ribbon(data = df, aes(x = years, ymin = case_model_low95,
                                                                   ymax = case_model_high95), fill = "red", alpha = 1.0)
    cases_graphs[[i]] <- cases_graphs[[i]]+geom_line(data = df, aes(x = years, y = case_model_low95), col = "red")
    cases_graphs[[i]] <- cases_graphs[[i]]+geom_line(data = df, aes(x = years, y = case_model_high95), col = "red")
    if(data_type == "multi"){
      cases_graphs[[i]] <- cases_graphs[[i]]+geom_ribbon(data = df, aes(x = years, ymin = case_model_low50,
                                                                     ymax = case_model_high50), fill = "orange", alpha = 1.0)
      # cases_graphs[[i]] <- cases_graphs[[i]]+geom_line(data = df, aes(x = years, y = case_model_low50), col = "orange")
      # cases_graphs[[i]] <- cases_graphs[[i]]+geom_line(data = df, aes(x = years, y = case_model_high50), col = "orange")
    }
    if(hide_observed == FALSE){
      cases_graphs[[i]] <- cases_graphs[[i]]+geom_point(data = df, aes(x = years, y = case_obs))
      cases_graphs[[i]] <- cases_graphs[[i]]+geom_errorbar(data = df, aes(x = years, ymin = case_obs_low, ymax = case_obs_high),
                                                           width = 0.5)
    }
    cases_graphs[[i]] <- cases_graphs[[i]]+scale_x_continuous(name = "", breaks = c(min(df$years):max(df$years)),
                                                              labels = c(min(df$years):max(df$years)))
    cases_graphs[[i]] <- cases_graphs[[i]]+scale_y_continuous(name = "Cases")
    cases_graphs[[i]] <- cases_graphs[[i]]+theme(axis.text.x = element_text(size = text_size),
                                                 axis.text.y = element_text(size = text_size),
                                                 title = element_text(size = text_size))


    deaths_graphs[[i]] <- ggplot(data = df) + theme_bw()+labs(title = region)
    deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_ribbon(data = df, aes(x = years, ymin = death_model_low95,
                                                                     ymax = death_model_high95), fill = "red", alpha = 1.0)
    deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_line(data = df, aes(x = years, y = death_model_low95), col = "red")
    deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_line(data = df, aes(x = years, y = death_model_high95), col = "red")
    if(data_type == "multi"){
      deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_ribbon(data = df, aes(x = years, ymin = death_model_low50,
                                                                       ymax = death_model_high50), fill = "orange", alpha = 1.0)
      # deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_line(data = df, aes(x = years, y = death_model_low50), col = "orange")
      # deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_line(data = df, aes(x = years, y = death_model_high50), col = "orange")
    }
    if(hide_observed == FALSE){
      deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_point(data = df, aes(x = years, y = death_obs))
      deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_errorbar(data = df, aes(x = years, ymin = death_obs_low,
                                                                         ymax = death_obs_high), width = 0.5)
    }
    deaths_graphs[[i]] <- deaths_graphs[[i]]+scale_x_continuous(name = "", breaks = c(min(df$years):max(df$years)),
                                                                labels = c(min(df$years):max(df$years)))
    deaths_graphs[[i]] <- deaths_graphs[[i]]+scale_y_continuous(name = "Deaths")
    deaths_graphs[[i]] <- deaths_graphs[[i]]+theme(axis.text.x = element_text(size = text_size),
                                                   axis.text.y = element_text(size = text_size),
                                                   title = element_text(size = text_size))
  }

  return(list(cases_graphs = cases_graphs, deaths_graphs = deaths_graphs))
}
