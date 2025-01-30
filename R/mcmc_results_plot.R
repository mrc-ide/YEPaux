#-------------------------------------------------------------------------------
#' @title plot_mcmc_FOI_R0_data
#'
#' @description Plot spillover force of infection (FOI) and reproduction number (R0) values from MCMC output data
#'
#' @details This function takes in spillover force of infection (FOI) and reproduction number (R0)
#' values calculated from MCMC output data using the get_mcmc_FOI_R0_data() function and plots the values in one of a
#' number of different ways to display the spread of values obtained (box plots, violin plots, simple plots with error
#' bars, and scatter plots of FOI vs R0).
#'
#' @param FOI_R0_values TBA
#' @param plot_type Type of plots to create (choose from "box", "violin", "error_bars", "scatter")
#' @param text_size1 Text size parameter for axis labels
#' '
#' @export
#'
plot_mcmc_FOI_R0_data <- function(FOI_R0_values=list(), plot_type="box", text_size1=10.0){

  assert_that(is.numeric(text_size1))
  if(is.null(data_frame$R0)==TRUE){
    assert_that(plot_type %in% c("box", "violin", "error_bars"))
  } else {
    assert_that(plot_type %in% c("box", "violin", "error_bars", "scatter"))
  }

  #Convert data to single time values if needed
  if(length(dim(FOI_R0_values$FOI))==3){
    FOI_R0_values=list(regions=FOI_R0_values$regions,FOI=rowMeans(FOI_R0_values$FOI,dims=2),
                       R0=rowMeans(FOI_R0_values$R0,dims=2))
  }

  #Convert data to frame
  n_regions=length(FOI_R0_values$regions)
  n_param_sets=dim(FOI_R0_values$FOI)[1]
  data_frame=data.frame(n_region=as.factor(rep(c(1:n_regions),n_param_sets)),
                        FOI=as.vector(t(FOI_R0_values$FOI)),R0=as.vector(t(FOI_R0_values$R0)))

  FOI_labels=10^c(-10:1)
  n_region=NULL

  if(plot_type %in% c("box", "violin")){
    FOI=R0=NULL
    p_FOI <- ggplot(data=data_frame, aes(x=n_region, y=log(FOI))) + theme_bw()
    if(plot_type=="box"){
      p_FOI <- p_FOI+geom_boxplot(outlier.size = 0)
    } else {
      p_FOI <- p_FOI+geom_violin(trim=FALSE, scale="width")
    }
    p_FOI <- p_FOI + scale_x_discrete(name="", breaks=c(1:n_regions), labels=FOI_R0_values$regions)
    p_FOI <- p_FOI + scale_y_continuous(name="FOI", breaks=log(FOI_labels), labels=FOI_labels)
    p_FOI <- p_FOI + theme(axis.text.x = element_text(angle = 90, hjust=1, size=text_size1),
                           axis.text.y = element_text(size = text_size1),
                           axis.title.y = element_text(size = text_size1))

    p_R0 <- ggplot(data=data_frame, aes(x=n_region, y=R0)) + theme_bw()
    if(plot_type=="box"){
      p_R0 <- p_R0+geom_boxplot(outlier.size = 0)
    } else {
      p_R0 <- p_R0+geom_violin(trim=FALSE, scale="width")
    }
    p_R0 <- p_R0 + scale_x_discrete(name="", breaks=c(1:n_regions), labels=FOI_R0_values$regions)
    p_R0 <- p_R0 + theme(axis.text.x = element_text(angle = 90, hjust=1, size=text_size1),
                         axis.text.y = element_text(size = text_size1),
                         axis.title.y = element_text(size = text_size1))

    output<-list(p_FOI=p_FOI, p_R0=p_R0)
  }

  if(plot_type=="scatter"){
    p_FOI_R0 <- ggplot() + theme_bw()
    p_FOI_R0 <- p_FOI_R0 + scale_x_continuous(name="FOI", breaks=log(FOI_labels), labels=FOI_labels)
    for(i in 1:n_regions){
      subset=data_frame[data_frame$n_region==i, ]
      names(subset)[names(subset)=="n_region"]="Region"
      Region=NULL
      p_FOI_R0 <- p_FOI_R0 + geom_point(data=subset, aes(x=log(FOI), y=R0, colour=Region))
      p_FOI_R0 <- p_FOI_R0 + theme(axis.text.x = element_text(size = text_size1),
                                   axis.text.y = element_text(size = text_size1),
                                   axis.title.y = element_text(size = text_size1))
    }
    output<-list(p_FOI_R0=p_FOI_R0)
  }

  if(plot_type=="error_bars"){
    blank=rep(NA, n_regions)
    lower=upper=NULL
    summary_frame_FOI=data.frame(n_region=c(1:n_regions), mean=blank, lower=blank, upper=blank)
    if(is.null(data_frame$R0)==FALSE){summary_frame_R0=summary_frame_FOI}
    for(i in 1:n_regions){
      subset=data_frame[data_frame$n_region==i, ]
      FOI_CI=exp(CI(log(subset$FOI)))
      summary_frame_FOI$mean[i]=FOI_CI[[2]]
      summary_frame_FOI$lower[i]=FOI_CI[[3]]
      summary_frame_FOI$upper[i]=FOI_CI[[1]]
      if(is.null(data_frame$R0)==FALSE){
        R0_CI=CI(subset$R0)
        summary_frame_R0$mean[i]=R0_CI[[2]]
        summary_frame_R0$lower[i]=R0_CI[[3]]
        summary_frame_R0$upper[i]=R0_CI[[1]]
      }
    }

    p_FOI <- ggplot(data=summary_frame_FOI, aes(x=n_region, y=log(mean))) + theme_bw()
    p_FOI <- p_FOI + scale_x_continuous(name="", breaks=c(1:n_regions), labels=FOI_R0_values$regions)
    p_FOI <- p_FOI + theme(axis.text.x = element_text(angle = 90, hjust=1, size=text_size1),
                           axis.text.y = element_text(size = text_size1))
    p_FOI <- p_FOI + scale_y_continuous(name="FOI", breaks=log(FOI_labels), labels=FOI_labels)
    p_FOI <- p_FOI + geom_line(data=summary_frame_FOI, aes(x=n_region, y=log(mean)))
    p_FOI <- p_FOI + geom_errorbar(data=summary_frame_FOI, aes(ymin=log(lower), ymax=log(upper)), width=0.5)

    p_R0 <- ggplot(data=summary_frame_R0, aes(x=n_region, y=mean)) + theme_bw()
    p_R0 <- p_R0 + scale_x_continuous(name="", breaks=c(1:n_regions), labels=FOI_R0_values$regions)
    p_R0 <- p_R0 + theme(axis.text.x = element_text(angle = 90, hjust=1, size=text_size1),
                         axis.text.y = element_text(size = text_size1))
    p_R0 <- p_R0 + geom_line(data=summary_frame_R0, aes(x=n_region, y=mean))
    p_R0 <- p_R0 + geom_errorbar(data=summary_frame_R0, aes(ymin=lower, ymax=upper), width=0.5)

    output<-list(p_FOI=p_FOI, p_R0=p_R0)
  }

  return(output)
}
#-------------------------------------------------------------------------------
#' @title plot_mcmc_enviro_coeff_data
#'
#' @description Plot environmental coefficients from MCMC results
#'
#' @details This function takes in a data frame of values (obtained using the get_mcmc_enviro_coeff_data() function) of
#'   the coefficients of environmental covariates and plots them on graphs
#'
#' @param data_frame Data frame of coefficient values obtained using get_mcmc_enviro_coeff_data()
#' @param env_vars List of environmental covariates
#' @param plot_type Type of plots to create (choose from "box", "violin")
#' @param text_size1 Text size parameter for axis labels
#' '
#' @export
#'
plot_mcmc_enviro_coeff_data <- function(data_frame=list(), env_vars=c(), plot_type="box", text_size1=10.0){
  assert_that(is.data.frame((data_frame)))
  assert_that(is.character(env_vars))
  assert_that(plot_type %in% c("box", "violin"))
  assert_that(is.numeric(text_size1))

  #TODO - Sort out environmental covariate numbers/labelling
  n_env_vars=length(env_vars)
  assert_that(n_env_vars==length(names(table(data_frame$n_env_var))))

  FOI_limits=c(floor(log(min(data_frame$FOI_coeffs), 10)), ceiling(log(max(data_frame$FOI_coeffs), 10)))
  FOI_labels=10^c(FOI_limits[1]:FOI_limits[2])
  R0_limits=c(floor(log(min(data_frame$R0_coeffs), 10)), ceiling(log(max(data_frame$R0_coeffs), 10)))
  R0_labels=10^c(R0_limits[1]:R0_limits[2])
  n_env_var=NULL

  if(plot_type %in% c("box", "violin")){
    FOI_coeffs=NULL
    p_FOI <- ggplot(data=data_frame, aes(x=n_env_var, y=log(FOI_coeffs, 10))) + theme_bw()
    if(plot_type=="box"){
      p_FOI <- p_FOI+geom_boxplot(outlier.size = 0)
    } else {
      p_FOI <- p_FOI+geom_violin(trim=FALSE, scale="width")}
    p_FOI <- p_FOI + scale_x_discrete(name="", breaks=c(1:n_env_vars), labels=env_vars)
    p_FOI <- p_FOI + scale_y_continuous(name="FOI coefficients", breaks=log(FOI_labels, 10), labels=FOI_labels,
                                        limits=FOI_limits)
    p_FOI <- p_FOI + theme(axis.text.x = element_text(angle = 90, hjust=1, size=text_size1),
                           axis.text.y = element_text(size = text_size1),
                           axis.title.y = element_text(size = text_size1))

    if(is.null(data_frame$R0_coeffs)==FALSE){
      R0_coeffs=NULL
      p_R0 <- ggplot(data=data_frame, aes(x=n_env_var, y=log(R0_coeffs, 10))) + theme_bw()
      if(plot_type=="box"){
        p_R0 <- p_R0+geom_boxplot(outlier.size = 0)
      } else {
        p_R0 <- p_R0+geom_violin(trim=FALSE, scale="width")
      }
      p_R0 <- p_R0 + scale_x_discrete(name="", breaks=c(1:n_env_vars), labels=env_vars)
      p_R0 <- p_R0 + scale_y_continuous(name="R0 coefficients", breaks=log(R0_labels, 10), labels=R0_labels,
                                        limits=R0_limits)
      p_R0 <- p_R0 + theme(axis.text.x = element_text(angle = 90, hjust=1, size=text_size1),
                           axis.text.y = element_text(size = text_size1),
                           axis.title.y = element_text(size = text_size1))
    } else {
      p_R0<-NULL
    }
    output<-list(p_FOI=p_FOI, p_R0=p_R0)
  }

  return(output)
}
#-------------------------------------------------------------------------------
#' @title plot_mcmc_prob_data
#'
#' @description Plot values of vaccine efficacy and/or reporting probabilities from MCMC output data
#'
#' @details This function takes in a data frame produced by functions get_mcmc_data() and/or truncate_mcmc_data(),
#' calculates the values of additional parameters fitted using the function besides spillover force of infection (FOI)
#' and reproduction number (R0) plots the values in one of a number of different ways to display the spread of values
#' obtained (box plots, violin plots, simple plots with error bars).
#'
#' @param input_frame Data frame of MCMC output data
#' @param plot_type Type of plots to create (choose from "box", "violin", "error_bars")
#' @param values List of names of parameters to plot (must be parameters appearing in the data)
#' @param text_size1 Text size parameter for axis labels
#'
#' @export
#'
plot_mcmc_prob_data <- function(input_frame=list(), plot_type="box", values=c("vaccine_efficacy"), text_size1=10.0){
  assert_that(is.data.frame((input_frame)))
  assert_that(plot_type %in% c("box", "violin", "error_bars"))
  assert_that(is.character(values))
  assert_that(is.numeric(text_size1))

  #TODO - Sort out variable naming/numbers
  n_values=length(values)
  for(i in 1:n_values){assert_that(values[i] %in% colnames(input_frame))}

  data<-input_frame[, colnames(input_frame) %in% values]
  n_lines=nrow(input_frame)
  labels=c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0)
  output_frame=data.frame(n_param=as.factor(rep(c(1:n_values), n_lines)), p=rep(NA, n_values*n_lines))
  for(i in 1:n_lines){
    output_lines=((i-1)*n_values)+c(1:n_values)
    output_frame$p[output_lines]=as.numeric(data[i, ])
  }

  if(plot_type %in% c("box", "violin")){
    n_param=p=NULL
    p_probs <- ggplot(data=output_frame, aes(x=n_param, y=log(p))) + theme_bw()
    if(plot_type=="box"){
      p_probs <- p_probs+geom_boxplot(outlier.size=0)
    } else {
      p_probs <- p_probs+geom_violin(trim=FALSE, scale="width")}
    p_probs <- p_probs + scale_x_discrete(name="", breaks=c(1:length(values)), labels=values)
    p_probs <- p_probs + scale_y_continuous(name="", breaks=log(labels), labels=labels)
    p_probs <- p_probs + theme(axis.text.x = element_text(angle = 90, hjust=1, size=text_size1),
                               axis.text.y = element_text(size = text_size1),
                               axis.title.y = element_text(size = text_size1))
  }
  if(plot_type=="error_bars"){
    blank=rep(NA, n_values)
    lower=upper=NULL
    summary_frame=data.frame(n_value=c(1:n_values), mean=blank, lower=blank, upper=blank)
    for(i in 1:n_values){
      subset=subset(output_frame, as.numeric(output_frame$n_param)==i)
      probs_CI=CI(subset$p)
      summary_frame$mean[i]=probs_CI[[2]]
      summary_frame$lower[i]=probs_CI[[3]]
      summary_frame$upper[i]=probs_CI[[1]]
    }
    n_value=NULL
    p_probs <- ggplot(data=summary_frame, aes(x=n_value, y=log(mean))) + theme_bw()
    p_probs <- p_probs + scale_x_continuous(name="", breaks=c(1:n_values), labels=values)
    p_probs <- p_probs + scale_y_continuous(name="", breaks=log(labels), labels=labels)
    p_probs <- p_probs + geom_line(data=summary_frame, aes(x=n_value, y=log(mean)))
    p_probs <- p_probs + geom_errorbar(data=summary_frame, aes(ymin=log(lower), ymax=log(upper)), width=0.5)
    p_probs <- p_probs + theme(axis.text.x = element_text(angle = 90, hjust=1, size=text_size1),
                               axis.text.y = element_text(size = text_size1),
                               axis.title.y = element_text(size = text_size1))

  }

  return(p_probs)
}
