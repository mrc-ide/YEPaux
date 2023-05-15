# R file for general functions in YFDaux package
#------------------------------------------------
#Global variables
p_severe_inf=0.12 #Probability that an infection is severe
p_death_severe_inf=0.39 #Probability that a severe infection becomes fatal
t_incubation <- 5 #Time for cases to incubate in mosquito
t_latent <- 5 #Latent period before cases become infectious
t_infectious <- 5 #Time cases remain infectious
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.
#' @importFrom assertthat assert_that
#' @importFrom coda gelman.diag gelman.plot mcmc mcmc.list
#' @importFrom ggplot2 aes element_blank element_text geom_boxplot geom_errorbar geom_line geom_point geom_ribbon geom_violin ggplot labs
#' @importFrom ggplot2 scale_x_continuous scale_x_discrete scale_y_continuous theme theme_bw
#' @importFrom graphics axis legend matplot par title
#' @importFrom grDevices dev.off png
#' @importFrom Rmisc CI
#' @importFrom sf read_sf st_bbox st_geometry
#' @importFrom stats cov dexp dnbinom dnorm prop.test rbinom runif
#' @importFrom utils read.csv write.csv
#' @import YEP
#-------------------------------------------------------------------------------
#' @title convert_model_output_combine_by_age # TODO - Amend to deal with different outputs
#'
#' @description Convert SEIRV output of other functions to combine values across ages
#'
#' @details Takes in the output of YEP::Model_Run and similar functions and sums together
#' S, E, I, R, and V values across age groups, outputting a list containing arrays of values by particle and time point
#'
#' @param model_output List of output data produced by appropriate functions
#' '
#' @export
#'
convert_model_output_combine_by_age <- function(model_output=list()){

  assert_that(is.list(model_output)) #TODO - msg
  assert_that(is.null(model_output$day)==FALSE) #TODO - msg

  N_age=dim(model_output$S)[1]
  n_particles=dim(model_output$S)[2]
  t_pts=dim(model_output$S)[3]
  blank=array(0,dim=c(n_particles,t_pts))
  S_sum=E_sum=I_sum=R_sum=V_sum=C_sum=blank
  for(i in 1:N_age){
    S_sum=S_sum+model_output$S[i,,]
    E_sum=E_sum+model_output$E[i,,]
    I_sum=I_sum+model_output$I[i,,]
    R_sum=R_sum+model_output$R[i,,]
    V_sum=V_sum+model_output$V[i,,]
    C_sum=C_sum+model_output$C[i,,]
  }

  return(list(day=model_output$day,year=model_output$year,S=S_sum,E=E_sum,I=I_sum,R=R_sum,V=V_sum,C=C_sum))
}
#-------------------------------------------------------------------------------
#' @title plot_model_output
#'
#' @description Plot outputs of SEIRV model combined across age groups
#'
#' @details Takes in SEIRV output from Full_Model_Run(), Basic_Model_Run() and similar functions and outputs a ggplot
#'   graph of S, R and V summed over all age groups, with error bars showing the spread of values from multiple
#'   particles. If data for a single particle is selected, E and I values are also shown.
#'
#' @param model_output List of SEIRV output data produced by appropriate functions
#' '
#' @export
#'
plot_model_output <- function(model_output=list()){

  model_output2=convert_model_output_combine_by_age(model_output)

  values=label=values_mean=values_low=values_high=NULL
  n_particles=dim(model_output$S)[2]
  t_pts=dim(model_output$S)[3]
  date_values=model_output$year[1]+((model_output$day-model_output$day[1]+1)/365.0)
  xlabels=c(model_output$year[1]:(model_output$year[t_pts]+1))
  ylabels=10^c(-1:10)

  if(n_particles==1){
    data_combined=data.frame(date_values=rep(date_values,5),
                             label=c(rep("S",t_pts),rep("E",t_pts),rep("I",t_pts),rep("R",t_pts),rep("V",t_pts)),
                             values=c(model_output2$S[1,],model_output2$E[1,],model_output2$I[1,],model_output2$R[1,],
                                      model_output2$V[1,]))
    data_combined$values[data_combined$values<0.1]=0
    plot1 <- ggplot(data=data_combined,aes(x=date_values,y=log(values),group=label))+theme_bw()
    plot1 <- plot1 + geom_line(aes(colour=label),size=1) + theme(legend.title=element_blank())
    plot1 <- plot1 + scale_x_continuous(name="",breaks=xlabels,labels=xlabels)
    plot1 <- plot1 + scale_y_continuous(name="",breaks=log(ylabels),labels=ylabels)
  } else {
    data_combined=data.frame(date_values=rep(date_values,3),label=c(rep("S",t_pts),rep("R",t_pts),rep("V",t_pts)),
                             values_mean=rep(NA,t_pts*3),values_low=rep(NA,t_pts*3),values_high=rep(NA,t_pts*3))
    for(i in 1:t_pts){
      S_CI=CI(model_output2$S[,i])
      R_CI=CI(model_output2$R[,i])
      V_CI=CI(model_output2$V[,i])
      data_combined$values_mean[i]=as.numeric(S_CI[2])
      data_combined$values_mean[i+t_pts]=as.numeric(R_CI[2])
      data_combined$values_mean[i+(2*t_pts)]=as.numeric(V_CI[2])
      data_combined$values_low[i]=as.numeric(S_CI[3])
      data_combined$values_low[i+t_pts]=as.numeric(R_CI[3])
      data_combined$values_low[i+(2*t_pts)]=as.numeric(V_CI[3])
      data_combined$values_high[i]=as.numeric(S_CI[1])
      data_combined$values_high[i+t_pts]=as.numeric(R_CI[1])
      data_combined$values_high[i+(2*t_pts)]=as.numeric(V_CI[1])
    }
    data_combined$values_low=log(data_combined$values_low)
    data_combined$values_low[is.infinite(data_combined$values_low)]=NA
    data_combined$values_high=log(data_combined$values_high)
    data_combined$values_high[is.infinite(data_combined$values_high)]=NA
    plot1 <- ggplot(data=data_combined,aes(x=date_values,y=log(values_mean),group=label))+theme_bw()
    plot1 <- plot1 + geom_line(aes(colour=label),size=1) + theme(legend.title=element_blank())
    plot1 <- plot1 + geom_errorbar(data=data_combined,aes(ymin=values_low,ymax=values_high,colour=label),
                                   width=(model_output$day[2]-model_output$day[1])/365.0)
    plot1 <- plot1 + scale_x_continuous(name="",breaks=xlabels,labels=xlabels)
    plot1 <- plot1 + scale_y_continuous(name="",breaks=log(ylabels),labels=ylabels)
  }

  return(plot1)
}
#-------------------------------------------------------------------------------
#' @title convert_model_output_tidy
#'
#' @description Convert output of Full_Model_Run() or Basic_Model_Run() functions to simple line list
#'
#' @details Takes in the output of the Full_Model_Run() or Basic_Model_Run() functions and outputs a data frame
#' with one set of SEIRV data per line, for use with tidyr functions
#'
#' @param model_output List of output data produced by appropriate functions
#' '
#' @export
#'
convert_model_output_tidy <- function(model_output=list()){

  assert_that(is.list(model_output)) #TODO - msg
  assert_that(is.null(model_output$day)==FALSE) #TODO - msg

  N_age=dim(model_output$S)[1]
  n_particles=dim(model_output$S)[2]
  t_pts=dim(model_output$S)[3]
  date_values1=model_output$year[1]+((model_output$day-model_output$day[1]+1)/365.0)
  particle_values=date_values2=model_output$S
  for(i in 1:n_particles){
    particle_values[,i,]=array(i,dim=c(N_age,t_pts))
    date_values2[,i,]=array(rep(date_values1,N_age),dim=c(N_age,t_pts))
  }

  output_frame=data.frame(age=rep(c(1:N_age),n_particles*t_pts),particle=as.vector(particle_values),
                          date=as.vector(date_values2),S=as.vector(model_output$S),E=as.vector(model_output$E),
                          I=as.vector(model_output$I), R=as.vector(model_output$R),V=as.vector(model_output$V))

  return(output_frame)
}
