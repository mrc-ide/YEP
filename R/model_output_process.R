#-------------------------------------------------------------------------------
#' @title convert_model_output_combine_by_age
#'
#' @description Convert SEIRV output to combine values across ages
#'
#' @details Takes in the output of the Model_Run function and sums together
#' S, E, I, R, and V values across age groups, outputting a list containing arrays of values by particle and time point
#'
#' @param model_output List of output data produced by appropriate functions
#' '
#' @export
#'
convert_model_output_combine_by_age <- function(model_output=list()){

  assert_that(is.list(model_output)) #TODO - msg
  assert_that(is.null(model_output$day)==FALSE) #TODO - msg

  N_age=dim(model_output$S)[2]
  t_pts=dim(model_output$S)[1]
  blank=rep(0,t_pts)
  S_sum=E_sum=I_sum=R_sum=V_sum=C_sum=blank
  for(i in 1:N_age){
    S_sum=S_sum+model_output$S[,i]
    E_sum=E_sum+model_output$E[,i]
    I_sum=I_sum+model_output$I[,i]
    R_sum=R_sum+model_output$R[,i]
    V_sum=V_sum+model_output$V[,i]
    C_sum=C_sum+model_output$C[,i]
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

  values=label=NULL
  t_pts=dim(model_output$S)[1]
  date_values=model_output$year[1]+((model_output$day-model_output$day[1]+1)/365.0)
  xlabels=c(model_output$year[1]:(model_output$year[t_pts]+1))
  ylabels=10^c(-1:10)

  data_combined=data.frame(date_values=rep(date_values,5),
                           label=c(rep("S",t_pts),rep("E",t_pts),rep("I",t_pts),rep("R",t_pts),rep("V",t_pts)),
                           values=c(model_output2$S,model_output2$E,model_output2$I,model_output2$R,model_output2$V))
  data_combined$values[data_combined$values<0.1]=0
  plot1 <- ggplot(data=data_combined,aes(x=date_values,y=log(values),group=label))+theme_bw()
  plot1 <- plot1 + geom_line(aes(colour=label),size=1) + theme(legend.title=element_blank())
  plot1 <- plot1 + scale_x_continuous(name="",breaks=xlabels,labels=xlabels)
  plot1 <- plot1 + scale_y_continuous(name="",breaks=log(ylabels),labels=ylabels)

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

  N_age=dim(model_output$S)[2]
  t_pts=dim(model_output$S)[1]
  date_values1=model_output$year[1]+((model_output$day-model_output$day[1]+1)/365.0)
  date_values2=array(rep(date_values1,N_age),dim=c(N_age,t_pts))

  return(data.frame(age=rep(c(1:N_age),t_pts),
                    date=as.vector(date_values2),S=as.vector(model_output$S),E=as.vector(model_output$E),
                    I=as.vector(model_output$I), R=as.vector(model_output$R),V=as.vector(model_output$V)))
}
