# R file for functions relating to the creation/formatting of input data for the yellow fever model
#-------------------------------------------------------------------------------
#' @title create_input_data
#'
#' @description Creates input data set in correct format for use by other functions
#'
#' @details Takes in vaccination and population data in data frames (in columns by age with columns showing the region
#'   and year for each row), extracts number of age groups (verifying that this is the same in each data frame),
#'   extracts data for specified regions and years, and creates list in format used by other functions (vectors of
#'   region names, years and age groups, 3-dimensional arrays of vaccination and population data).
#'
#' @param vacc_data Data frame containing vaccination coverage data with region in column 1, year in column 2 and
#'   coverage values by age in remaining columns
#' @param pop_data Data frame containing population data with region in column 1, year in column 2 and
#'   population values by age in remaining columns
#' @param regions Vector of regions for which to extract data from vacc_data and pop_data (in alphabetical order)
#' @param years Vector of years for which to extract data from vacc_data and pop_data
#' '
#' @export
#'
create_input_data <- function(vacc_data=list(),pop_data=list(),regions=c(),years=c()){

  #Check data
  assert_that(is.data.frame(vacc_data))
  assert_that(is.data.frame(pop_data))
  assert_that(is.character(regions))
  assert_that(all(regions==sort(regions)),msg="Regions must be ordered by name")
  assert_that(is.numeric(years))
  assert_that(ncol(pop_data)==ncol(vacc_data),msg="Numbers of columns in vaccination and population data must match")
  vacc_regions=unique(vacc_data[,1])
  assert_that(all(vacc_regions==sort(vacc_regions)),msg="Regions must be ordered by name in vaccination data")
  vacc_years=unique(vacc_data[,2])
  assert_that(all(regions %in% vacc_regions),msg="All specified regions must appear in vaccination data")
  assert_that(all(years %in% vacc_years),msg="All specified years must appear in vaccination data")
  pop_regions=unique(pop_data[,1])
  assert_that(all(pop_regions==sort(pop_regions)),msg="Regions must be ordered by name in population data")
  pop_years=unique(pop_data[,2])
  assert_that(all(regions %in% pop_regions),msg="All specified regions must appear in population data")
  assert_that(all(years %in% pop_years),msg="All specified years must appear in population data")

  N_age=ncol(vacc_data)-2
  n_regions=length(regions)
  n_years=length(years)

  #Subset data and sort into correct order
  vacc_data_subset=subset(vacc_data,vacc_data[,1] %in% regions)
  vacc_data_subset=subset(vacc_data_subset,vacc_data_subset[,2] %in% years)
  vacc_data_subset=vacc_data_subset[order(vacc_data_subset[,1]),]
  vacc_data_subset=vacc_data_subset[order(vacc_data_subset[,2]),]
  pop_data_subset=subset(pop_data,pop_data[,1] %in% regions)
  pop_data_subset=subset(pop_data_subset,pop_data_subset[,2] %in% years)
  pop_data_subset=pop_data_subset[order(pop_data_subset[,1]),]
  pop_data_subset=pop_data_subset[order(pop_data_subset[,2]),]

  #Organize vaccination and population data into arrays
  vacc_coverage_array=pop_array=array(data=rep(0,n_regions*n_years*N_age),dim=c(n_regions,n_years,N_age))
  for(n_region in 1:n_regions){
    region=regions[n_region]
    lines=vacc_data_subset[,1]==region
    vacc_data_subset2=vacc_data_subset[lines,]
    pop_data_subset2=pop_data_subset[lines,]
    vacc_coverage_array[n_region,,]=as.matrix(vacc_data_subset2[,c(3:(2+N_age))])
    pop_array[n_region,,]=as.matrix(pop_data_subset2[,c(3:(2+N_age))])
  }

  #Create and output dataset
  return(list(region_labels=regions,years_labels=years,age_labels=c(0:(N_age-1)),
              vacc_data=vacc_coverage_array,pop_data=pop_array))
}
#-------------------------------------------------------------------------------
#' @title input_data_check
#'
#' @description Check that input data is correctly formatted
#'
#' @details This function takes in a list of input data for use with other functions and checks that it is correctly
#'   formatted, including containing all necessary elements and having years and ages in sequence
#'
#' @param input_data List of population and vaccination data for multiple regions
#'
#' @export
#'
input_data_check <- function(input_data=list()){

  result=TRUE

  assert_that(is.list(input_data))
  assert_that(is.vector(input_data$region_labels))
  n_regions=length(input_data$region_labels)

  assert_that(is.vector(input_data$years_labels))
  n_years=length(input_data$years_labels)
  for(i in 2:n_years){
    assert_that(input_data$years_labels[i]==input_data$years_labels[i-1]+1,
                msg="Year labels in input data must be continuous series of consecutive years")
    }

  assert_that(is.vector(input_data$age_labels))
  N_age=length(input_data$age_labels)
  for(i in 2:N_age){
    assert_that(input_data$age_labels[i]==input_data$age_labels[i-1]+1,
                msg="Age labels in input data must be continuous series of consecutive ages")
    }

  assert_that(length(dim(input_data$vacc_data))==3,msg="Vaccination data must be 3-dimensional array")
  assert_that(dim(input_data$vacc_data)[1]==n_regions,
              msg="First dimension of vaccination data must equal number of region labels")
  assert_that(dim(input_data$vacc_data)[2]==n_years,
              msg="Second dimension of vaccination data must equal number of year labels")
  assert_that(dim(input_data$vacc_data)[3]==N_age,
              msg="Third dimension of vaccination data must equal number of age labels")

  assert_that(length(dim(input_data$pop_data))==3,msg="Population data must be 3-dimensional array")
  assert_that(dim(input_data$pop_data)[1]==n_regions,
              msg="First dimension of population data must equal number of region labels")
  assert_that(dim(input_data$pop_data)[2]==n_years,
              msg="Second dimension of population data must equal number of year labels")
  assert_that(dim(input_data$pop_data)[3]==N_age,
              msg="Third dimension of population data must equal number of age labels")

  return(result)
}
#-------------------------------------------------------------------------------
#' @title input_data_process
#'
#' @description Cross-reference input data with serological and/or annual case/death data for comparison
#'
#' @details This function, used to prepare input data for functions used to calculate the likelihood of observed data,
#'   amends a list of population and vaccination data used as input for other functions,
#'   cross-referencing it with seroprevalence and/or case data and adding connection information.
#'
#' @param input_data List of population and vaccination data for multiple regions (created using create_input_data()
#'   function and usually loaded from an RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#'
#' @export
#'
input_data_process <- function(input_data=list(),obs_sero_data=NULL,obs_case_data=NULL){

  assert_that(input_data_check(input_data),msg="Input data must be in correct format")
  if(is.null(obs_sero_data)==FALSE){
    assert_that(max(obs_sero_data$year)<max(input_data$years_labels),msg="Input data years must encompass sero data range")
  }
  if(is.null(obs_case_data)==FALSE){
    assert_that(max(obs_case_data$year)<max(input_data$years_labels),msg="Input data years must encompass case data range")
  }
  N_age=length(input_data$age_labels)
  n_years=length(input_data$years_labels)

  regions1=input_data$region_labels
  n_regions1=length(regions1)
  assert_that(all(regions1==sort(regions1)),msg="Region labels must be in alphabetical order")
  regions_sero_com=unique(obs_sero_data$region)
  regions_case_com=unique(obs_case_data$region)
  regions_sero_unc=regions_case_unc=c()

  for(region in regions_sero_com){regions_sero_unc=append(regions_sero_unc,strsplit(region,",")[[1]])}
  for(region in regions_case_com){regions_case_unc=append(regions_case_unc,strsplit(region,",")[[1]])}
  regions_sero_unc=unique(regions_sero_unc)
  regions_case_unc=unique(regions_case_unc)
  regions_all_data_unc=unique(c(regions_sero_unc,regions_case_unc))

  for(region in regions_all_data_unc){
    if(region %in% regions1==FALSE){
      cat("\nInput data error - ",region," does not appear in input data\n",sep="")
      stop()
    }
  }

  input_regions_check=regions1 %in% regions_all_data_unc
  regions2=regions1[input_regions_check]
  n_regions2=length(regions2)

  # #TODO - fix mechanism for skipping if update not needed
  # if(n_regions2==n_regions1 && is.null(input_data$year_data_begin)==FALSE){
  #   #Skip steps if the input data already has the same set of regions and the cross-reference tables
  #   #assert_that() #Check all tables present
  #   return(input_data)
  # }else{

    blank=rep(0,n_regions2)
    flag_sero=flag_case=year_end=blank
    year_data_begin=rep(Inf,n_regions2)
    sero_line_list=case_line_list=list()
    for(i in 1:n_regions2){
      region=regions2[i]
      if(region %in% regions_sero_unc){
        flag_sero[i]=1
        sero_line_list[[i]]=c(0)
        for(j in 1:nrow(obs_sero_data)){
          if(grepl(region,obs_sero_data$region[j])==TRUE){
            sero_line_list[[i]]=append(sero_line_list[[i]],j)
            year_data_begin[i]=min(obs_sero_data$year[j],year_data_begin[i])
            year_end[i]=max(obs_sero_data$year[j],year_end[i])
          }
        }
        sero_line_list[[i]]=sero_line_list[[i]][c(2:length(sero_line_list[[i]]))]
      }
      if(region %in% regions_case_unc){
        flag_case[i]=1
        case_line_list[[i]]=c(0)
        for(j in 1:nrow(obs_case_data)){
          if(grepl(region,obs_case_data$region[j])==TRUE){
            case_line_list[[i]]=append(case_line_list[[i]],j)
            year_data_begin[i]=min(obs_case_data$year[j],year_data_begin[i])
            year_end[i]=max(obs_case_data$year[j],year_end[i])
          }
        }
        case_line_list[[i]]=case_line_list[[i]][c(2:length(case_line_list[[i]]))]
      }
    }

    return(list(region_labels=input_data$region_labels[input_regions_check],
                years_labels=input_data$years_labels,age_labels=input_data$age_labels,
                vacc_data=array(input_data$vacc_data[input_regions_check,,],dim=c(n_regions2,n_years,N_age)),
                pop_data=array(input_data$pop_data[input_regions_check,,],dim=c(n_regions2,n_years,N_age)),
                year_data_begin=year_data_begin,year_end=year_end,flag_sero=flag_sero,flag_case=flag_case,
                sero_line_list=sero_line_list,case_line_list=case_line_list))
  #}
}

#-------------------------------------------------------------------------------
#' @title input_data_truncate
#'
#' @description Truncate input data list for shorter set of regions
#'
#' @details TBA
#'
#' @param input_data List of population and vaccination data for multiple regions (created using create_input_data()
#'   function and usually loaded from an RDS file)
#' @param regions_new Vector of regions (subset of input_data$region_labels) for which to create new, shorter dataset
#'
#' @export
#'
input_data_truncate <- function(input_data=list(),regions_new=c()){

  assert_that(input_data_check(input_data))
  N_age=length(input_data$age_labels)
  n_years=length(input_data$years_labels)

  assert_that(all(input_data$region_labels==sort(input_data$region_labels)),msg="Region labels must be in alphabetical order")
  assert_that(all(regions_new==sort(regions_new)),msg="Specified regions must be in alphabetical order")
  assert_that(all(regions_new %in% input_data$region_labels),msg="Specified regions must be present in dataset")

  input_regions_check=input_data$region_labels %in% regions_new
  n_regions=length(input_regions_check[input_regions_check==TRUE])

  input_data_new=list(region_labels=input_data$region_labels[input_regions_check],
                      years_labels=input_data$years_labels,age_labels=input_data$age_labels,
                      vacc_data=array(input_data$vacc_data[input_regions_check,,],dim=c(n_regions,n_years,N_age)),
                      pop_data=array(input_data$pop_data[input_regions_check,,],dim=c(n_regions,n_years,N_age)))

  return(input_data_new)
}
#-------------------------------------------------------------------------------
#' @title regions_breakdown
#'
#' @description Break down regions from datasets to get list of unique regions
#'
#' @details Takes in vector of region labels, potentially including labels which contain more than one region (e.g.
#'   labels for countrywide data with all adm1 level regions as one label separated by comma) and produces an alphabetically
#'   ordered list of region labels, breaking down comma-separated groups.
#'
#'   For example, supplying a vector where all labels are "BRA.1_1,BRA.2_1,BRA.3_1" will return a vector of length 3
#'   - c("BRA.1_1","BRA.2_1","BRA.3_1").
#'
#' @param region_labels Vector of region labels
#'
#' @export
#'
regions_breakdown <- function(region_labels=c()){
  #TODO - assert_that

  region_labels_unique=unique(region_labels)

  regions_out=c()
  for(region in region_labels_unique){regions_out=append(regions_out,strsplit(region,",")[[1]])}
  regions_out=unique(sort(regions_out))

  return(regions_out)
}
#-------------------------------------------------------------------------------
#' @title template_region_xref
#'
#' @description Cross-reference template data with individual regions
#'
#' @details Examines template data (serological, case or burden) and compares with vector of region names in order
#'  to check which lines in each set of template data require model data from which region(s). For example, if a
#'  line in a set of serological data has its region given as "AGO.1_1,AGO.2_1,AGO.3_1" and the compared vector of
#'  regions is c("AGO.1_1","AGO.2_1","AGO.3_1",...), then that line requires data from regions 1, 2 and 3.
#'
#'  This function is used when generating a dataset from one or more templates; it is normally used by the functions
#'  Generate_Dataset, Generate_Sero_Dataset, Generate_Case_Dataset, Generate_VIMC_Burden_Dataset and
#'  Generate_Multiple_Datasets. It returns a list containing [TBA].
#'
#' @param template List containing one or more sets of template data (serological data, case data or burden data)
#' @param regions Vector of individual regions
#'
#' @export
#'
template_region_xref <- function(template=list(),regions=c()){
  #TODO - assert_that functions

  n_lines=nrow(template)
  n_regions=length(regions)
  output=list(line_list=list(),year_data_begin=rep(Inf,n_regions),year_end=rep(-Inf,n_regions))

  for(n_region in 1:n_regions){
    region=regions[n_region]
    output$line_list[[n_region]]=NA
    flag_started=FALSE
    for(n_line in 1:n_lines){
      if(grepl(region,template$region[n_line])==TRUE){
        if(flag_started){output$line_list[[n_region]]=append(output$line_list[[n_region]],n_line)} else {
          output$line_list[[n_region]]=c(n_line)
          flag_started=TRUE
        }
        output$year_data_begin[n_region]=min(template$year[n_line],output$year_data_begin[n_region])
        output$year_end[n_region]=max(template$year[n_line],output$year_end[n_region])
      }
    }
  }

  return(output)
}
