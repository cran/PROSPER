#' @title Checking the main parameters
#' \code{link{stuc_preparation}}

#' @description Before the preparation (\code{struc_preparation}) it is necessary to check, whether the necessary information is provided to principally conduct a simulation run. Furthermore the crop.list is brought to the correct form.

#' @template param.weed
#' @template area
#' @template duration
#' @template repetitions
#' @template crop_list
#' @template max_vec_length


#' @details If no genetics are required \code{af} has to be set to \code{NA}. The value of \code{dom} is adjusted to the length of \code{af}. The \code{dom} is cut or the first given value is repeated to fit the number of given alleles in \code{af}. The variables af and \code{dom} are corrected and \code{n_loci} is created, which is 0 when no genetics is included.



mod_check <- function(param.weed, area, duration, repetitions, crop_list, max_vec_length){
cat("Check general model settings... ")
if(!is.data.frame(param.weed)){stop("param.weed need to be defined as data.frame")}
if(is.na(area)){stop("area need to be defined")}
if(is.na(duration)){stop("duration need to be defined")}
if(is.na(crop_list)){stop("crop_list need to be defined")}
if(is.na(max_vec_length)) {warning("max_vec_length is set to 1e+07")
                          assign("max_vec_length", 1e+07, pos = -1, envir=parent.frame(n = 1))
                          assign("max_vec_length", 1e+07, pos = -1, envir=parent.frame(n = 2))
                          }

#-- creating the crop.list object ------------------------------------------------------
crop.list <- rep(crop_list,ceiling(duration/length(crop_list)))[1:duration] #setting for crop rotation
assign("crop.list", value=crop.list, pos = -1, envir=parent.frame(n = 1))
cat(" finished!\n")
#return()
}