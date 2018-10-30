#' @title Check for some stop conditions.
#' @seealso \code{\link{struc_nextCycle}} \code{\link{struc_saveSimData}}

#' @description Conditions that are necessary to continue the simulation are checked.
#' @export
#' @param simcycle The count of simulation runs. The  value is \code{duration} \eqn{\ge  0}.
#' @param break_col_name The name of the column in \code{dfgenotype}, which is the first in the next simulation cycle. 


#' @details This function combines two conditions to terminate the simulation. The first is the \code{duration} that is defined for the simulation. The second is the the number of individuals for the next simulation cycle: when the population is extinct the simulation ends.

#' @return \code{logical} TRUE if the simulation has to stop. 



struc_endSim <- function (simcycle=year, break_col_name="SB_autumn") {

cat("The simulation cycle ended.\n")
dfgenotype <- get0("dfgenotype", envir = parent.frame(n = 1))
duration <- get0("duration", envir = parent.frame(n = 1))
break_sim <- FALSE
if(simcycle == duration){break_sim <- TRUE}
beak_col_sum <- eval(parse(text=paste("sum(dfgenotype$",break_col_name,")",sep="")))
if(beak_col_sum == 0) {break_sim <- TRUE}

ifelse(break_sim==TRUE,return(TRUE),return(FALSE))
}
