#' @title Preparing the data for the next simulation run
#' @seealso \code{\link{struc_endSim}} \code{\link{struc_saveSimData}}
#' @description \code{struc_nextCycle} finishes the simulation cycle. Therefore, the \code{dfgenotype} is reset and prepared for the next cycle.
#' @export

#' @param start_name name of the column with the start population for the next simulation cycle. If missing the column name after the column 'resist' will be used. \code{character}.
#' @param end_name name of the column with the population at the end of the simulation cycle. This will be used as the start values in the next simulation cycle (in the column with \code{start_name}). If missing the last column name is used. \code{character}.

#' @return A \code{data.table} object.


struc_nextCycle <- 
function(start_name, end_name){
simdataob <- get0("sim_result", envir = parent.frame(n = 1))
    
    rep    <- "repetition" 
    simcyc <- "simcycle"
    GTcol  <- "genotype"
    namestokeep <- c("genotype", "resist")
    simd <- intern_SimData(simdataob)
    simn <- names(simd)
cat("\n struc_nextCycle 1...")
#- setting the stnadart start and end columnname
   if(missing(start_name)) start_name <- simn[which(simn=="resist")+1]
   if(missing(end_name))   end_name   <- simn[length(simn)]
#- selection of the last computed simulation cycle    
    latest_rep <- simd[,max(.SD),.SDcols=rep]
    latest_cyc <- simd[eval(parse(text=rep))==latest_rep,  max(eval(parse(text=simcyc)))]           
cat(" 2...")    
      nextsim <-      simd[eval(parse(text=simcyc))==latest_cyc & eval(parse(text=rep))==latest_rep & !is.na(eval(parse(text=GTcol))),][,.SD, .SDcols=c(namestokeep, end_name)] 

    names(nextsim) <- c(namestokeep, start_name)
    assign("dfgenotype", value=nextsim, pos = -1, envir=parent.frame(n = 1))   
cat(" finished\n")                 
    invisible(nextsim)
}




