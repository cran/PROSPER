#' @title Surviving a non-selective process
#' @seealso \code{\link{quanti}} \code{\link{pop_germc}} 

#' @export
#' @description \code{pop_step} picks the individuals that will pass to the next development stage. This is a random process for every individual, which does not exert any selection pressure. 

#  @template dfgenotype
#' @template start
#' @template result
#' @template max_vec_length  
#' @template start_comb
#' @param stepname name of the new column of dfgenotype added by this function. \code{character}.
#' @param surv_prob probability to survive this step and reach the next growth stage. \code{numeric}. 

#' @details Individuals that reach the next growth stage are picked by using \code{\link{rbinom}}. In contrast to \code{\link{sel_herb}}, \code{pop_step} does not exert any evolutionary selection pressure.
#' When more than one column is selected with \code{start}, they are summed and the result is passed to the picking process. By setting \code{start_comb} the sum is added as a column to \code{dfgenotype}.


#' @return A new column is added to \code{dfgenotype} containing the surviving individuals of the different genotypes.

#' @examples
#' struc_preparation(Rmx=10, n_loci=2, epis=0, dom=1)
#' gen_freq(af=c(0.01,0.8), n_seeds=10000)
#' #How many individuals of each genotype will reach the next growth stage?
#' pop_step(start="initialSB",  stepname="survivingthewinter",
#'                             surv_prob=0.4)

pop_step <-
function(start, start_comb = NA, result = NA, stepname = NA, surv_prob, max_vec_length=1e+07){
cat("pop_step starts...")
if(is.na(result) & is.na(stepname)){warning("In pop_step() either result or stepname should be given.")}
if(is.na(result)){result <- paste("surv_",stepname,sep="")}         # standard name

dfgenotype <- get0("dfgenotype", envir = parent.frame(n = 1))


### --- name flexibel
first_amount <- data.frame(matrix(nrow=nrow(dfgenotype),ncol=length(start)))
names(first_amount) <- start
second_amount <- rep(0,nrow(dfgenotype))                                          #returning amounts of GT

for(cohort in seq_along(start)){ 
 first_amount[cohort] <- eval(parse(text=paste("dfgenotype$",start[cohort],sep="")))   #the amount that gives the start
}#END for(cohort)
first_amount <- rowSums(first_amount)
if(length(start)>1 & !is.na(start_comb)){
                   eval(parse(text=paste("dfgenotype$",start_comb,"<-first_amount",sep="")))                
                   }#END if(length(start))

for(pres_GT in which(first_amount > 0)){
  i1 <- first_amount[pres_GT] %/% max_vec_length
  i2 <- first_amount[pres_GT] %% max_vec_length
    for(i in seq(len=i1)){
      second_amount[pres_GT] <- second_amount[pres_GT] + rbinom(1, max_vec_length, prob = surv_prob)
    }#END for(i)
  second_amount[pres_GT] <- second_amount[pres_GT] + rbinom(1, i2, prob = surv_prob)   #sum up all surviving weeds
}#END for(j)


eval(parse(text=paste("dfgenotype$",result,"<-","second_amount",sep="")))  
assign("dfgenotype", value=dfgenotype, pos = -1, envir=parent.frame(n = 1)) 
cat("finished!\n")
return(dfgenotype)
}
