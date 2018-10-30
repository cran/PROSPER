#' @title Surviving the herbicide
#' @seealso \code{\link{sel_resist}}

#' @description \code{sel_herb} calculates the surviving number of each genotype. \code{sel_herb} selects for resistant individuals.
#' @export

#' @template put
#' @template sdrate
#' @template rate
#' @template thresh
#' @template start
#' @template result

#' @template max_vec_length


#' @details For every genotype \code{\link{intern_herbicide}} is called.


#' @examples
#' struc_preparation(Rmx=10, n_loci=2, epis=0, dom=1)
#' gen_freq( af=c(0.01,0.8), n_seeds=10000)
#' sel_herb(start="initialSB", result="winter", 
#'                        thresh=20, sdrate=0.4, rate=100, put=0.04)

sel_herb <-
function(start, result, thresh, sdrate, rate, put, max_vec_length=1e+07){
###------- value check of genotype
dfgenotype <- get0("dfgenotype", envir = parent.frame(n = 1))
if(is.null(dfgenotype$resist)){
                             eval(parse(text=paste("dfgenotype$",result,"<-", "dfgenotype$",start,sep="")))
                             return(dfgenotype)
                             }#END if(is.null)

###------- preparing the usage of herbicide()
 first_amount <- eval(parse(text=paste("dfgenotype$",start,sep="")))   #the amount that gives the start
 return_amount <- rep(0,nrow(dfgenotype))
       for(pres_GT in which(first_amount > 0)){
            i1 <- first_amount[pres_GT] %/% max_vec_length
            i2 <- first_amount[pres_GT] %% max_vec_length
            for(i in seq_len(i1)){
                     return_amount[pres_GT] <- intern_herbicide(resist=dfgenotype$resist[pres_GT], n_samples=max_vec_length, put=put, rate= rate, sdrate=sdrate, thresh=thresh)
                     }#END                  
            return_amount[pres_GT] <- intern_herbicide(resist=dfgenotype$resist[pres_GT], n_samples=i2, put=put, rate= rate, sdrate=sdrate, thresh=thresh)
               }#END for(pres_GT)
 eval(parse(text=paste("dfgenotype$",result,"<-", data.frame(return_amount),sep="")))            #for saving the result
 assign("dfgenotype", value=dfgenotype, pos = -1, envir=parent.frame(n = 1))            
 return(dfgenotype)
}
