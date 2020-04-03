#' @title Germination
#' @seealso \code{\link{pop_reprod}} \code{\link{pop_step}} \code{\link{quanti}}  
#' @description \code{pop_germc} describes germination as random event for individual seeds. The function considers different cohorts and dormancy.
#' @export
#' @param germ \code{numeric} germination probabilities for the different cohorts. The sum of \code{germ} shall be \eqn{\ge 0} and \eqn{\le 1}. The difference of 1 and the sum of \code{germ} is the probability of dormancy. \code{numeric vector} 
#' @template max_vec_length
#' @param init_sb column name of the initial seed bank in the data.frame dfgenotype that is delivered by struc_preparation. \code{character}. 



#' @details Each individual has a chance to germinate or to stay dormant. In case of germination, it emerges in one of the cohorts. The distribution of individual seeds to cohorts or dormancy is random. The function uses the column \code{initialSB} of the data.frame \code{dfgenotype} as input. The output values are the numbers of seedlings of each genotype and each cohort.



#' @return Columns are added to \code{dfgenotype}: "germ_dorm" contains the numbers of each genotype that remain dormant, "germ1" to "germX" contain the numbers for each of X cohorts.

#' @examples 
#' struc_preparation2(Rmx=10, af=c(0.01,0.8), epis=0, dom=1)
#' ls()
#' gen_freq( af=c(0.01,0.8), n_seeds=10000)

#' #Distribute the individuals to three cohorts with the germination 
#' #probabilities 0.2, 0.4 and 0.4.
#' pop_germc( init_sb="initialSB", germ=c(0.2,0.4,0.4))

#'rm(mf, dfgenotype, xprobab) 

pop_germc <-
function(init_sb, germ, max_vec_length=1e+07){
cat("pop_germc starts...")
if(sum(germ)>1|sum(germ)<0){stop("pop_germc(germ): germ follows 0 <= germ <= 1.")}
dfgenotype <- get0("dfgenotype", envir = parent.frame(n = 1))


### --- name flexibel
first_amount <- eval(parse(text=paste("dfgenotype$",init_sb,sep="")))       # the amount that gives the start

#####----- structure for all cohorts
germ_tmp <- data.frame(matrix(nrow=nrow(dfgenotype), ncol=length(germ)+1))        #for saving the germination result
germ_tmp[is.na(germ_tmp)] <-0
names(germ_tmp)[-length(germ_tmp)] <- paste("germ",c(1:(length(germ))),sep="")
names(germ_tmp)[length(germ_tmp)] <- "germ_dorm"                                    #the last is named "dorm"



for(pres_gt in which(first_amount > 0)){                                  #just work with the genotypes that are present
    i1 <- first_amount[pres_gt] %/% max_vec_length
    i2 <- first_amount[pres_gt] %% max_vec_length
    if(i1 > 0){     #for every i1 chunk the same like the i2 rest (see below)
      for(i in 1:i1){
            tmp <- 0
            for(cohort in seq_along(germ)){
                       tmp0 <- rbinom(1, max_vec_length-tmp, prob = germ[cohort])
                       germ_tmp[,cohort][pres_gt] <- germ_tmp[,cohort][pres_gt] + tmp0
                       tmp <- tmp + tmp0                                        #every chunck of Max_cev is devidet
                            }#END for(cohort)
           germ_tmp[,length(germ)][pres_gt] <- germ_tmp[,cohort][pres_gt] + max_vec_length - tmp# the last cohort is just a difference
           
      }#END for(i)
    }#END if(i1)
  ###----- i2 part
    tmp <- 0
    for(cohort in seq_along(germ)){
            tmp0 <- rbinom(1, i2-tmp, prob=germ[cohort])
            germ_tmp[,cohort][pres_gt] <- germ_tmp[,cohort][pres_gt] + tmp0
            tmp <- tmp + tmp0                                                   #every chunch of Max_cev is devidet
            }#END for(cohort)
    germ_tmp[,length(germ)+1][pres_gt] <- germ_tmp[,length(germ)+1][pres_gt] + i2 - tmp# the last cohort is just a difference
}#END for(pres_gt)
dfgenotype <- data.frame(dfgenotype,germ_tmp)
cat("finished!\n")
assign("dfgenotype", value=dfgenotype, pos = -1, envir=parent.frame(n = 1))            
return(dfgenotype)     
}
