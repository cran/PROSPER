
#setwd("D:/prosp/")



#file.create("D:\\prosp\\output.txt")
#fileConn<-file("output.txt")


create_prosper <- function(fname=NA, path=NA){

#--- create the file ------
path <- ifelse(is.na(path), getwd(), path)
file_name <- ifelse(is.na(fname), "new_prosper", fname)

file.create("D:\\prosp\\output.txt")
file_name <- paste(file_name, ".txt", sep="")
file.create(paste(path,"/", file_name, sep=""))


file_name <- paste(path,"/",file_name, sep="")

#--- create the text parts ------
start_text <- c("
prosper.LOLRI <- function(
      param.weed     = PROSPER::param.LOLRI,
      area           = 100,
      af             = c(0.005, 0.01),
      dom            = 0.5,
      epis           = 0,
      put            = 0.05,
      sdrate         = 0.4,
      thresh         = 20,
      Rmx            = 10,
      dc             = 150,
      kc             = 1/11,
      kw             = 1/33,
      SSmax          = 30000,
      rate           = 100,             
      pen_co         = c(1,0.5),
      duration       = 10,
      repetitions    = 1,
      crop_list      = c(\"wheat\"),
      max_vec_length = 1e+07)
{\n
")

prep <- c("
#-----Simulation preparations --------------------------------------------------
#things that you do only want to do ONCE in your simulation run
cat(\"SIMULATION START\")

#### --- Controll of the given data ----
crop.list <- NA # not ideal
mod_check(param.weed=param.weed, area=area, duration=duration, repetitions=repetitions, crop_list=crop_list, max_vec_length=max_vec_length)

#- generating necessary objects: dfgenotype, mf, xprobab -----------------------
struc_preparation2(Rmx=Rmx, af=af, epis=epis, dom=dom)
cat(\"--------preparation ready------\")\n
")

repet_prep <- c("
for(rep_counter in seq(repetitions)){
cat(\"---------------------------------------------\")
cat(\"starting simulation run\t\", rep_counter)
crop <- crop.list[1]
year <- 0

# -- the starting weed density per m^2
weed_den          <- quanti(step_name=\"den0\", crop=crop)   #getting the density
gen_freq(result=\"SB_autumn\", af=af, n_seeds=weed_den*area) 

#Here come the things that you want to do every year in duration
cat(\"---------------------------------------------\")\n
")

start_time_loop <- c("
cat(\"---------------------------------------------\")
cat(\"starting the timeloop\")
repeat{
year <- year + 1
crop <- crop.list[year]
cat(\"time step \", year, \"starts\")\n
")

quanti_step <- c("
### winter reduces the SB_autumn  -----------------------------------------------
# and forms the SB_spring
dieprob     <- quanti(step_name=\"Pbetweebdeath\", crop=crop, area=area, res_max=1, res_min=0)
pop_step(start=\"SB_autumn\", result=\"SB_spring\",surv_prob=1-dieprob)\n
")


quanti_germc <- c("
### germination of cohorts -----------------------------------------------------
# cohorts should be defined here
# calculates the probability of germination in the first and second cohort
germ1 <- quanti(crop=crop, area=area, step_name=\"germ1\")  #calculates the probability of germination in the
germ2 <- quanti(crop=crop, area=area, step_name=\"germ2\")  #calculates the probability of germination in the second cohort
pop_germc(init_sb=\"SB_spring\", germ=c(germ1,germ2))  #defines which individual goes which way\n
")


herb <- c("
### surviving of herbicide -----------------------------------------------------                        
sel_herb(start=\"surv01\", result=\"survherb1\", thresh=thresh, sdrate=sdrate, rate=rate, put=put)
sel_herb(start=\"germ2\", result=\"survherb2\", thresh=thresh, sdrate=sdrate, rate=rate, put=put)\n                        
")


reprod <- c("
### production of new seeds ----------------------------------------------------
pop_reprod(start=c(\"survherb1\",\"survherb2\"), area=area, kw=kw, pen_co=pen_co,
                             kc=kc, dc=dc, crop_inr=which(crop_list==crop), SSmax=SSmax, yield=FALSE)\n 
")
                             

diploid <- c("
# -- Genotypes of new seeds  -------
gen_diploid(start=c(\"survherb1\", \"survherb2\"), start_comb=\"weedatharvest\", result=\"new_seeds\")\n
")                             

saveSimdata <- c("
struc_saveSimData(rep_counter=rep_counter, simcycle=year, 
                               start_name=\"SB_autumn\", end_name=\"SB_autumn_end\",
                                simstruc=c(repetitions, duration))\n
")

end_text <- c("
#-------------------------------------------------------------------------------
### some rules for ending the repetition.
simend <- struc_endSim(year,\"SB_autumn\")
if(simend==TRUE) {break}
}#END repeat (TIMELOOP)
}#END for(rep_counter)
invisible(get0(\"sim_result\"))
}\n     
")


#--- writing the text parts to the file ----

cat(start_text, file=file_name)
cat(prep, file=file_name, append=TRUE)
cat(repet_prep, file=file_name, append=TRUE)
cat(start_time_loop, file=file_name, append=TRUE)
cat(quanti_step, file=file_name, append=TRUE)
cat(quanti_germc, file=file_name, append=TRUE)
cat(herb, file=file_name, append=TRUE)
cat(reprod, file=file_name, append=TRUE)
cat(diploid, file=file_name, append=TRUE)
cat(saveSimdata, file=file_name, append=TRUE)
cat(end_text, file=file_name, append=TRUE)


#close(fileConn)
}#END function
