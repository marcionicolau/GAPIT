`GAPIT.Memory` <- function(Memory =NULL,Infor) {
  #Object: To report memory usage
  #Output: Memory 
  #Authors: Zhiwu Zhang
  # Last update: June 6, 2011 
  ##############################################################################################
  gc()
  size <- memory.size()
  #print(paste("Memory usage: ",size," for", Infor))
  if(is.null(Memory)) {
    Increased=0
    Memory =cbind(Infor,size ,Increased)
  }else{
    Increased=0
    Memory.current=cbind(Infor,size ,Increased)
    Memory=rbind(Memory,Memory.current)
    Memory[nrow(Memory),3]=
      as.numeric(as.matrix(Memory[nrow(Memory),2]))-as.numeric(as.matrix(Memory[nrow(Memory)-1,2]))
  }    
  return (Memory)
}#end of GAPIT.Memory function

`GAPIT.Timmer` <- function(Timmer=NULL,Infor) {
  #Object: To report current time
  #Output: Timmer
  #Authors: Zhiwu Zhang
  # Last update: may 8, 2011 
  ##############################################################################################
  
  Time<-Sys.time()
  if(is.null(Timmer)) {
    Elapsed=0
    Timmer=cbind(Infor,Time,Elapsed)
  }else{
    Elapsed=0
    Timmer.current=cbind(Infor,Time,Elapsed)
    Timmer=rbind(Timmer,Timmer.current)
    Timmer[nrow(Timmer),3]=
      as.numeric(as.matrix(Timmer[nrow(Timmer),2]))-as.numeric(as.matrix(Timmer[nrow(Timmer)-1,2]))
  }
  
  #print(paste('Time used: ', Timmer[nrow(Timmer),3], ' seconds for ',Infor,sep="" )) 
  return (Timmer)
}#end of GAPIT.Timmer function

`GAPIT.RemoveDuplicate` <- function(Y) {
  #Object: NA
  #Output: NA
  #Authors: Zhiwu Zhang 
  # Last update: Augus 30, 2011 
  ##############################################################################################
  return (Y[match(unique(Y[,1]), Y[,1], nomatch = 0), ] )
}

`GAPIT.replaceNaN` <- function(LL) {
  #handler of grids with NaN log
  #Authors: Zhiwu Zhang
  # Last update: may 12, 2011 
  ##############################################################################################
  
  #handler of grids with NaN log 
  index <- (LL == "NaN")
  if(length(index) > 0) theMin <- min(LL[!index])
  if(length(index) < 1) theMin <- "NaN"
  LL[index] <- theMin
  return(LL)    
}

`GAPIT.Memory.Object` <- function(name.of.trait="Trait"){
  # Object: To report memoery usage
  # Authors: Heuristic Andrew
  # http://heuristically.wordpress.com/2010/01/04/r-memory-usage-statistics-variable/
  # Modified by Zhiwu Zhang
  # Last update: may 29, 2011 
  ##############################################################################################
  
  # print aggregate memory usage statistics 
  print(paste('R is using', memory.size(), 'MB out of limit', memory.limit(), 'MB')) 
  
  # create function to return matrix of memory consumption 
  object.sizes <- function() { 
    return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name) 
      object.size(get(object.name)))))) 
  } 
  
  # export file in table format 
  memory=object.sizes() 
  file=paste("GAPIT.", name.of.trait,".Memory.Object.csv" ,sep = "")
  write.table(memory, file, quote = FALSE, sep = ",", row.names = TRUE,col.names = TRUE)
  
  # export file in PDF format 
  pdf(paste("GAPIT.", name.of.trait,".Memory.Object.pdf" ,sep = ""))
  # draw bar plot 
  barplot(object.sizes(), 
          main="Memory usage by object", ylab="Bytes", xlab="Variable name", 
          col=heat.colors(length(object.sizes()))) 
  # draw dot chart 
  dotchart(object.sizes(), main="Memory usage by object", xlab="Bytes") 
  # draw pie chart 
  pie(object.sizes(), main="Memory usage by object")
  dev.off()  
}

`GAPIT.0000` <- function() {
  ################################################################################
  #GAPIT: Genome Association and Prediction Integrated Tool
  #This is an R package that performs Genome Wide Association Study (GWAS) and 
  # genome prediction (or genomic selection). #This program uses state-of-the-art methods 
  #developed for statistical genetics, such as the unified mixed model, EMMA, 
  #the compressed mixed linear model, and P3D/EMMAx.
  
  #Designed by Zhiwu Zhang
  #Writen by Alex Lipka, Feng Tian and Zhiwu Zhang
  GAPIT.Version="2.22, December 7, 2012 (GS indicator: 1-Phe, 1.5-noPhe&GRPPhe, 2-rest)"
  return(GAPIT.Version)
}
