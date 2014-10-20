library("ggplot2")
toSeconds <- function(x){
    # convert time from H:M:S or HH:MM:SS to seconds
    if (!is.character(x)) stop("x must be a character string of the form H:M:S")
    if (length(x)<=0)return(x)
    
    unlist(
        lapply(x,
               function(i){
                   i <- as.numeric(strsplit(i,':',fixed=TRUE)[[1]])
                   if (length(i) == 3) 
                       i[1]*3600 + i[2]*60 + i[3]
                   else if (length(i) == 2) 
                       i[1]*60 + i[2]
                   else if (length(i) == 1) 
                       i[1]
               }  
        )  
    )  
} 
# testing conditions: 2 CPUs on SLURM, as much memory as they need so SLURM doesn't kill the job
setwd("C:/Users/User1/Dropbox/Research/Huttenhower/analysis")
# taking the max of the time for human RNA db and Silva DB
rna_stats_raw <- read.csv("rna_stats.csv", header=TRUE)
times <- toSeconds(as.character(rna_stats_raw$CPUTime))
print(times)
# memory was originally in KB. Converting to GB here
# converting time from seconds to minutes
plotData <- data.frame(Tool=rna_stats_raw$Tool, CPUTime=times/60, MaxRSS=rna_stats_raw$MaxRSS/1000000, MaxVMSize=rna_stats_raw$MaxVMSize/1000000)

# don't include the stuff that was run without trimming (optional)
plotData = subset(plotData, Tool == 'bowtie2-trim' | Tool == 'bwa-trim' | Tool == 'bmtagger')

p_runtime <- qplot(data=plotData, x=Tool, y=CPUTime, geom="boxplot", ylab="CPUTime (min)", xlab="Tool", color=Tool) + 
    ggtitle("Running time, Taking into Account Parallelization")
print(p_runtime)

p_maxrss <- qplot(data=plotData, x=Tool, y=MaxRSS, color=Tool, geom="boxplot", xlab="Tool", ylab="Max Memory Usage (GB)") +
    ggtitle("Max Memory Usage, Taking into Account Parallelization")
print(p_maxrss)

p_maxvmsize <- qplot(data=plotData, x=Tool, y=MaxVMSize, color=Tool, geom="boxplot", xlab="Tool", ylab="Max Virtual Memory Usage (GB)") +
    ggtitle("Max Virtual Memory Usage, Taking Into Account Parallelization")
print(p_maxvmsize)

# including the human RNA db and Silva DB separately
rna_stats_expanded <- read.csv("rna_stats_expanded.csv", header=TRUE)
times_exp <- toSeconds(as.character(rna_stats_expanded$CPUTime))
print(times_exp)
plotData_exp <- data.frame(Tool=rna_stats_expanded$Tool, CPUTime=times_exp/60, MaxRSS=rna_stats_expanded$MaxRSS/1000000, MaxVMSize=rna_stats_expanded$MaxVMSize/1000000)

# don't include the stuff that was run without trimming (optional)
plotData_exp = subset(plotData_exp, Tool != 'bowtie2-notrim-human' & Tool != 'bowtie2-notrim-silva' & Tool != 'bwa-notrim-human' & Tool != 'bwa-notrim-silva')

p_runtime_exp <- qplot(data=plotData_exp, x=Tool, y=CPUTime, color=Tool, geom="boxplot", xlab="Tool", ylab="CPUTime (min)") +
    ggtitle("Running time")
# make the labels vertical
p_runtime_exp <- p_runtime_exp + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_runtime_exp)

p_maxrss_exp <- qplot(data=plotData_exp, x=Tool, y=MaxRSS, color=Tool, geom="boxplot", xlab="Tool", ylab="Max Memory Usage (GB)") + 
    ggtitle("Max Memory Usage")
p_maxrss_exp <- p_maxrss_exp + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_maxrss_exp)

p_maxvmsize_exp <- qplot(data=plotData_exp, x=Tool, y=MaxVMSize, color=Tool, geom="boxplot", xlab="Tool", ylab="Max Virtual Memory Usage (GB)") +
    ggtitle("Max Virtual Memory Usage")
p_maxvmsize_exp <- p_maxvmsize_exp + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_maxvmsize_exp)