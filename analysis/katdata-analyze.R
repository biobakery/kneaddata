library("ggplot2")

setwd("C:/Users/User1/Dropbox/Research/Huttenhower/analysis")
datafile = "katdata.csv"
dataset <- read.csv(datafile, header=T)
totalSensitivity <- (dataset$MergedHuman + dataset$MergedSilva)/dataset$MergedOut
totalSpecificity <- 1 - (dataset$MergedOut - dataset$MergedHuman - dataset$MergedSilva) / (dataset$OrigTotal - dataset$OrigHuman + dataset$OrigSilva)
print(totalSensitivity)
print(totalSpecificity)
humanSensitivity <- dataset$HumanHuman/dataset$OrigHuman
humanSpecificity <- 1 - (dataset$HumanTotal - dataset$HumanHuman)/(dataset$OrigTotal - dataset$OrigHuman)

silvaSensitivity <- dataset$SilvaSilva/dataset$OrigSilva
silvaSpecificity <- 1 - (dataset$SilvaTotal - dataset$SilvaSilva)/(dataset$OrigTotal - dataset$OrigSilva)

print(humanSensitivity)
print(humanSpecificity)
print(silvaSensitivity)
print(silvaSpecificity)
plotData <- data.frame(Aligner=rep(dataset$Aligner,6), measurement=c(totalSensitivity, totalSpecificity, humanSensitivity, humanSpecificity, silvaSensitivity, silvaSpecificity), type=c(rep("Total Sensitivity", 5), rep("Total Specificity", 5), rep("Human Sensitivity", 5), rep("Human Specificity", 5), rep("rRNA Sensitivity",5), rep("rRNA Specificity", 5)))

# don't include the stuff that was run without trimming (optional)
plotData = subset(plotData, Aligner == 'bowtie2' | Aligner == 'bwa' | Aligner == 'bmtagger')
print(plotData)

# plot the performance metrics
ggplot(data = plotData, aes(x = Aligner, y = measurement, fill = Aligner)) + geom_bar(stat="identity") + facet_wrap(~type) + scale_y_continuous(limits=c(0.975, 1), oob = rescale_none) + 
    ggtitle("Accuracy on Transcriptomic Data") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot zoomed in version
ggplot(data = plotData, aes(x = Aligner, y = measurement, fill = Aligner)) + geom_bar(stat="identity") + facet_wrap(~type) + scale_y_continuous(limits=c(0.999, 1), oob = rescale_none) + 
    ggtitle("Accuracy on Transcriptomic Data, Zoomed In") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))