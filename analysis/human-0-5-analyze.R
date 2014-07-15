library("ggplot2")
data <- read.csv("human-0-5.csv")

# Calculate specificity and sensitivity
totalHuman <- data$HumanTrimmedPE + data$HumanTrimmedSE1 + data$HumanTrimmedSE2
removedHuman <- data$HumanRemovedPE + data$HumanRemovedSE1 + data$HumanRemovedSE2
sensitivity <- removedHuman / totalHuman 
sensitivity
sensitivity[1] <- 1
sensitivity

prop.contam <- data$HumanOrig/data$NumReadsOrig
prop.contam

totalTrimmedReads <- data$NumReadsTrimmedPE + data$NumReadsTrimmedSE1 +
data$NumReadsTrimmedSE2

totalNonHuman <- totalTrimmedReads - totalHuman
removedNonHuman <- data$NonHumanRemovedPE + data$NonHumanRemovedSE1 + 
    data$NonHumanRemovedSE2
specificity <- 1 - (removedNonHuman/totalNonHuman)

plot_data <- data.frame(prop.contam, specificity, sensitivity)
plt <- ggplot(plot_data, aes(x = prop.contam)) + 
    geom_point(aes(y = sensitivity, color="P(human|human)")) +
    geom_point(aes(y = specificity, color="P(nonhuman|nonhuman)")) +
    geom_line(aes(y = sensitivity, color="P(human|human)")) +
    geom_line(aes(y = specificity, color="P(nonhuman|nonhuman)"))

plt <- plt + xlab("Human Reads / Total Reads") + ylab("Proportion") + 
    theme(legend.title=element_blank())
plt <- plt + 
    ggtitle("Specificity and Sensitivity\nas a Function of Human Contamination")
plt
