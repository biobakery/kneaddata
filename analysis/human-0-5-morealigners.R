library("ggplot2")

# load data from hand-curated csv file
data <- read.csv("human-0-5-morealigners-aggregated.csv", header=T)

data <- subset(data, Aligner != 'bwa-notrim' & Aligner != 'bowtie2-notrim')
#data <- subset(data, Aligner != 'bwa-notrim')
#data <- subset(data, Aligner != 'bowtie2-notrim')
data
len_data = length(data$TrimHuman)
sensitivity <- rep(NA, len_data)
for (i in 1:len_data) {
    if (data$OutHuman[i] == 0 && data$TrimHuman[i] == 0) {
        sensitivity[i] <- 1
    }
    else {
        sensitivity[i] <- data$OutHuman[i] / data$TrimHuman[i]
    }
}
sensitivity
prop.contam <- data$OrigHuman / data$OrigTotal
prop.contam

totalNonHuman <- data$TrimTotal - data$TrimHuman
removedNonHuman <- data$OutTotal - data$OutHuman
specificity <- 1 - (removedNonHuman/totalNonHuman)
specificity

plot_data <- data.frame(Dataset=data$Dataset, Aligner=data$Aligner, prop.contam,
                        sensitivity, specificity)

plt <- ggplot(plot_data, aes(x = prop.contam)) +
    geom_point(aes(y = sensitivity, color = Aligner, shape = "P(human|human)")) + 
    geom_point(aes(y = specificity, color = Aligner, shape =
                   "P(nonhuman|nonhuman)")) + 
    geom_line(aes(y = sensitivity, color = Aligner, linetype = "P(human|human)")) + 
    geom_line(aes(y = specificity, color = Aligner, linetype =
                   "P(nonhuman|nonhuman)"))
plt <- plt + xlab("Human Reads / Total Reads") + ylab("Proportion") + 
    theme(legend.title=element_blank())
plt <- plt + 
    ggtitle("Specificity and Sensitivity\nas a Function of Human Contamination")
plt
