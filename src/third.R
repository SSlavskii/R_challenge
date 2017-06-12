library(cumstats)
load("fdr.RData")

fdr <- function(labels, scores) {
  # create a dataframe with 2 columns: value of FDR and corr threshold
  labels <- labels[order(scores, decreasing = T)]
  data.frame(FDR = (1 - cummean(labels)) * 100, value = sort(scores, decreasing = T))
}

get_threshold <- function(df, percentage) {
  min(subset(df, FDR < percentage)$value) 
}

start <- Sys.time()
# converst type to boolean variable
fdr.df$label <- fdr.df$type == "control"
test <- fdr(fdr.df$label, fdr.df$value)
get_threshold(test, 5)

Sys.time() - start

ggplot(fdr.df, aes(x = value, fill = type)) + 
  geom_histogram(aes(y=..density..), alpha = 0.5, position="identity", col="black") + 
  theme_bw()
