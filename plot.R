###### plot ###


## 每輪的點圖

plot_result <-function(x, title=NULL)
{
  df <- data.frame(result = unlist(x))
  
  ggplot(df, aes(y = result, x = c(1:length(result)))) +
    geom_point() +
    geom_hline(yintercept = mean(df$result), linetype = "dashed", color = 'red') +
    labs(title = paste("Empirical Likelihood Simulation", title), x = 'outer iterations') +
    scale_x_continuous(breaks = seq(0, length(df$c), 1))
  
}