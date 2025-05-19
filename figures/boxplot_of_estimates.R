cols = c("#787ff6", "#ffbd59", "#8bdddb", "#ff99ff",  "#7dd5f6") ## color palette
boxplot_of_estimates = function(data) {
  data |> 
    ggplot(aes(x = method, y = est, fill = method)) + 
    geom_hline(aes(yintercept = truth), linetype = 2, color = "black") + 
    geom_boxplot() + 
    xlab("Method") + 
    ylab("Estimate") + 
    facet_grid(cols = vars(parametric), 
               rows = vars(quantity), 
               scales = "free", 
               labeller = labeller(parametric = label_value, 
                                   quantity = label_parsed)) + 
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
    scale_fill_manual(name = "Method:", values = cols, guide = "none") + 
    theme_minimal() + 
    theme(legend.position = "top", 
          strip.background = element_rect(fill = "black"), 
          strip.text = element_text(color = "white"))
}