# Be reproducible queens
set.seed(1) 

# Functions to generate data (surrogates non-overlapping)
gen.data = function(setting, n1, n0) {
  s1 = g.1(n1)
  y1 = f.cond.1(s1)
  s0 = g.0(n0)
  y0 = f.cond.0(s0)
  return(data.frame("s1" = s1, "y1" = y1, "s0" = s0, "y0" = y0))
}
f.cond.1 = function(s.vector) {
  eps1 = rnorm(length(s.vector),0,3)
  y1 = 2+5*s.vector+1 + 1*s.vector + eps1
  return(y1)		
}
f.cond.0 = function(s.vector) {
  eps0 = rnorm(length(s.vector),0,3)
  y0 = 2+5*s.vector+ eps0
  return(y0)		
}
g.1 = function(n, alpha0=5) { return(rnorm(n, alpha0 + 1,2))}
g.0 = function(n, alpha0=5) { return(rnorm(n, alpha0,1))}
## Simulate data (10,000 patients per group for demonstration)
data_overlap = gen.data(n1=1000, n0=1000) 

# Functions to generate data (surrogates non-overlapping)
## Only g.1 changes 
g.1 = function(n, alpha0=5) { return(rnorm(n, alpha0 + 1,1/2))}
## Simulate data (10,000 patients per group for demonstration)
data_nonoverlap = gen.data(n1=1000, n0=1000) 

# Load packages
cols = c("#787ff6", "#ffbd59", "#8bdddb", "#ff99ff",  "#7dd5f6") ## color palette

# Plot histograms of surrogates in control/treatment groups 
data_overlap |> 
  mutate(Overlap = "Surrogate Distributions Overlap") |> 
  bind_rows(
    data_nonoverlap |> 
      mutate(Overlap = "Surrogate Distributions Do Not Overlap")
    ) |> 
  mutate(Overlap = factor(x = Overlap, 
                          levels = c("Surrogate Distributions Overlap", 
                                     "Surrogate Distributions Do Not Overlap"))) |> 
  ggplot() + 
  geom_histogram(aes(x = s0, fill = "Control"), 
                 alpha = 0.7) + 
  geom_histogram(aes(x = s1, fill = "Treatment"), 
                 alpha = 0.5) + 
  facet_wrap(~Overlap) + 
  labs(x = "Surrogate Marker", 
       y = "Number of Patients") + 
  scale_fill_manual(
    values = c("Control" = cols[2], "Treatment" = cols[1]), # Specify colors for each level
    labels = c("Control" = "Control Group", "Treatment" = "Treatment Group"), # Specify custom labels
    name = ""
  ) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white"))

ggsave(filename = "figures/surrogate_density_overlap.pdf", 
       device = "pdf", width = 7, height = 5)
