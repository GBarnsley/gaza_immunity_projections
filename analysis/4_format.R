#summarise susceptible to infection for expert elicitation
susceptible_plots <- projections %>% 
    split(~vaccine_type) %>%
    map(summarise_susceptible)
ggsave(
   filename = "plots/susceptibles.pdf", 
   plot = gridExtra::marrangeGrob(susceptible_plots, nrow=1, ncol=1), 
   width = 9, height = 9
)
#now format S and V for francesco

#What just want S and V? or prop immune

#Reff Calculations (use kevins contact matrix)
#?epimixr::adjust_immunity
