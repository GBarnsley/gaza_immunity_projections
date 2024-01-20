#summarise susceptible to infection for expert elicitation
susceptible_plots <- projections %>% 
    split(~vaccine_type) %>%
    map(summarise_susceptible)
ggsave(
   filename = "plots/susceptibles.pdf", 
   plot = gridExtra::marrangeGrob(susceptible_plots, nrow=1, ncol=1), 
   width = 9, height = 9
)
#csv file (do weekly)
projections %>% 
    mutate(week = floor_date(date, "week")) %>%
    group_by(scenario, vaccine_type, week) %>%
    summarise(
        date = median(date),
        infection = mean(1 - (sum(immune)/sum(population))),
        disease = mean(1 - (sum(immune_disease)/sum(population))),
        .groups = "drop"
    ) %>% 
    transmute(
        scenario = scenario,
        disease_type = vaccine_type,
        date = date,
        immune_to_infection = round(infection, 3),
        immune_to_disease = round(disease, 3)
    ) %>%
    write_csv("data/output/immunity_projections.csv")
#now format S and V for francesco

#What just want S and V? or prop immune

#Reff Calculations (use kevins contact matrix)
#?epimixr::adjust_immunity
