#summarise susceptible to infection for expert elicitation
susceptible_plots <- projections_full %>% 
    split(~vaccine_type) %>%
    map(summarise_susceptible)
ggsave(
   filename = "plots/susceptibles.pdf", 
   plot = gridExtra::marrangeGrob(susceptible_plots, nrow=1, ncol=1), 
   width = 9, height = 9
)
susceptible_plots <- projections_children %>% 
    split(~vaccine_type) %>%
    map(summarise_susceptible)
ggsave(
   filename = "plots/susceptibles_children.pdf", 
   plot = gridExtra::marrangeGrob(susceptible_plots, nrow=1, ncol=1), 
   width = 9, height = 9
)

#csv file (do weekly)
projections_full %>% 
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
        susceptible_to_infection = round(infection, 3),
        susceptible_to_disease = round(disease, 3)
    ) %>%
    write_csv("data/output/immunity_projections.csv")

projections_children %>% 
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
        susceptible_to_infection = round(infection, 3),
        susceptible_to_disease = round(disease, 3)
    ) %>%
    write_csv("data/output/immunity_projections_children.csv")
#now format S and V for francesco

#What just want S and V? or prop immune

#Reff Calculations using Somaliland Digaale Contact Matrix and R0 estimates
digaale <- socialmixr::get_survey("https://zenodo.org/doi/10.5281/zenodo.5226280")
digaale_pop <- read.csv("https://zenodo.org/records/7071876/files/espicc_somaliland_digaale_survey_population.csv?download=1")
contact_data <- socialmixr::contact_matrix(
    digaale,
    survey.pop = digaale_pop,
    age.limits = digaale_pop$lower.age.limit,
    symmetric = TRUE,
    weigh.dayofweek = TRUE,
    weights = "sample_weight",
    per.capita = TRUE
)

#calculate proportion immune to infection in each age group of the contact matrix (at the start of projections)
contact_data$demography$age.group
contact_age_groups <- c(seq(10, 60, 10), 100)
age_group_map <- map_int(output_age_group_end_adults, ~min(which(.x <= contact_age_groups)))

immunity <- projections_full %>% 
    filter(date == min(date)) %>%
    mutate(age_group = age_group_map[as.numeric(age_group)]) %>%
    group_by(scenario, vaccine_type, age_group) %>%
    summarise(
        immunity = sum(immune)/sum(population),
        .groups = "drop"
    ) %>%
    arrange(age_group) %>% 
    split(.$scenario) %>%
    map(~split(.x, .$vaccine_type) %>% map(function(x){x$immunity}))

R0 <- read_csv("data/raw/infection_parameters.csv") %>%
    filter(parameter %in% c("r0_max", "r0_min")) %>%
    select(disease, parameter, value_gen) %>%
    pivot_wider(names_from = parameter, values_from = value_gen) %>%
    rowwise() %>%
    transmute(
        disease = disease,
        R0 = list(list(min = r0_min, max = r0_max))
    ) %>%
    pull(R0, disease)

R0$`Hib disease` <- list(min = 3, max = 3.5)  #https://karger.com/neo/article-abstract/50/2/114/367324/Colonisation-of-Haemophilus-influenzae-and?redirectedFrom=PDF
R0$`pneumococcal disease` <- list(min = 3, max = 3.5)  #https://karger.com/neo/article-abstract/50/2/114/367324/Colonisation-of-Haemophilus-influenzae-and?redirectedFrom=PDF
R0$rotavirus <- list(min = 28, max = 32) #https://elischolar.library.yale.edu/cgi/viewcontent.cgi?article=1071&context=ysphtdl

Reff <- map(immunity, function(scenario) {
    imap(scenario, function(disease, disease_name) {
        prop_immune <- epimixr::adjust_immunity(
            contact_data$matrix,
            disease
        )
        list(
            min = R0[[disease_name]]$min * (1 - prop_immune),
            max = R0[[disease_name]]$max * (1 - prop_immune)
        )
    })
})
saveRDS(Reff, "data/output/Reff_at_start.rds")
