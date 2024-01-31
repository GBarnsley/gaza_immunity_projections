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

#output full results for the start of each month of simulation for francesco's simulations
output_dates <- accumulate(1:projection_months, function(x, y) {
    if(is.na(x)) {
        date_projection_start
    } else {
        month(x) <- month(x) + 1
        x
    }
}, .init = NA)[-1] #we're just ticking up one month each time

final_outputs <- projections_full %>%
    filter(
        date %in% output_dates
    ) %>% 
    mutate(
        susceptible = 1 - (immune_disease / population),
        immune_disease_only = immune_disease / population,
        immune = immune / population,
    ) %>%
    convert_to_output_format_compartments(c("susceptible", "immune_disease_only", "immune")) %>%
    map(~list_transpose(.x)) %>%
    list_transpose()

saveRDS(final_outputs$susceptible, "data/output/output_susceptible.rds")
saveRDS(final_outputs$immune_disease_only, "data/output/output_immune_disease_only.rds")
saveRDS(final_outputs$immune, "data/output/output_immune.rds")

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

R0 <- read_csv("data/raw/infection_parameters.csv", show_col_types = FALSE) %>%
    filter(parameter %in% c("r0_max", "r0_min")) %>%
    select(disease, parameter, value_gen) %>%
    pivot_wider(names_from = parameter, values_from = value_gen) %>%
    rowwise() %>%
    transmute(
        vaccine_type = disease,
        R0_min = r0_min,
        R0_max = r0_max
    ) %>% 
    rbind(
        tribble(
            ~vaccine_type, ~R0_min, ~R0_max,
            "Hib disease", 3, 3.5, #https://karger.com/neo/article-abstract/50/2/114/367324/Colonisation-of-Haemophilus-influenzae-and?redirectedFrom=PDF
            "pneumococcal disease", 3, 3.5, #https://karger.com/neo/article-abstract/50/2/114/367324/Colonisation-of-Haemophilus-influenzae-and?redirectedFrom=PDF
            "rotavirus", 28, 32 #https://elischolar.library.yale.edu/cgi/viewcontent.cgi?article=1071&context=ysphtdl
        )
    )

Reff <- projections_full %>% 
    filter(date %in% output_dates) %>%
    mutate(age_group = age_group_map[as.numeric(age_group)]) %>%
    group_by(scenario, vaccine_type, date, age_group) %>%
    summarise(
        immunity = sum(immune)/sum(population),
        .groups = "drop"
    ) %>%
    arrange(age_group) %>%
    left_join(R0, by = "vaccine_type") %>%
    group_by(scenario, vaccine_type, date) %>%
    summarise(
        Reff = list({
            prop_immune <- epimixr::adjust_immunity(
                contact_data$matrix,
                immunity
            )
            R0s <- sample_uniform(N, unique(R0_min), unique(R0_max))
            ecdf(R0s)
        }),
        .groups = "drop"
    ) %>%
    convert_to_output_format_Reff()

saveRDS(Reff, "data/output/output_Reff.rds")
