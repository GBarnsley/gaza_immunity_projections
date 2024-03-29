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

#bar graph susceptibility

scenario_cols <- c("Escalation" = palette_periods$escalation, "Status Quo" = palette_periods$`status quo`, "Ceasefire" = palette_periods$ceasefire)

plot_df <- projections_children %>%
    filter(age_group <= "6th") %>%
    mutate(col = "Children under 59mo") %>%
    rbind(mutate(projections_full, col = "All ages")) %>%
    group_by(col, scenario, vaccine_type, date) %>%
    summarise(
        susceptibility = 1 - (sum(immune)/sum(population)),
        susceptibility_disease = 1 - (sum(immune_disease)/sum(population)),
        .groups = "drop_last"
    ) %>%
    summarise(
        susceptibility = mean(susceptibility),
        susceptibility_disease = mean(susceptibility_disease),
        .groups = "drop"
    ) %>%
    mutate(
        scenario = factor(scenario, levels = c("ceasefire", "status_quo", "escalation")),
        scenario = fct_recode(scenario, "Escalation" = "escalation", "Status Quo" = "status_quo", "Ceasefire" = "ceasefire"),
        col = factor(col, levels = c("Children under 59mo", "All ages")),
        vaccine_type = str_replace(str_to_title(vaccine_type), " ", "\n"),
        vaccine_type = case_when(
            vaccine_type == "Polio\n- Vaccine-Derived" ~ "Polio\n(Vaccine)",
            vaccine_type == "Polio\n- Wildtype" ~ "Polio\n(Wildtype)",
            TRUE ~ vaccine_type
        )
    )

col_plot <- plot_df %>%
    pivot_longer(c(susceptibility, susceptibility_disease), names_to = "type", values_to = "value") %>%
    mutate(
        type = if_else(type == "susceptibility", "Against Infection", "Against Disease"),
        type = factor(type, levels = c("Against Infection", "Against Disease"))
    ) %>%
    ggplot(aes(x = vaccine_type, y = value, fill = scenario)) +
        geom_col(position = "dodge", alpha = 0.5) +
        ggpubr::theme_pubr() +
        facet_grid(rows = vars(type), cols = vars(col)) +
        labs(fill = "", x = "", y = "Susceptibility (%)", title = paste0("Average estimated susceptibility to vaccine preventable diseases")) +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values = scenario_cols)
        
ggsave(
    "plots/output_susceptibility_cols.png", col_plot, height = 10, width = 17
)

#prewar and current susceptibility (hib, pcv and rota)#
pre_war_output_dates <- seq(t_record_from + date_start, date_projection_start, by = "1 month")
pre_war_output_dates <- pre_war_output_dates[-length(pre_war_output_dates)]
pre_war_outputs <- projections_full %>%
    filter(
        date %in% pre_war_output_dates
    ) %>% 
    mutate(
        susceptible = 1 - (immune_disease / population),
        immune_disease_only = immune_disease / population,
        immune = immune / population,
    ) %>%
    convert_to_output_format_compartments(c("susceptible", "immune_disease_only", "immune")) %>%
    map(~list_transpose(.x)) %>%
    list_transpose()

pre_war_outputs <- map(
    pre_war_outputs,
    ~map(.x, function(x) {
        map(x, function(y) {
            rownames(y) <- c("m-1", paste0("intervening", 1:(length(rownames(y)) - 1)))
            y
        })
    })
)
saveRDS(pre_war_outputs$susceptible, "data/output/output_pre_war_susceptible.rds")
saveRDS(pre_war_outputs$immune_disease_only, "data/output/output_pre_war_immune_disease_only.rds")
saveRDS(pre_war_outputs$immune, "data/output/output_pre_war_immune.rds")