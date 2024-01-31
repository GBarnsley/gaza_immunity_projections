#What type of model
type <- "static_model"

t_projection_starts <- as.numeric(date_projection_start - date_start)
t_projection_ends <- as.numeric(date_projection_end - date_start)

#collect parameters
baseline_parameters <- list(
    type = type,
    t_projection_starts = t_projection_starts,
    t_projection_ends = t_projection_ends,
    age_group_sizes = age_group_sizes,
    death_rates = death_rates,
    tt_death_rates = tt_death_rates,
    birth_rates = birth_rates,
    tt_birth_rates = tt_birth_rates,
    duration_of_infectious = NULL,
    duration_of_pre_infectious = NULL,
    duration_of_maternal_immunity = duration_of_maternal_immunity,
    additional_parameters = additional_parameters,
    vaccine_names = vaccine_names,
    force_of_infection = foi,
    tt_force_of_infection = tt_foi,
    vaccine_efficacy = vaccine_efficacy,
    vaccine_efficacy_disease = vaccine_efficacy_disease,
    vaccinations = vaccinations,
    tt_vaccinations = tt_vaccinations,
    duration_of_immunity = duration_of_immunity,
    S_0 = S_0,
    R_0 = R_0,
    M_0 = M_0,
    I_0 = NULL
)

source("analysis/scenario_funcs.R")

vaccine_formulations <- unique(disease_map)
names(vaccine_formulations) <- vaccine_formulations

date_vaccine_coverage <- map(vaccine_formulations, ~c(date_crisis_start, date_crisis_start + 90))

vaccine_coverage_escalation <- map(vaccine_formulations, ~c(0.05, 0.05))
vaccine_coverage_status_quo <- map(vaccine_formulations, ~c(0.1, 0.2))
vaccine_coverage_ceasefire <- map(vaccine_formulations, ~c(0.35, 0.55))

#convert to per disease
date_vaccine_coverage <- map(disease_map, ~date_vaccine_coverage[[.x]])
vaccine_coverage_escalation <- map(disease_map, ~vaccine_coverage_escalation[[.x]])
vaccine_coverage_status_quo <- map(disease_map, ~vaccine_coverage_status_quo[[.x]])
vaccine_coverage_ceasefire <- map(disease_map, ~vaccine_coverage_ceasefire[[.x]])

maintain_foi <- TRUE

escalation_parameters <- apply_scenario(
    baseline_parameters = baseline_parameters,
    vaccine_coverage = vaccine_coverage_escalation,
    date_vaccine_coverage = date_vaccine_coverage,
    date_start = date_start,
    date_crisis_start = date_crisis_start,
    maintain_foi = maintain_foi
)

status_quo_parameters <- apply_scenario(
    baseline_parameters = baseline_parameters,
    vaccine_coverage = vaccine_coverage_status_quo,
    date_vaccine_coverage = date_vaccine_coverage,
    date_start = date_start,
    date_crisis_start = date_crisis_start,
    maintain_foi = maintain_foi
)

ceasefire_parameters <- apply_scenario(
    baseline_parameters = baseline_parameters,
    vaccine_coverage = vaccine_coverage_ceasefire,
    date_vaccine_coverage = date_vaccine_coverage,
    date_start = date_start,
    date_crisis_start = date_crisis_start,
    maintain_foi = maintain_foi
)

res <- list(
    escalation = do.call(project_point_estimate, escalation_parameters),
    status_quo = do.call(project_point_estimate, status_quo_parameters),
    ceasefire = do.call(project_point_estimate, ceasefire_parameters)
) %>%
    transpose() %>%
    map(~map_dfr(.x, function(x) x, .id = "scenario"))

#reformat age groups
output_age_group_names <- c(
    "0mo", "1to11mo", "12to59mo", "5to9yo", "10to14yo", "15to19yo", "20to29yo", "30to39yo", "40to49yo", "50to59yo", "60to100yo"
)
output_age_group_end_adults <- c(1/12, 1, 5, 10, 15, 20, 30, 40, 50, 60, 100)

model_age_group_start <- cumsum(c(0, age_group_sizes))/365
model_age_group_end <- cumsum(c(age_group_sizes, 100*365 - sum(age_group_sizes)))/365

index <- map_int(model_age_group_end, ~min(which(.x <= output_age_group_end_adults)))

res$demographics %>%
    mutate(
        age_group = factor(output_age_group_names[index[age_group]], levels = output_age_group_names, ordered = TRUE),
        date = date_start + t
    ) %>%
    group_by(scenario, age_group, date) %>%
    summarise(
        population = sum(value),
        .groups = "drop"
    ) %>%
    arrange(date, age_group, scenario) %>%
    saveRDS("data/derived/demographics.rds")

projections_full <- res$projections %>%
    mutate(
        age_group = factor(output_age_group_names[index[age_group]], levels = output_age_group_names, ordered = TRUE),
        date = date_start + t
    ) %>%
    group_by(scenario, age_group, vaccine_type, date) %>%
    summarise(
        population = sum(Population),
        immune = sum(Immune),
        immune_disease = sum(`Immune(Disease)`),
        .groups = "drop"
    ) %>%
    arrange(date, age_group, vaccine_type, scenario)

saveRDS(projections_full, "data/derived/projections_full.rds")

output_age_group_names <- c(
    "0mo", "1to11mo", "2nd", "3rd", "4th", "5th", "6th", "extra"
)

output_age_group_end <- c(1/12, 1, 2, 3, 4, 5, 6, 100)

model_age_group_start <- cumsum(c(0, age_group_sizes))/365
model_age_group_end <- cumsum(c(age_group_sizes, 100*365 - sum(age_group_sizes)))/365

index <- map_int(model_age_group_end, ~min(which(.x <= output_age_group_end)))

projections_children <- res$projections %>%
    mutate(
        age_group = factor(output_age_group_names[index[age_group]], levels = output_age_group_names, ordered = TRUE),
        date = date_start + t
    ) %>%
    group_by(scenario, age_group, vaccine_type, date) %>%
    summarise(
        population = sum(Population),
        immune = sum(Immune),
        immune_disease = sum(`Immune(Disease)`),
        .groups = "drop"
    ) %>%
    arrange(date, age_group, vaccine_type, scenario) %>%
    filter(age_group != "extra")

saveRDS(projections_children, "data/derived/projections_children.rds")

#projections_2 <- res$projections %>%
#    mutate(
#        date = date_start + t
#    ) %>%
#    group_by(scenario, age_group, vaccine_type, date) %>%
#    summarise(
#        population = sum(Population),
#        immune = sum(Immune),
#        immune_disease = sum(`Immune(Disease)`),
#        .groups = "drop"
#    ) %>%
#    arrange(date, age_group, vaccine_type, scenario) %>%
#    filter(date == min(date))
#
#susceptible_plots <- projections_2 %>%
#    filter(age_group <= 27) %>%
#    split(~vaccine_type) %>%
#    map(summarise_susceptible)
#ggsave(
#   filename = "plots/susceptibles_og.pdf", 
#   plot = gridExtra::marrangeGrob(susceptible_plots, nrow=1, ncol=1), 
#   width = 15, height = 9
#)
#
#projections_2 %>% filter(vaccine_type == "diphtheria" & scenario == "status_quo" & age_group %in% seq(4,6))
#
#
#pars <- status_quo_parameters
#pars$t_projection_starts <- NULL
#pars$t_projection_ends <- NULL
#pars$t <- c(0, seq(status_quo_parameters$t_projection_starts, status_quo_parameters$t_projection_ends, 1))
#formatted_pars <- do.call(IVODE:::collate_parameters, pars)[[1]]
#do.call(simulate, formatted_pars) %>%
#    format_output(c("Immune", "Immune(Vaccine)", "Immune(Maternal)", "Population")) %>%
#    filter(age_group %in% c(4, 5, 6)) %>%
#    pivot_wider(names_from = compartment, values_from = value) %>% 
#    transmute(
#      t = t, age_group = age_group, vacc_protect = `Immune(Vaccine)`/Population, maternal_protect = `Immune(Maternal)`/Population, total_protect = Immune/Population
#    )
#
#