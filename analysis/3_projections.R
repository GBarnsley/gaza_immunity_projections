#What type of model
type <- "deterministic_gz"

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
    M_0 = M_0
)

source("analysis/scenario_funcs.R")

vaccine_formulations <- unique(disease_map)
names(vaccine_formulations) <- vaccine_formulations

date_vaccine_coverage <- map(vaccine_formulations, ~c(date_projection_start, date_projection_start + 90))

vaccine_coverage_pessimistic <- map(vaccine_formulations, ~c(0.05, 0.05))
vaccine_coverage_central <- map(vaccine_formulations, ~c(0.1, 0.2))
vaccine_coverage_pessimistic <- map(vaccine_formulations, ~c(0.35, 0.55))

#convert to per disease
date_vaccine_coverage <- map(disease_map, ~date_vaccine_coverage[[.x]])
vaccine_coverage_pessimistic <- map(disease_map, ~vaccine_coverage_central[[.x]])
vaccine_coverage_central <- map(disease_map, ~vaccine_coverage_central[[.x]])
vaccine_coverage_optimistic <- map(disease_map, ~vaccine_coverage_optimistic[[.x]])

pessimistic_parameters <- apply_scenario(
    baseline_parameters = baseline_parameters,
    vaccine_coverage = vaccine_coverage_pessimistic,
    date_vaccine_coverage = date_vaccine_coverage,
    date_start = date_start,
    date_projection_start = date_projection_start
)

central_parameters <- apply_scenario(
    baseline_parameters = baseline_parameters,
    vaccine_coverage = vaccine_coverage_central,
    date_vaccine_coverage = date_vaccine_coverage,
    date_start = date_start,
    date_projection_start = date_projection_start
)

optimistic_parameters <- apply_scenario(
    baseline_parameters = baseline_parameters,
    vaccine_coverage = vaccine_coverage_optimistic,
    date_vaccine_coverage = date_vaccine_coverage,
    date_start = date_start,
    date_projection_start = date_projection_start
)

res <- list(
    pessimistic = do.call(project_point_estimate, pessimistic_parameters),
    central = do.call(project_point_estimate, central_parameters),
    optimistic = do.call(project_point_estimate, optimistic_parameters)
) %>%
    transpose() %>%
    map(~map_dfr(.x, function(x) x, .id = "scenario"))

#reformat age groups
output_age_group_names <- c(
    "0mo", "1to11mo", "12to59mo", "5to9yo", "10to14yo", "15to19yo", "20to29yo", "30to39yo", "40to49yo", "50to59yo", "60to100yo"
)
output_age_group_end <- c(1/12, 1, 5, 10, 15, 20, 30, 40, 50, 60, 100)

model_age_group_start <- cumsum(c(0, age_group_sizes))/365
model_age_group_end <- cumsum(c(age_group_sizes, 100*365 - sum(age_group_sizes)))/365

index <- map_int(model_age_group_end, ~min(which(.x <= output_age_group_end)))

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

projections <- res$projections %>%
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

saveRDS(projections, "data/derived/projections.rds")
