#What type of model
type <- "deterministic_gz"

t_projection_starts <- as.numeric(date_projection_start - date_start)
t_projection_ends <- as.numeric(date_projection_end - date_start)

#initial state of the model

res <- project_point_estimate(
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
    vaccinations = vaccinations,
    tt_vaccinations = tt_vaccinations,
    duration_of_immunity = duration_of_immunity,
    S_0 = S_0,
    R_0 = R_0,
    M_0 = M_0
)

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
    group_by(age_group, date) %>%
    summarise(
        population = sum(value),
        .groups = "drop"
    ) %>%
    arrange(date, age_group) %>%
    saveRDS("data/derived/demographics.rds")

projections <- res$projections %>%
    mutate(
        age_group = factor(output_age_group_names[index[age_group]], levels = output_age_group_names, ordered = TRUE),
        date = date_start + t
    ) %>%
    group_by(age_group, vaccine_type, date) %>%
    summarise(
        population = sum(Population),
        immune = sum(Immune),
        .groups = "drop"
    ) %>%
    arrange(date, age_group, vaccine_type)

saveRDS(projections, "data/derived/projections.rds")