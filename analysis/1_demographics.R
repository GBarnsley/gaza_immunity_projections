
age_group_sizes <- c(
    1/12, #0-1 month
    3/12, #1-4 months
    8/12, #4-12 months
    6/12, #12-18 months
    6/12, #18-24 months
    1, #2-3 years
    1, #3-4 years
    1, #4-5 years
    1, #5-6 years
    4, #6-10 years
    5, #10-15 years
    5, #15-20 years
    10, #20-30 years
    10, #30-40 years
    10, #40-50 years
    10 #50-60 years
    #60+ don't define a size
) * 365 #in days

age_group_sizes <- c(
    rep(1/12, 12), #0-1 year in months
    rep(1/12, 12), #1-2 years in months
    1, #2-3 years
    1, #3-4 years
    1, #4-5 years
    1, #5-6 years
    4, #6-10 years
    5, #10-15 years
    5, #15-20 years
    10, #20-30 years
    10, #30-40 years
    10, #40-50 years
    10 #50-60 years
    #60+ don't define a size
) * 365 #in days

#now demographic data
#crude birth rates are per total population
#for now we just use the rates per total pop and split based on proportion immune
crude_birth_rate <- read_csv("data/raw/birth_rate.csv", show_col_types = FALSE) #CIA
crude_birth_rate <- tibble(
    year = as.numeric(crude_birth_rate[1,-1]),
    crude_birth_rate = as.numeric(crude_birth_rate[2,-1])
)
tt_birth_rates <- as.numeric(ymd(paste0(crude_birth_rate$year, "-01-01")) - date_start)
birth_rates <- crude_birth_rate$crude_birth_rate/1000/365 #convert to days per person

tibble(
    Date = date_start + tt_birth_rates,
    `Crude Birth Rate` = birth_rates
) %>%
    complete(Date = seq(date_start, date_projection_end, by = "day")) %>%
    fill( `Crude Birth Rate`, .direction = "down") %>%
    mutate(
        projection_period = Date >= date_crisis_start
    ) %>%
    saveRDS("data/derived/crude_birth_rate.rds")
rm(crude_birth_rate)
#Death rates
#death rates by age
death_rate_age <- read_csv("data/raw/crude_death_by_age.csv", show_col_types = FALSE) %>% #MoH
    mutate(
        age_group_start = if_else(
            `Age group` == "<1Y",
            0,
            as.numeric(str_extract(`Age group`, "\\d+"))
        )*365,
        age_group_end = lead(age_group_start, 1, default = NA)
    ) %>%
    select(!c(`Age group`, `2020`, `2021`)) %>% #covid years
    pivot_longer(
        cols = -c(age_group_start, age_group_end),
        names_to = "year",
        values_to = "death_rate"
    ) %>%
    group_by(age_group_start, age_group_end) %>%
    summarise(
        death_rate = mean(death_rate, na.rm = TRUE),
        .groups = "drop"
    )
#really we should fit a mixture distribution

#instead assume that within these age groups the distribution is uniform and the age distibution is uniform
#simplifies things
model_age_group_start <- cumsum(c(0, age_group_sizes))
model_age_group_end <- c(cumsum(age_group_sizes), NA)

#add missing start times
missing_starts <- setdiff(model_age_group_start, death_rate_age$age_group_start)
death_rate_age <- death_rate_age %>%
    rbind(
        tibble(
            age_group_start = missing_starts,
            age_group_end = NA,
            death_rate = NA
        )
    ) %>%
    arrange(age_group_start) %>%
    mutate(
        age_group_end = lead(age_group_start, 1)
    ) %>% 
    fill(death_rate, .direction = "down")
#now recombine into only the model age groups
death_rate_age$group <- map2_int(
    death_rate_age$age_group_start, death_rate_age$age_group_end, function(start, end) {
        if(is.na(end)) return(NA)
        if(sum(model_age_group_start >= end) == 0) return(NA)
        which(
            model_age_group_start <= start & model_age_group_end >= end
        )
    })
death_rate_age$age_group_end[is.na(death_rate_age$age_group_end)] <- death_rate_age$age_group_start[is.na(death_rate_age$age_group_end)] + (10*365)
death_rate_age$group[is.na(death_rate_age$group)] <- max(death_rate_age$group, na.rm = TRUE) + 1

#now recombine weighting based on width of age groups
death_rate_age <- death_rate_age %>%
    group_by(group) %>%
    summarise(
        death_rate = sum(death_rate * (age_group_end - age_group_start), na.rm = TRUE) / sum(age_group_end - age_group_start, na.rm = TRUE),
        age_group_start = min(age_group_start, na.rm = TRUE),
        age_group_end = max(age_group_end, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    select(-group)

additional_parameters <- list(prop_death = death_rate_age$death_rate/sum(death_rate_age$death_rate))
#save for plotting
tibble(
    age_group_start = cumsum(c(0, age_group_sizes))/365,
    age_group_end = c(cumsum(age_group_sizes)/365, 100),
    proportional_death_rate = additional_parameters$prop_death
) %>%
    saveRDS("data/derived/age_proportional_death_rate.rds")

rm(death_rate_age, model_age_group_start, model_age_group_end, missing_starts)

death_rates <- read_csv("data/raw/crude_death_rate.csv", show_col_types = FALSE) %>% #WB
    arrange(Year)
tt_death_rates <- as.numeric(ymd(paste0(death_rates$Year, "-01-01")) - date_start)
death_rates <- death_rates$Value/1000/365 #convert to days per person
death_rates <- death_rates[tt_death_rates >= 0]
tt_death_rates <- tt_death_rates[tt_death_rates >= 0]

tibble(
    Date = date_start + tt_death_rates,
    `Crude Death Rate` = death_rates
) %>%
    complete(Date = seq(date_start, date_projection_end, by = "day")) %>%
    fill(`Crude Death Rate`, .direction = "down") %>%
    mutate(
        projection_period = Date >= date_crisis_start
    ) %>%
    saveRDS("data/derived/crude_death_rates.rds")

#duration of maternal immunity
duration_of_maternal_immunity <- (6/12) * 365 #6 months, will be convert into a number of compartments in the model

#Calculate total population in each age group (currently take Palestine breakdown from WPP and split by proportion of population in each area)
age_vars <- c("0-4", "5-9", "10-14",
"15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49",
"50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84",
"85-89", "90-94", "95-99", "100+")
total_population <- read_csv("data/raw/wpp.csv", skip = 16) %>%
    filter(Year == 2000 & `Region, subregion, country or area *` == "State of Palestine" ) %>%
    mutate(
        across(all_of(age_vars), ~as.numeric(.x) * 1000)
    ) %>%
    select(all_of(age_vars)) %>%
    pivot_longer(
        cols = all_of(age_vars),
        names_to = "Age",
        values_to = "Population"
    ) %>%
    filter(Population > 0) %>%
    mutate(
        age_group_start = as.numeric(str_split_i(Age, "-", 1)) * 365,
        age_group_end = (as.numeric(str_split_i(Age, "-", 2)) + 1) * 365
    )
#should really fit an expoential to reformat to the new age groups
model_age_group_start <- c(0, cumsum(age_group_sizes))
model_age_group_end <- cumsum(c(age_group_sizes,  365*10))
missing_starts <- setdiff(model_age_group_end, total_population$age_group_start)

#assume that the population is evenly distributed across the age group 
model_total_population <- total_population %>%
    mutate(old_age_group_size = age_group_end - age_group_start) %>%
    complete(age_group_start = c(model_age_group_start, total_population$age_group_start)) %>%
    mutate(
        age_group_end = lead(age_group_start, 1, default = max(age_group_start) + (5*365)),
        age_group_size = age_group_end - age_group_start,
        Population_to_split = Population
    ) %>%
    arrange(age_group_start) %>%
    fill(Population_to_split, old_age_group_size, .direction = "down") %>%
    mutate(
        Population = Population_to_split * age_group_size / old_age_group_size
    ) %>%
    select(age_group_start, age_group_end, Population)

model_total_population$group <- map2_int(
    model_total_population$age_group_start, model_total_population$age_group_end, function(start, end) {
        if(is.na(end)) return(NA)
        if(sum(model_age_group_start >= end) == 0) return(NA)
        which(
            model_age_group_start <= start & model_age_group_end >= end
        )
    })
model_total_population$group[is.na(model_total_population$group)] <- max(model_total_population$group, na.rm = TRUE) + 1

model_total_population <- model_total_population %>%
    group_by(group) %>%
    summarise(
        Population = sum(Population),
        age_group_start = min(age_group_start, na.rm = TRUE),
        age_group_end = max(age_group_end, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    select(-group)
total_pop <- pull(model_total_population, Population)
rm(model_total_population, total_population, age_vars, missing_starts, model_age_group_start, model_age_group_end)

#adjust for 2000 Westbank Gaza population PCBS
gaza_2000_pop <- 1109677

total_pop <- gaza_2000_pop * total_pop/sum(total_pop)
rm(gaza_2000_pop)
