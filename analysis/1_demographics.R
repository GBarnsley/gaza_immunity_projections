#define age groups
#age_group_sizes <- c(
#    1/12, #0-1 month
#    11/12, #1-12 months
#    4, #1-5 years
#    10, #5-15 years
#    15, #15-30 years
#    15, #30-45 years
#    5, #45-50 years
#    5, #50-55 years
#    5, #55-60 years
#    5, #60-65 years
#    5, #65-70 years
#    5 #70-75 years
#    #75 + don't define a size
#) * 365 #in days

age_group_sizes <- c(
    1/12, #0-1 month
    11/12, #1-12 months #we need vaccination years to be yearly
    1, #1-2 years
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
##we will weight them for our age groups using the 2006 fertility rates for gaza see palestine_fertility_rates.pdf
#fertility_rates <- tibble(
#    age_group_size = c(15, rep(5, 7))*365,
#    fertility_rate = c(0, 60.6, 242.9, 273.6, 241.1, 161.8, 64.2, 7.2)
#) %>%
#    mutate(
#        age_group_start = cumsum(c(0, age_group_size[-length(age_group_size)])),
#        age_group_end = cumsum(age_group_size)
#    ) %>%
#    filter(
#        age_group_start > 0
#    )
#
##fit a gamma distribution (via mle)
#ggplot(fertility_rates, aes(x = age_group_start, xend = age_group_end, y = fertility_rate, yend = fertility_rate)) +
#    geom_segment()

#for now we just use the rates per total pop and split based on proportion immune
crude_birth_rate <- read_csv("data/raw/birth_rate.csv", show_col_types = FALSE) #CIA
crude_birth_rate <- tibble(
    year = as.numeric(crude_birth_rate[1,-1]),
    crude_birth_rate = as.numeric(crude_birth_rate[2,-1])
)
tt_birth_rate <- as.numeric(ymd(paste0(crude_birth_rate$year, "-01-01")) - date_start)
birth_rate <- crude_birth_rate$crude_birth_rate/1000/365 #convert to days per person

tibble(
    Date = date_start + tt_birth_rate,
    `Crude Birth Rate` = birth_rate
) %>%
    complete(Date = seq(date_start, date_projection_end, by = "day")) %>%
    fill( `Crude Birth Rate`, .direction = "down") %>%
    mutate(
        projection_period = Date >= date_projection_start
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
    proportional_death_rate = additional_parameters$prob_death
) %>%
    saveRDS("data/derived/age_proportional_death_rate.rds")

rm(death_rate_age, model_age_group_start, model_age_group_end, missing_starts)

#for now we just use this as scaling assuming underlying population distribution remains unchanged, come back to this
#also need an underlying age distribution, wait until we have this
death_rate <- read_csv("data/raw/crude_death_rate.csv", show_col_types = FALSE) %>% #WB
    arrange(Year)
tt_death_rate <- as.numeric(ymd(paste0(death_rate$Year, "-01-01")) - date_start)
death_rate <- death_rate$Value/1000/365 #convert to days per person
death_rate <- death_rate[tt_death_rate >= 0]
tt_death_rate <- tt_death_rate[tt_death_rate >= 0]

tibble(
    Date = date_start + tt_death_rate,
    `Crude Death Rate` = death_rate
) %>%
    complete(Date = seq(date_start, date_projection_end, by = "day")) %>%
    fill(`Crude Death Rate`, .direction = "down") %>%
    mutate(
        projection_period = Date >= date_projection_start
    ) %>%
    saveRDS("data/derived/crude_death_rate.rds")

#duration of maternal immunity
duration_maternal_immunity <- (6/12) * 365 - age_group_sizes[1] #6 months, adjusted for only occuring in the 2nd compartment

