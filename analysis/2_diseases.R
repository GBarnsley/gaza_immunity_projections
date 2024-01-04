#what diseases/vaccines are we modelling (meningitis not routinely given)
vaccine_names <- c("diphtheria", "measles", "pertussis", "polio - wildtype", "polio - wildtype - vaccine-derived", "Hib disease", "pneumococcal disease", "rotavirus")

#Force of Infection data
#really we need to adjust this FoI for vaccinated/immune i.e. susceptibles actually have a much higher FoI than the average person

#notes
#diphtheria - no cases
#measles - few cases in 2020/2021
#pertussis - few cases in 2022 but no data before seems to be tranmission in lebanon
#polio - wildtype - no cases
#Tuberculosis - cases from 2010 - 2022

foi <- list(
    diphtheria = 0, polio - wildtype = 0
)
tt_foi <- list(
    diphtheria = 0, polio - wildtype = 0
)

IR <- read_csv("data/raw/disease_IRs.csv", show_col_types = FALSE)

locations <- str_split_i(as.character(IR[1,-c(1:3)]), "\n", 1)
years <- as.numeric(names(IR)[-c(1:3)])

IR <- tibble(
    location = locations,
    year = years,
    measles = as.numeric(IR[which(IR$`Scientific name` == "measles Virus Infection"), -c(1:3)]),
    pertussis = as.numeric(IR[which(IR$`Scientific name` == "Bordetella pertussis Infection"), -c(1:3)]),
    meningitis = as.numeric(IR[which(IR$`Scientific name` == "Neisseria meningitidis Infection"), -c(1:3)]),
    tuberculosis_p = as.numeric(IR[which(IR$`Scientific name` == "Mycobacterium tuberculosis Infection (Pulmonary)"), -c(1:3)]),
    tuberculosis_ep = as.numeric(IR[which(IR$`Scientific name` == "Mycobacterium tuberculosis Infection (Extra-pulmonary)"), -c(1:3)])
) %>% 
    fill(year, .direction = "down")
rm(locations, years)
IR <- IR %>%
    pivot_longer(
        cols = -c(location, year),
        names_to = "disease",
        values_to = "IR"
    ) %>%
    mutate(
        location = factor(location, levels = c("Gaza", "Palestine", "West Bank"), ordered = TRUE) #order so that when we sort, Gaza is in priority order
    ) %>%
    group_by(year, disease) %>%
    arrange(
        year, disease, location
    ) %>%
    summarise(
        IR = IR[!is.na(IR)][1],
        .groups = "drop"
    ) %>%
    mutate(
        IR = IR/1000/365 #make per day per person
    ) %>%
    pivot_wider(
        names_from = disease,
        values_from = IR
    ) %>%
    mutate(t = as.numeric(ymd(paste0(year, "-01-01")) - date_start)) %>% 
    arrange(t)
#set missing pertussis to 0 (maybe set to lebanon numbers)
IR$pertussis[is.na(IR$pertussis)] <- 0

foi$measles <- IR$measles
foi$pertussis <- IR$pertussis
foi$Meningitis <- IR$meningitis
foi$Tuberculosis <- IR$tuberculosis_p + IR$tuberculosis_ep
tt_foi$measles <- tt_foi$pertussis <- tt_foi$Meningitis <- tt_foi$Tuberculosis <- IR$t
rm(IR)

map2_dfr(foi, tt_foi, .id = "Pathogen", function(par, tt) {
    tibble(
        Date = date_start + tt,
        `Force of Infection` = par
    )
}) %>%
    group_by(Pathogen) %>%
    complete(Date = seq(date_start, date_projection_end, by = "day")) %>%
    fill(`Force of Infection`, .direction = "down") %>%
    mutate(
        projection_period = Date >= date_projection_start
    ) %>%
    saveRDS("data/derived/foi.rds")

#Vaccine efficacy data (for infection)
vaccine_efficacy <- list(
    diphtheria = list(central = 0.6, lower = 0.51, upper = 0.68),
    measles = list(central = 0.85, lower = 0.25, upper = 0.97),
    pertussis = list(central = 0.83, lower = 0.6, upper = 0.92),
    Meningitis = list(central = NA, lower = NA, upper = NA), #No data?
    polio - wildtype = list(central = 0.28, lower = 0.22, upper = 0.29),
    Tuberculosis = list(central = NA, lower = NA, upper = NA) #No data?
)
vaccine_efficacy <- map(vaccine_efficacy, ~.x$central)

#Vaccine coverage
schedule <- list(
    pentavalent = c(2/12, 4/12, 6/12, 6) * 365,
    BCG = c(1/12) * 365,
    MMR = c(1, 1.5) * 365,
    `IPV-OPV` = c(1/12, 2/12, 4/12, 6/12, 18/12, 6) * 365
)
#convert to age group indicators
age_groups <- cumsum(c(age_group_sizes, 365*10))
schedule <- map(schedule, ~unique(map_int(.x, function(x) min(which(x < age_groups)))))

initial_vaccine_coverage <- list(
    pentavalent = c(0.958, 0.994),
    BCG = c(0.995),
    MMR = c(0.986), #somewhere between meales and MMR
    `IPV-OPV` = c(0.983, 0.997, 0.995) #based on OPV3
)

#where we have just a single value we assume thats the coverage for both doses
vaccine_coverage <- tibble(
    year = c(2010, 2014, 2019),
    pentavalent = c(0.925, 0.999, 0.973),
    BCG = c(0.997, 0.998, 1),
    MMR = c(0.959, 0.994, 0.965),
    `IPV-OPV` = c(0.951, 0.998, 0.942)
) %>%
    rbind(
        tibble(
            year = c(2022, 2021, 2020, 2019, 2018, 2017, 2016, 2015),
            pentavalent = c(98.1, 96.8, 98.4, 99.3, 98.2, 98.1, 98.8, 98.3)/100,
            BCG = c(98.7, 99.9, 99.0, 99.5, 98.8, 99.4, 99.8, 99)/100,
            MMR = c(94.3, 98.4, 98.3, 99.9, 98.9, 99.6, 99.0, 99.8)/100,
            `IPV-OPV` = c(98.8, 98.3, 98.6, 99.8, 98.9, 98.3, 98.4, 98.3)/100
        )
    ) %>%
    group_by(year) %>%
    summarise_all(mean) %>%
    arrange(year) %>%
    mutate(
        year = year - 1 #assume each one represents the coverage at the end of the year
    ) %>%
    complete(year = seq(year(date_start), max(year))) %>%
    fill(everything(), .direction = "up") %>%
    group_by(pentavalent, BCG, MMR, `IPV-OPV`) %>%
    summarise(
        year = min(year),
        .groups = "drop"
    ) %>%
    arrange(year)
tt_vaccine_coverage <- list()
tt_vaccine_coverage$measles <- tt_vaccine_coverage$pertussis <-
    tt_vaccine_coverage$Meningitis <- tt_vaccine_coverage$Tuberculosis <-
        as.numeric(ymd(paste0(vaccine_coverage$year, "-01-01")) - date_start)
vaccine_coverage <- vaccine_coverage  %>%
    select(!year) %>% 
    as.list() %>%
    imap(function(coverage, name) {
        vaccine_schedule <- schedule[[name]]
        coverage_matrix <- matrix(0, nrow = length(coverage), ncol = length(age_groups))
        coverage_matrix[,vaccine_schedule] <- coverage
        coverage_matrix
    })
#convert to per disease
disease_map <- list(
    Dipetheria = "pentavalent",
    pertussis = "pentavalent",
    Tuberculosis = "BCG",
    measles = "MMR",
    polio - wildtype = "IPV-OPV"
)
vaccine_coverage <- map(disease_map, ~vaccine_coverage[[.x]])
rm(disease_map, schedule, age_groups)

map2_dfr(vaccine_coverage, tt_foi, .id = "Pathogen", function(par, tt) {
    tibble(
        Date = date_start + tt,
        `Force of Infection` = par
    )
}) %>%
    group_by(Pathogen) %>%
    complete(Date = seq(date_start, date_projection_end, by = "day")) %>%
    fill(`Force of Infection`, .direction = "down") %>%
    mutate(
        projection_period = Date >= date_projection_start
    ) %>%
    saveRDS("data/derived/foi.rds")

#Duration of Immunity
duration_of_immunity <- list(
    diphtheria = list(central = 0.4, lower = 1.111, upper = 0.25),
    measles = list(central = 0.009, lower = 0.016, upper = 0.005),
    pertussis = list(central = 0.043, lower = 0.059, upper = 0.034),
    Meningitis = list(central = NA, lower = NA, upper = NA), #No data?
    polio - wildtype = list(central = 0.05, lower = 0.067, upper = 0.033),
    Tuberculosis = list(central = NA, lower = NA, upper = NA) #No data?
)
duration_of_immunity <- map(duration_of_immunity, ~1/.x$central*365)
