#what diseases/vaccines are we modelling (meningitis not routinely given)
vaccine_names <- c(
    "diphtheria", "measles", "pertussis", "polio - wildtype",
    "polio - vaccine-derived", "Hib disease", "pneumococcal disease", "rotavirus"
)

disease_map <- list(
    diphtheria = "pentavalent", measles = "MMR", pertussis = "pentavalent",
    `polio - wildtype` = "IPV-OPV", `polio - vaccine-derived` = "IPV-OPV",
    `Hib disease` = "pentavalent", `pneumococcal disease` = "PCV",
    rotavirus = "ROTAVAC"
)

names(vaccine_names) <- vaccine_names

disease_parameters <- map(vaccine_names, IVODE::sample_parameters)

#Force of Infection data
#really we need to adjust this FoI for vaccinated/immune i.e. susceptibles actually have a much higher FoI than the average person

#notes
#diphtheria - no cases - medium reproduction number, likely large background immunity - minor importations
#measles - few cases in 2020/2021 - high reproduction number, likely large background immunity - minor importations
#pertussis - few cases in 2022 but no data before seems to be tranmission in lebanon - high reproduction number, likely large background immunity - endemic (some evidence that vaccine has no impact on tranmission)
#polio - wildtype - no cases - high reproduction number, likely large background immunity - eliminated no tranmission since mid-90s
#polio - vaccine-derived - no cases - high reproduction number, likely large background immunity? - tranmission in Israel
#Hib disease - cases of meningitis suggests transmission - likely low as protects against carriage (need to adjust for this, potentially first carriage -> disease)
#pneumococcal disease - will be carriage in the population, how does this translate to immunity
#rotavirus - nearly everyone gets infected at an early age high reproduction number, likely large background immunity - endemic

high_background_immunity <- c("diphtheria", "measles", "pertussis", "polio - wildtype", "rotavirus")
names(high_background_immunity) <- high_background_immunity
moderate_background_immunity <- c("polio - vaccine-derived", "Hib disease", "pneumococcal disease")
names(moderate_background_immunity) <- moderate_background_immunity

use_data <- c("measles", "Hib disease")
names(use_data) <- use_data

endemic <- c("pertussis", "rotavirus", "pneumococcal disease")
names(endemic) <- endemic

importations <- c("diphtheria")
names(importations) <- importations

eliminated <- c("polio - wildtype", "polio - vaccine-derived")
names(eliminated) <- eliminated

foi <- list()
foi <- c(foi, map(importations, ~disease_parameters[[.x]]$R0)) #(need to think about this too)
foi <- c(foi, map(endemic, ~disease_parameters[[.x]]$R0))
foi <- c(foi, map(use_data, ~disease_parameters[[.x]]$R0)) #need to think about how to include the data

tt_foi <- map(foi, ~0)

additional_parameters <- map(foi, ~list(
    prop_death = additional_parameters$prop_death
))
rm(use_data, endemic, importations)

map2_dfr(foi, tt_foi, .id = "Pathogen", function(par, tt) {
    tibble(
        Date = date_start + tt,
        `Basic Reproduction Number` = par
    )
}) %>%
    group_by(Pathogen) %>%
    complete(Date = seq(date_start, date_projection_end, by = "day")) %>%
    fill(`Basic Reproduction Number`, .direction = "down") %>%
    mutate(
        projection_period = Date >= date_crisis_start
    ) %>%
    saveRDS("data/derived/foi.rds")

#Vaccine efficacy data (for infection)
vaccine_efficacy <- read_csv("data/raw/vaccine_parameters.csv", show_col_types = FALSE) %>%
    filter(
        `Outcome of vaccination (direct/individual protection only)` == "Vaccine effectiveness against infection (% of fully vaccinated people who are protected against acquiring infection)"
    ) %>%
    rowwise() %>%
    transmute(
        disease = c(
            "Diphtheria (as part of pentavalent vaccine)" = "diphtheria",
            "Pertussis (as part of pentavalent vaccine)" = "pertussis",
            "Measles (as part of MMR)" = "measles",
            "Polio wild-type 1 or 3 / IPV-OPV sequential" = "polio - wildtype",
            "Polio vaccine-derived cVPD2 / IPV-OPV sequential" = "polio - vaccine-derived",
            "Haemophilus influenzae type B (as part of pentavalent vaccine)" = "Hib disease", 
            "Pneumococcus (conjugate vaccine)" = "pneumococcal disease",
            "Rotavirus (1-valent)" = "rotavirus"
        )[`Infectious disease / vaccine`],
        vaccine_efficacy = list(list(central = `Central estimate`, lower = `Lower bound`, upper = `Upper bound`))
    ) %>%
    filter(disease %in% vaccine_names)
vaccine_efficacy <- setNames(vaccine_efficacy$vaccine_efficacy, vaccine_efficacy$disease)

vaccine_efficacy <- map(vaccine_efficacy, ~.x$central)


#Vaccine efficacy data (for disease)
vaccine_efficacy_disease <- read_csv("data/raw/vaccine_parameters.csv", show_col_types = FALSE) %>%
    filter(
        `Outcome of vaccination (direct/individual protection only)` %in% c(
            "Vaccine effectiveness against severe disease (% of fully vaccinated people who are protected against acquiring severe disease)",
            #"Vaccine effectiveness against pneumonia (% of fully vaccinated people who are protected against acquiring severe disease)",
            "Vaccine effectiveness against meningitis or other invasive disease (% of fully vaccinated people who are protected against acquiring severe disease)"
        )
    ) %>%
    rowwise() %>%
    transmute(
        disease = c(
            "Diphtheria (as part of pentavalent vaccine)" = "diphtheria",
            "Pertussis (as part of pentavalent vaccine)" = "pertussis",
            "Measles (as part of MMR)" = "measles",
            "Polio wild-type 1 or 3 / IPV-OPV sequential" = "polio - wildtype",
            "Polio vaccine-derived cVPD2 / IPV-OPV sequential" = "polio - vaccine-derived",
            "Haemophilus influenzae type B (as part of pentavalent vaccine)" = "Hib disease", 
            "Pneumococcus (conjugate vaccine)" = "pneumococcal disease",
            "Rotavirus (1-valent)" = "rotavirus"
        )[`Infectious disease / vaccine`],
        vaccine_efficacy = list(list(central = `Central estimate`, lower = `Lower bound`, upper = `Upper bound`))
    ) %>%
    filter(disease %in% vaccine_names)
vaccine_efficacy_disease <- setNames(vaccine_efficacy_disease$vaccine_efficacy, vaccine_efficacy_disease$disease)

vaccine_efficacy_disease <- map(vaccine_efficacy_disease, ~.x$central)

#fix polio-wildtype
vaccine_efficacy_disease$`polio - wildtype` <- max(vaccine_efficacy$`polio - wildtype`, vaccine_efficacy_disease$`polio - wildtype`)

#Vaccine coverage
schedule <- list(
    pentavalent = c(2/12, 4/12, 6/12, 6),
    BCG = c(1/12),
    MMR = c(1, 1.5),
    `IPV-OPV` = c(1/12, 2/12, 4/12, 6/12, 18/12, 6),
    PCV = c(2/12, 4/12, 1),
    ROTAVAC = c(2/12, 4/12, 6/12)
)
#simplify this so that first 6 months doses happen at 4 months, 1 year at 1 year, 1.5 years at 1.5 years
schedule <- map(schedule, ~unique(case_when(
    .x %in% c(2/12, 4/12, 6/12) ~ 4/12,
    .x == 1/12 ~ 0,
    TRUE ~ .x
)) * 365)
#convert to age group indicators
age_groups <- c(0, cumsum(age_group_sizes))
schedule <- map(schedule, ~unique(map_int(.x, function(x) min(which(x < age_groups)) - 1))) #shift back by one since vaccinations happens upon entering an age group

initial_vaccine_coverage <- list( #from 2000
    pentavalent = c(0.958, 0.994),
    MMR = c(0.986), #somewhere between measles and MMR
    `IPV-OPV` = c(0.983, 0.997, 0.995), #based on OPV3
    PCV = c(0, 0),
    ROTAVAC = c(0)
)

#where we have just a single value we assume thats the coverage for both doses
vaccine_coverage <- tibble(
    year = c(2010, 2014, 2019),
    pentavalent = c(0.925, 0.999, 0.973),
    BCG = c(0.997, 0.998, 1),
    MMR = c(0.959, 0.994, 0.965),
    `IPV-OPV` = c(0.951, 0.998, 0.942),
    PCV = c(NA, NA, 0.951),
    ROTAVAC = c(NA, NA, 0.916)
) %>%
    rbind(
        tibble(
            year = c(2022, 2021, 2020, 2019, 2018, 2017, 2016, 2015),
            pentavalent = c(98.1, 96.8, 98.4, 99.3, 98.2, 98.1, 98.8, 98.3)/100,
            BCG = c(98.7, 99.9, 99.0, 99.5, 98.8, 99.4, 99.8, 99)/100,
            MMR = c(94.3, 98.4, 98.3, 99.9, 98.9, 99.6, 99.0, 99.8)/100,
            `IPV-OPV` = c(98.8, 98.3, 98.6, 99.8, 98.9, 98.3, 98.4, 98.3)/100,
            PCV = NA,
            ROTAVAC = NA
        )
    ) %>%
    group_by(year) %>%
    summarise_all(~mean(.x, na.rm = TRUE)) %>%
    arrange(year) %>%
    mutate(
        year = year - 1 #assume each one represents the coverage at the end of the year
    ) %>%
    complete(year = seq(year(date_start), max(year))) %>%
    fill(everything(), .direction = "updown") %>%
    group_by(across(c(-year))) %>%
    summarise(
        year = min(year),
        .groups = "drop"
    ) %>%
    arrange(year)

#rota and PCV not introduced at the start
vaccine_coverage$PCV[vaccine_coverage$year < 2012] <- 0
vaccine_coverage$ROTAVAC[vaccine_coverage$year < 2016] <- 0

vaccine_coverage %>%
    mutate(Date = ymd(paste0(year, "-01-01"))) %>%
    complete(Date = seq(date_start, date_projection_end, by = "day")) %>%
    select(-year) %>% 
    fill(-Date, .direction = "down") %>%
    pivot_longer(
        cols = -Date,
        names_to = "Vaccine",
        values_to = "Coverage"
    ) %>%
    mutate(
        projection_period = Date >= date_crisis_start
    ) %>%
    saveRDS("data/derived/vaccine_coverage.rds")

tt_vaccinations <- map(vaccine_names, ~as.numeric(ymd(paste0(vaccine_coverage$year, "-01-01")) - date_start))

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
vaccinations <- map(disease_map, ~vaccine_coverage[[.x]])
rm(schedule, age_groups, vaccine_coverage)

#Duration of Immunity
duration_of_immunity <- read_csv("data/raw/vaccine_parameters.csv", show_col_types = FALSE) %>%
    filter(
        `Outcome of vaccination (direct/individual protection only)` == "Immunity waning rate per year (1 / mean duration of functional vaccine protection against infection and/or disease)"
    ) %>%
    rowwise() %>%
    transmute(
        disease = c(
            "Diphtheria (as part of pentavalent vaccine)" = "diphtheria",
            "Pertussis (as part of pentavalent vaccine)" = "pertussis",
            "Measles (as part of MMR)" = "measles",
            "Polio wild-type 1 or 3 / IPV-OPV sequential" = "polio - wildtype",
            "Polio vaccine-derived cVPD2 / IPV-OPV sequential" = "polio - vaccine-derived",
            "Haemophilus influenzae type B (as part of pentavalent vaccine)" = "Hib disease", 
            "Pneumococcus (conjugate vaccine)" = "pneumococcal disease",
            "Rotavirus (1-valent)" = "rotavirus"
        )[`Infectious disease / vaccine`],
        duration_of_immunity = list(list(central = `Central estimate`, lower = `Lower bound`, upper = `Upper bound`))
    ) %>%
    filter(disease %in% vaccine_names)
duration_of_immunity <- setNames(duration_of_immunity$duration_of_immunity, duration_of_immunity$disease)
duration_of_immunity <- map(duration_of_immunity, ~(1/.x$central) * 365)

#duration of pre-infectious period
duration_of_pre_infectious <- map(disease_parameters, ~.x$duration_of_pre_infectious)

#duration of infectious period
duration_of_infectious <- map(disease_parameters, ~.x$duration_of_infectious)

#initial states
#no vaccinated, shouldn't matter that much (explore)

#if a disease has infections at the start then we assume some maternal and acquired immunity
M_0 <- map(vaccine_names, ~NULL) #can't assume any atm
p_infected <- 0.001
p_moderate <- 0.75
p_high <- 0.95
S_0 <- map(moderate_background_immunity, ~(1-p_moderate) * total_pop) %>% 
    c(
        map(high_background_immunity, ~(1-p_high) * total_pop)
    )

I_0 <- imap(S_0, \(x, disease) if(disease %in% eliminated) {x*0} else {x * p_infected})

S_0 <- imap(S_0, \(x, disease) if(disease %in% eliminated) {x} else {x * (1 - p_infected)})

R_0 <- map(moderate_background_immunity, ~p_moderate * total_pop) %>% 
    c(
        map(high_background_immunity, ~p_high * total_pop)
    )
