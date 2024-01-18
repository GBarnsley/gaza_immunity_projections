apply_scenario <- function(baseline_parameters, vaccine_coverage, date_vaccine_coverage, date_start, date_projection_start, maintain_foi = TRUE) {
    parameters <- baseline_parameters
    diseases <- names(parameters$force_of_infection)
    n_age <- length(parameters$age_group_sizes) + 1
    t_projection <- as.numeric(date_projection_start - date_start)

    #remove all vaccinations from the projection start date
    for (disease in diseases) {
        parameters$vaccinations[[disease]] <- rbind(
            parameters$vaccinations[[disease]][parameters$tt_vaccinations[[disease]] < t_projection, ],
            rep(0, n_age)
        )
        parameters$tt_vaccinations[[disease]] <- c(
            parameters$tt_vaccinations[[disease]][parameters$tt_vaccinations[[disease]] < t_projection],
            t_projection
        )
    }

    if(!is.null(vaccine_coverage)) {
        #update according to the vaccine coverage
        for (disease in diseases) {
            existing_coverage <- parameters$vaccinations[[disease]]

            t_coverage <- which(colSums(existing_coverage) > 0)

            dates_to_cover <- date_vaccine_coverage[[disease]][date_vaccine_coverage[[disease]] >= date_projection_start]
            t_new_coverage <- as.numeric(dates_to_cover - date_start)

            new_coverage <- matrix(0, nrow = length(t_new_coverage), ncol = n_age)

            for(i in 1:length(date_vaccine_coverage[[disease]])) {
                new_coverage[i, t_coverage] <- vaccine_coverage[[disease]][i]
            }

            parameters$vaccinations[[disease]] <- rbind(
                existing_coverage, new_coverage
            )
            parameters$tt_vaccinations[[disease]] <- c(
                parameters$tt_vaccinations[[disease]],
                t_new_coverage
            )
        }
    }
    if (!maintain_foi) {
        #remove all force of infection from the projection start date
        for (disease in diseases) {
            parameters$force_of_infection[[disease]] <- c(
                parameters$force_of_infection[[disease]][parameters$tt_force_of_infection[[i]] < t_projection],
                0
            )
            parameters$tt_force_of_infection[[disease]] <- c(
                parameters$tt_force_of_infection[[disease]][parameters$tt_force_of_infection[[i]] < t_projection],
                t_projection
            )
        }
    }
    return(parameters)
}