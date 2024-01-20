apply_scenario <- function(baseline_parameters, vaccine_coverage, date_vaccine_coverage, date_start, date_crisis_start, maintain_foi = TRUE) {
    parameters <- baseline_parameters
    diseases <- names(parameters$force_of_infection)
    n_age <- length(parameters$age_group_sizes) + 1
    t_projection <- as.numeric(date_crisis_start - date_start)

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
            t_existing_coverage <- parameters$tt_vaccinations[[disease]]
            
            t_coverage <- which(colSums(existing_coverage) > 0)

            dates_to_cover <- date_vaccine_coverage[[disease]][date_vaccine_coverage[[disease]] >= date_crisis_start]
            t_new_coverage <- as.numeric(dates_to_cover - date_start)

            existing_coverage <- existing_coverage[t_existing_coverage < t_projection, ]
            t_existing_coverage <- t_existing_coverage[t_existing_coverage < t_projection]

            new_coverage <- matrix(0, nrow = length(t_new_coverage), ncol = n_age)

            for(i in seq_along(t_new_coverage)) {
                new_coverage[i, t_coverage] <- vaccine_coverage[[disease]][i]
            }

            parameters$vaccinations[[disease]] <- rbind(
                existing_coverage, new_coverage
            )
            parameters$tt_vaccinations[[disease]] <- c(
                t_existing_coverage,
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

summarise_susceptible <- function(df) {
    disease <- str_to_title(unique(df$vaccine_type))
    scenario_map1 <- c(central = "Central", pessimistic = "Pessimistic", optimistic = "Optimistic")
    scenario_map <- c(Central = 2, Pessimistic = 1, Optimistic = 3)

    cols <- ghibli::ghibli_palette("PonyoMedium", 3, type = "discrete")
    names(cols) <- names(scenario_map)
    
    #plot of total immunity
    total_plot <- df %>%
        group_by(scenario, date) %>%
        summarise(
            infection = 1 - (sum(immune)/sum(population)),
            disease = 1 - (sum(immune_disease)/sum(population)),
            .groups = "drop"
        ) %>%
        mutate(
            `Scenario:` = scenario_map1[scenario]
        ) %>%
        pivot_longer(
            cols = c(infection, disease),
            names_to = "type",
            values_to = "immunity"
        ) %>%
        mutate(
            type = str_to_title(type),
            type = factor(type, levels = c("Infection", "Disease"))
        ) %>%
        ggplot(aes(x = date, y = immunity, color = `Scenario:`, linetype = type)) +
            geom_line() +
            labs(x = "Date", y = "Susceptible %\n(Whole Population)", title = disease, type = "Protection against:") +
            ggpubr::theme_pubclean() +
            scale_y_continuous(labels = scales::percent) +
            scale_x_date(date_breaks = "1 month", date_labels = "%m-%Y") +
            scale_color_manual(values = cols) +
            coord_cartesian(ylim = c(0, 1)) +
            theme(plot.title = element_text(hjust = 0.5))

    text_df <- df %>%
        group_by(scenario, age_group, date) %>%
        summarise(
            infection = 1 - (sum(immune)/sum(population)),
            disease = 1 - (sum(immune_disease)/sum(population)),
            .groups = "drop_last"
        ) %>%
        summarise(
            infection = mean(infection),
            disease = mean(disease),
            .groups = "drop"
        ) %>%
        mutate(
            scenario = str_to_title(scenario),
            infection = scales::percent(infection, accuracy = 0.1),
            disease = scales::percent(disease, accuracy = 0.1)
        ) %>%
        mutate(
            x_pos = as.numeric(age_group),
            y_pos = scenario_map[scenario]
        )
    cols_2 <- as.list(cols)
    names(cols_2) <- scenario_map[names(cols_2)]
    cols_2$`4` <- "black"
    
    text_plot <- text_df %>%
        transmute(
            x_pos = x_pos, y_pos = y_pos, text = paste0(infection, "\n", disease)
        ) %>%
        rbind(
            tibble(
                x_pos = 0,
                y_pos = scenario_map,
                text = names(scenario_map)
            )
        ) %>% 
        rbind(
            tibble(
                x_pos = seq_along(unique(df$age_group)),
                y_pos = max(text_df$y_pos) + 1,
                text = unique(df$age_group)
            )
        ) %>%
        rbind(
            tibble(
                x_pos = 0,
                y_pos = max(text_df$y_pos) + 1,
                text = paste0("Infection %\nDisease %")
            )
        ) %>%
        mutate(
            x_pos = as.factor(x_pos),
            y_pos = as.factor(y_pos)
        ) %>%
        ggplot(aes(x = x_pos, y = y_pos, label = text, colour = as.character(y_pos))) +
            geom_text() +
            labs(x = "", y = "", title = "Average Susceptibility by age-group") +
            ggpubr::theme_pubclean() +
            theme(
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                panel.grid.major.y  = element_blank(),
                plot.title = element_text(hjust = 0.5),
                legend.position = "none"
            ) +
            scale_color_manual(values = cols_2)
    ggpubr::ggarrange(
        total_plot, text_plot,
        ncol = 1, nrow = 2,
        heights = c(4, 2)
    )
}
