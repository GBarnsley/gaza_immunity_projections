fois <- setdiff(unique(unlist(foi)), 0)
fois <- exp(seq(log(min(fois)), log(max(fois)), length.out = 5))
fois <- 0.1*10^-(0:5)
diseases <- names(keep(status_quo_parameters$force_of_infection, ~tail(.x, 1) >0))

sensitivity_parameters <- status_quo_parameters
sensitivity_parameters <- map(sensitivity_parameters, function(x){
    if(all(diseases %in% names(x))){
        x[diseases]
    } else {
        x
    }
})

under_6 <- which(cumsum(age_group_sizes)/365 == 6)

fois <- c(fois, NA)

results <- map_dfr(fois, function(new_foi){
    if(is.na(new_foi)){
        sensitivity_parameters$force_of_infection <- status_quo_parameters$force_of_infection[diseases]
        new_foi <- "Model Value"
    } else {
        sensitivity_parameters$force_of_infection <- map(sensitivity_parameters$force_of_infection, ~new_foi)
    }
    output <- do.call(project_point_estimate, sensitivity_parameters)
    #limit to children
    output$projections %>%
        filter(age_group <= under_6) %>%
        group_by(t, vaccine_type) %>%
        summarise(
            immunity = sum(Immune)/sum(Population),
            foi = as.character(new_foi),
            .groups = "drop"
        )
})

colour_scale <- c(viridis::viridis(length(fois)-1, "Blues"), "black")

names(colour_scale) <- c(fois[-length(fois)], "Model Value")

foi_plot <- results %>%
    mutate(
        disease = paste0(str_to_title(vaccine_type), ": ", signif(unlist(status_quo_parameters$force_of_infection[vaccine_type]), 3))
    ) %>%
    ggplot(aes(x = t, y = immunity, colour = foi, linetype = foi == "Model Value")) +
        geom_line(linewidth = 1) +
        facet_wrap(~disease, ncol = 2, scales = "free_y") +
        scale_y_continuous(labels = scales::percent) +
        scale_colour_manual(values = colour_scale) +
        ggpubr::theme_pubclean() +
        labs(x = "Time (days)", y = "Proportion immune (Under 6)", title = "Sensitivity to assumed FoI in the status_quo scenario", linetype = "Modelled Value:", colour = "Force of\nInfection") +
        theme(
            legend.position = "bottom"
        )

ggsave("plots/sensitivity_foi.pdf", foi_plot)

#sensitivity to assumed immunity
diseases <- names(c(moderate_background_immunity, high_background_immunity))
p_prior_immunity <- seq(0, 1, length.out = 10)

sensitivity_parameters <- status_quo_parameters
sensitivity_parameters <- map(sensitivity_parameters, function(x) {
    if(all(diseases %in% names(x))) {
        x[diseases]
    } else {
        x
    }
})

convert_to_percentages <- function(x) {
    paste0(as.character(round(x * 100)), "%")
}

p_prior_immunity <- c(p_prior_immunity, NA)

results <- map_dfr(p_prior_immunity, function(prior_immunity) {
    temp_pars <- sensitivity_parameters
    if(is.na(prior_immunity)){
        prior_immunity <- "Model Value"
    } else {
        total_pop <- sensitivity_parameters$S_0[[1]] + sensitivity_parameters$R_0[[1]]
        temp_pars$S_0 <- map(sensitivity_parameters$S_0, ~(1 - prior_immunity) * total_pop)
        temp_pars$R_0 <- map(sensitivity_parameters$S_0, ~prior_immunity * total_pop)
        prior_immunity <- convert_to_percentages(prior_immunity)
    }
    output <- do.call(project_point_estimate, temp_pars)
    #limit to children
    output$projections %>%
        group_by(t, vaccine_type) %>%
        summarise(
            immunity = sum(Immune)/sum(Population),
            prior_immunity = prior_immunity,
            .groups = "drop"
        )
})

results$prior_immunity

colour_scale <- c(viridis::viridis(length(p_prior_immunity)-1), "black")

names(colour_scale) <- c(convert_to_percentages(p_prior_immunity[-length(p_prior_immunity)]), "Model Value")

true_prop <- convert_to_percentages(if_else(diseases %in% moderate_background_immunity, p_moderate, p_high))
names(true_prop) <- diseases

prior_immunity_plot <- results %>%
    mutate(
        disease = paste0(str_to_title(vaccine_type), ": ", true_prop[vaccine_type])
    ) %>%
    ggplot(aes(x = t, y = immunity, colour = prior_immunity, linetype = prior_immunity == "Model Value")) +
        geom_line(linewidth = 1) +
        facet_wrap(~disease, ncol = 2, scales = "free_y") +
        scale_y_continuous(labels = scales::percent) +
        scale_colour_manual(values = colour_scale) +
        ggpubr::theme_pubclean() +
        labs(x = "Time (days)", y = "Proportion immune", title = "Sensitivity to assumed level of immunity in 2000", linetype = "Modelled Value:", colour = "Prior\nImmunity") +
        theme(
            legend.position = "bottom"
        )

ggsave("plots/sensitivity_prior_immunity.pdf", prior_immunity_plot)