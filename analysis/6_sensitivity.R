fois <- setdiff(unique(unlist(foi)), 0)
fois <- exp(seq(log(min(fois)), log(max(fois)), length.out = 5))
fois <- 0.1*10^-(0:5)
diseases <- names(keep(central_parameters$force_of_infection, ~tail(.x, 1) >0))

sensitivity_parameters <- central_parameters
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
        sensitivity_parameters$force_of_infection <- central_parameters$force_of_infection[diseases]
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

colour_scale <- c(RColorBrewer::brewer.pal(length(fois)-1, "Blues"), "black")

names(colour_scale) <- c(fois[-length(fois)], "Model Value")

foi_plot <- results %>%
    mutate(
        disease = paste0(str_to_title(vaccine_type), ": ", signif(unlist(central_parameters$force_of_infection[vaccine_type]), 3))
    ) %>%
    ggplot(aes(x = t, y = immunity, colour = foi, linetype = foi == "Model Value")) +
        geom_line(linewidth = 1) +
        facet_wrap(~disease, ncol = 2, scales = "free_y") +
        scale_y_continuous(labels = scales::percent) +
        scale_colour_manual(values = colour_scale) +
        ggpubr::theme_pubclean() +
        labs(x = "Time (days)", y = "Proportion immune (Under 6)", title = "Sensitivity to assumed FoI in the central scenario", linetype = "Modelled Value:", colour = "Force of\nInfection") +
        theme(
            legend.position = "bottom"
        )

ggsave("plots/sensitivity_foi.pdf", foi_plot)