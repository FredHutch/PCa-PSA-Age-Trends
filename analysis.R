################################################################################
# Analysis supporting "Trends in Age and Prostate-Specific Antigen at 
# Prostate Cancer Diagnosis between 2010 and 2019"
# Authors: Lukas Owens
#          Ojas Brahme 
#          Roman Gulati 
#          Ruth Etzioni
# Fred Hutchinson Cancer Center
# Description: This study examines trends in age and PSA at diagnosis among 
#     men diagnosed with prostate cancer between 2010 and 2019 in the SEER-17 
#     database. Empirical distributions of these two quantities are plotted and
#     quantile regressions of the 25th, 85th, 75th and 85th percentiles by year
#     and race/ethinic origin are estimated.
# Data: The SEER queries creating the study data are available in the `data` 
#     folder. The output of these queries must be saved in correspondingly 
#     named CSV files. These queries were created with SEER*Stat v8.4.3
################################################################################

########################################
# Set up
########################################

library(tidyverse)
library(ggrepel)

# Reference year for standardizing population
ref_year <- 2010

gg_save <- function(filepath, object, width = 500, height = 500) {
    png(filepath, width = width, height = height, res = 96)
    print(g)
    dev.off()
}

########################################
# Read SEER case data
########################################
parse_stage <- function(stage) {
    stage <- case_when(stage %in% c('Localized', 'Regional') ~ 'Localized / Regional',
                       stage == 'Distant' ~ 'Distant',
                       TRUE ~ NA)
    factor(stage, levels = c('Localized / Regional', 'Distant'))
}
age_groups <- function(age, width = 5, max_age = 85) {
    age_group <-  5*floor(age/5)
    ifelse(age >= max_age, paste0(max_age, '+'), str_glue('{age_group}-{age_group+4}'))
}
read_seer <- function() {
    dset <- read_csv('data/pca-case-listing.csv',
                     col_names =  c('year', 'race', 'age', 'stage', 'psa'),
                     col_types = 'iccc',
                     skip = 1)
    # Parse PSAs and keep missing value info
    dset <- dset |> mutate(psa_value = parse_number(psa))
    dset <- dset |> mutate(psa_str = ifelse(is.na(psa_value), psa, NA_character_))
    dset <- dset |> select(-psa)
    
    # Parse other fields
    dset <- dset |> mutate(age = as.integer(str_sub(age, 1, 2)))
    dset <- dset |> mutate(age_group = factor(age_groups(age)))
    dset <- dset |> mutate(stage = parse_stage(stage))
    
    # Parse race
    dset <- dset |> mutate(race = str_remove(race, ' \\(All Races\\)'),
                           race = str_replace(race, 'Non-Hispanic', 'NH'),
                           race = str_replace(race, 'American Indian/Alaska Native', 'AI/AN'),
                           race = str_replace(race, 'Asian or Pacific Islander', 'API'),
                           race = ifelse(race == 'NH Unknown Race', NA, race),
                           race = fct_infreq(factor(race)))
    # Cut off at age 85
    dset <- dset |> filter(age <= 85)
    
    dset <- dset |> select(year, age, race, age_group, stage, psa_value, psa_str)
    dset
}
case_data <- read_seer()
case_data

########################################
# US population data for standardization
########################################
read_pop <- function() {
    dset <- read_csv('data/us-population.csv',
                     col_names = c('age', 'race', 'year', 'count'),
                     skip = 1)
    dset <- dset |> mutate(age = as.integer(str_sub(age, 1, 2)))
    dset <- dset |> mutate(age_group = factor(age_groups(age)))
    dset <- dset |> filter(str_length(year) == 4)
    dset <- dset |> mutate(year = as.integer(year))
    dset
}
pop_data_race <- read_pop()
pop_data <- pop_data_race |> group_by(age, year) |> summarize(count = sum(count))
ref_data <- pop_data |> filter(year == ref_year) |> select(age, count)

########################################
# Data quality checks
########################################

# Check missing PSA data by year
g <- case_data |> mutate(ind = as.integer(is.na(psa_value))) |> 
    group_by(year) |> 
    summarize(missing = sum(ind), n = n()) |> 
    mutate(percent_missing = missing / n) |> 
    ggplot(aes(x = year, y = percent_missing))
g <- g + geom_col(fill = 'lightgray', color = 'black')
g <- g + scale_y_continuous(limits = c(0,1), labels = scales::percent)
g <- g + scale_x_continuous(breaks = seq(2010, 2019, by = 2))
g <- g + labs(x = 'Year', y = 'Percent missing PSA')
g <- g + theme_bw()
g


########################################
# Baseline characteristics
########################################
case_data |> select(age, race, stage, psa_value) |> 
    gtsummary::tbl_summary() |>
    gtsummary::as_flex_table() |> 
    flextable::save_as_docx(path = 'output/table-1.docx')


########################################
# Joint contour plot
########################################
# Plot parameters
contour_plots <- function(case_data,
                          n_bins = 5,
                          smoothing = 2) {
    
    # Plot data
    gset <- case_data |> filter(!is.na(psa_value), !is.na(race))
    gset <- gset |> mutate(year_start = (year == ref_year))
    
    # Set bandwidth
    h <- c(MASS::bandwidth.nrd(gset$age), MASS::bandwidth.nrd(gset$psa_value))
    h <- smoothing*h # Increase to smooth contours
    
    # 2010 distribution for reference
    oset <- gset |> filter(year == ref_year) |> select(-year) 
    oset <- cross_join(oset, tibble(year = 2011:2019))
    
    # Calculate mode
    dens <- function(x, y, x0, y0, h) {
        n <- length(x)
        if (missing(h)) h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y)) else h
        h <- h/4
        ax <- (x - x0) / h[1L]
        ay <- (y - y0) / h[2L]
        sum(dnorm(ax) * dnorm(ay))/ (n * h[1L] * h[2L])
    }
    f <- function(data) {
        dens_fn <- function(v) dens(data$age, data$psa_value, v[1], v[2], h)
        o <- optim(c(70, 7.5), dens_fn, control = list(fnscale = -1))
        list(age = o$par[1],
             psa_value = o$par[2])
    }
    tset <- gset |> select(year, age, psa_value) |> nest(data = c(age, psa_value))
    tset <- tset |> mutate(mode = map(data, f))
    tset <- tset |> unnest_wider(mode) |> arrange(year) |> select(-data)
    tset <- tset |> mutate(label = str_glue('Mode Age: {scales::number(age, 0.1)} years\nMode PSA: {scales::number(psa_value, 0.1)} ng/mL'))
    
    # Annotations
    annotate_data <- tibble(
        year = 2011,
        text_x = 65,
        text_y = 10.1,
        text_label = '2010 contour',
    )
    
    # Plot
    g <- ggplot(gset)
    g <- g + geom_density_2d_filled(aes(x = age, y = psa_value), alpha = 0.85, color = 'black',
                                    h = h, contour_var = 'ndensity', bins = n_bins)
    g <- g + geom_density_2d(data = oset, aes(x = age, y = psa_value), contour_var = 'ndensity', 
                             bins = n_bins, alpha = 1, color = 'red', h = h)
    g <- g + geom_text_repel(data = annotate_data, aes(x = text_x, y = text_y, label = text_label),
                             color = 'red', size = 3, nudge_x = -5, nudge_y = 1.15)
    g <- g + geom_text(data = tset, aes(x = 46, y = 0.3, label = label), size = 3, hjust = 'left', vjust = 'bottom')
    g <- g + geom_hline(aes(yintercept = 12), color = 'white')
    g <- g + facet_wrap(~year, nrow = 2)
    g <- g + scale_x_continuous(expand = c(0,0), breaks = seq(45,85, by = 10), limits = c(45, 85))
    g <- g + scale_y_continuous(expand = c(0,0), limits = c(0, 12))
    g <- g + scale_fill_manual(values = c('#FFFFFF', '#A8E1BCFF', '#3AAEADFF', '#3671A0FF', '#3E356BFF'))
    g <- g + scale_color_manual(values = c('black', 'red'))
    g <- g + theme_bw()
    g <- g + theme(panel.spacing.x = unit(1, 'lines'),
                   panel.grid.major = element_blank(),
                   axis.ticks.length = unit(0.2, 'cm'),
                   legend.position = 'bottom',
                   strip.background = element_blank(),
                   legend.key = element_rect(color = 'black'))
    g <- g + guides(color = 'none', fill = 'none')
    g <- g + labs(x = 'Age at Diagnosis (Years)', y = 'PSA at Diagnosis (ng/mL)', fill = 'Density')
    g
    
}
g <- contour_plots(case_data)
gg_save('output/figure-1.png', g, 800, 500)

##############################
# Plot age and PSA at diagnosis
##############################
regression_plots <- function(case_data,
                             ref_data) {
    # Age
    gset <- case_data |> filter(!is.na(psa_value), !is.na(race))
    gset <- rbind(
        gset |> mutate(race = as.character(race)),
        gset |> mutate(race = 'All Races')
    )
    gset <- gset |> mutate(race = factor(race, levels = c('All Races', levels(case_data$race))))
    gset <- gset |> left_join(ref_data)
    gset <- gset |> group_by(year, race) |> summarize(q25 = modi::weighted.quantile(age, count, prob = 0.25),
                                                      q50 = modi::weighted.quantile(age, count, prob = 0.50),
                                                      q75 = modi::weighted.quantile(age, count, prob = 0.75),
                                                      q85 = modi::weighted.quantile(age, count, prob = 0.85))
    gset <- gset |> pivot_longer(cols = starts_with('q'), names_prefix = 'q')
    g1 <- ggplot(gset)
    g1 <- g1 + geom_smooth(aes(x = year, y = value, group = name, color = name), method = 'lm', 
                           se = FALSE, linewidth = 0.2)
    g1 <- g1 + geom_point(aes(x = year, y = value, color = name))
    g1 <- g1 + facet_wrap(~race, nrow = 1)
    g1 <- g1 + theme_bw()
    g1 <- g1 + scale_x_continuous(breaks = seq(2010, 2018, by = 2))
    g1 <- g1 + scale_y_continuous(limits = c(50, 75))
    g1 <- g1 + scale_color_viridis_d(option = 'D', end = 0.75)
    g1 <- g1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     panel.grid = element_blank(),
                     legend.position = 'none',
                     axis.ticks.length = unit(0.2, 'cm'),
                     strip.background = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(colour = 'black'))
    g1 <- g1 + labs(x = '', y = 'Age at Diagnosis (Years)', color = 'Percentile', shape = 'Percentile')
    g1
    
    # PSA
    gset <- case_data |> filter(!is.na(psa_value), !is.na(race))
    gset <- rbind(
        gset |> mutate(race = as.character(race)),
        gset |> mutate(race = 'All Races')
    )
    gset <- gset |> mutate(race = factor(race, levels = c('All Races', levels(case_data$race))))
    gset <- gset |> left_join(ref_data)
    gset <- gset |> group_by(year, race) |> summarize(q25 = modi::weighted.quantile(psa_value, count, prob = 0.25),
                                                      q50 = modi::weighted.quantile(psa_value, count, prob = 0.50),
                                                      q75 = modi::weighted.quantile(psa_value, count, prob = 0.75),
                                                      q85 = modi::weighted.quantile(psa_value, count, prob = 0.85))
    gset <- gset |> pivot_longer(cols = starts_with('q'), names_prefix = 'q')
    gset <- gset |> mutate(name = paste0(name, 'th'))
    g2 <- ggplot(gset)
    g2 <- g2 + geom_smooth(aes(x = year, y = value, group = name, color = name), method = 'lm', 
                           se = FALSE, linewidth = 0.2)
    g2 <- g2 + geom_point(aes(x = year, y = value, color = name))
    g2 <- g2 + facet_wrap(~race, nrow = 1)
    g2 <- g2 + theme_bw()
    g2 <- g2 + scale_x_continuous(breaks = seq(2010, 2018, by = 2))
    g2 <- g2 + scale_y_continuous(limits = c(0, NA))
    g2 <- g2 + scale_color_viridis_d(option = 'D', end = 0.75)
    g2 <- g2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     panel.grid = element_blank(),
                     legend.position = 'bottom',
                     legend.box.spacing = unit(0, 'pt'),
                     axis.ticks.length = unit(0.2, 'cm'),
                     strip.background = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line(colour = 'black'))
    g2 <- g2 + labs(x = 'Year of Diagnosis', y = 'PSA at Diagnosis (ng/mL)', color = 'Percentile', shape = 'Percentile')
    g2
    
    g <- cowplot::plot_grid(g1, g2, labels = c('A', 'B'), ncol = 1, rel_heights = c(1, 1.15))
    g
    
}
g <- regression_plots(case_data, ref_data)
gg_save('output/figure-2.png', g, 725, 550)

########################################
# Trends in age and PSA at diagnosis by race
########################################
ci <- function(mean, std_error, alpha=0.05) {
    str_glue('{scales::number(mean+qnorm(alpha/2)*std_error, accuracy = 0.01)}, {scales::number(mean+qnorm(1-alpha/2)*std_error, accuracy = 0.01)}')
}

dset <- case_data |> left_join(ref_data)
dset <- dset |> filter(!is.na(race), !is.na(psa_value))
dset <- dset |> mutate(year = year-ref_year)

quantile_regression <- function(formula, savepath = NULL) {
    quantiles <- c(0.25, 0.50, 0.75, 0.85)
    fit <- quantreg::rq(formula,
                        tau = quantiles,
                        data = dset,
                        weights = dset$count,
                        method = 'pfn')
    s <- summary(fit)
    s <- do.call(rbind, lapply(1:length(quantiles), 
                               \(x) as_tibble(s[[x]]$coefficients, rownames = 'coef') |> mutate(tau = x)))
    names(s) <- c('term', 'estimate', 'std.error', 't_value', 'p', 'tau')
    s <- s |> mutate(ci = ci(estimate, std.error),
                     estimate = scales::number(estimate, accuracy = 0.01),
                     p = scales::number(p, accuracy = 0.001),
                     quantity = as.character(formula)[2])
    if (!is.null(savepath)) {
        s |> select(quantity, term, estimate, ci, p) |> 
            flextable::as_flextable(max_row = Inf) |> 
            flextable::save_as_docx(path = savepath)
    }
    s
}
quantile_regression(formula = age ~ year, 
                    savepath = 'output/table-1a.docx')
quantile_regression(formula = psa_value ~ year,  
                    savepath = 'output/table-2a.docx')
quantile_regression(formula = psa_value ~ year + race + year:race,  
                    savepath = 'output/table-2b.docx')
