library(tidyverse)
library(parallel)
library(jsonlite)
library(data.table)
library(see)
library(latex2exp)
source('rscripts/utils.R')

####################################### loading data #####################
# chosen_disease <- 'SIR_with_Ho'
chosen_disease <- 'SEIR_with_Ho'
# convert strategy names to their correct once in paper
swap_tracing_types <- list('forward tracing', 'back tracing', 'global oracle', 'oracle tracer', 'random', 'oracle backward tracer') %>% #
  set_names(c('forward', 'backward', 'globaloracle', 'oracleTracer', 'random', 'backOracleTracer'))
if (chosen_disease == 'SIR_with_Ho') {
  simulation_parameter_names <- c('disease', 'N', 'num_tracers', 'type', 'init_infects',
                                'rate_IS2II' ,'rate_I2R', 'R0', 'dispersion', 'HoP',
                                'rate_I2Ho', 'org_rate_I2R', 'pCT')
} else {
  simulation_parameter_names <- c('disease', 'N', 'num_tracers', 'type', 'init_infects',
                                  'rate_IS2IE' ,'rate_I2R', 'rate_E2I', 'R0', 'dispersion', 'HoP',
                                  'rate_I2Ho', 'org_rate_I2R', 'pCT')
}

output_path <- str_interp('data/${chosen_disease}_exhaustive_parameter')

used_params <- mclapply(list.files(output_path), function(file) {
  readRDS(file.path(output_path, file)) %>%
    imap(function(sim, json_name) {
      sim <- process_sim_data_format(sim)
      # obtain the parameters from the name
      ret <- str_split(str_replace(json_name, '.json', ''), '-')[[1]] %>%
        as.list() %>%
        map(convert_to_numeric) %>%
        magrittr::set_names(simulation_parameter_names)
      ret[['type']] <- swap_tracing_types[[ret[['type']]]]
      ret2_names <- c('cor_random_daily_possitive', 'average_random_test_ratio', 'average_daily_tested_positive', 'average_daily_in_infected')

      if (ret[['num_tracers']] > 0 & ret[['type']] %in% c('forward tracing', 'back tracing', 'random')) {
        # compute the correlations
        ret2 <- compute_average_ratio(sim$tests_correlation_records, num_tracers = {if (all(is.null(sim$tracers_per_day))) ret[['num_tracers']] else sim$tracers_per_day}, N = sim$network_sizes, type = '', sim = sim) %>%
          magrittr::set_names(ret2_names)
      } else {
        ret2 <- as.list(rep(NA, length(ret2_names))) %>% magrittr::set_names(ret2_names)
      }

      if (chosen_disease == 'SEIR_with_Ho') {
        current_infected <- map2(sim$daily_counts$I, sim$daily_counts$E, `+`)
      } else {
        current_infected <- sim$daily_counts$I
      }
      mean_ignore_inf <- function(x) mean(discard(x, is.infinite))
      ret3 <- list(
        total_infected = mean(sim$infected_totals / sim$network_sizes),
        ave_days2end = mean(sim$days2end),
        top_5_communities = if (!is.null(sim$top_5_communities_infection)) aggregate_top_5_comminuty_infection(sim$top_5_communities_infection) else NA_real_,
        highest_alert_levels = mean_ignore_inf(map_int(current_infected, function(y) max(process_alert_levels(y), na.rm = TRUE)))
      )

      c(ret, ret2, ret3)
    })
}, mc.cores = 30) %>%
  map(~rbindlist(file, fill=TRUE)) %>%
  rbindlist(fill=TRUE) %>%
  mutate(dispersion = if_else(dispersion == 'NA', 'ER', as.character(dispersion)))

used_params[num_tracers > 0 & type %in% c('forward tracing', 'back tracing', 'random') & is.na(cor_random_daily_possitive), cor_random_daily_possitive := 0]
used_params[type == 'random', pCT:=0]

# for num_tracers=0, we may need to duplicate the results
if (length(unique(used_params[num_tracers==0]$type)) == 1) {
  used_params <- rbind(
    used_params[num_tracers!=0],
    crossing(used_params[num_tracers==0, -'type'], tibble(type = unique(used_params$type)))
  )
}


## get the results and rename columnes with the right names
if (chosen_disease == 'SIR_with_Ho') {
  res <- used_params %>%
    select(`Tracing type` = type, `Initial infections` = init_infects, beta = rate_IS2II, gamma = org_rate_I2R, R0, k = dispersion, `Probability of hopitalization` = HoP, `Average final infectioned population` = total_infected, `Average days to end` = ave_days2end, `Average final infection proportion of top 5 largest communities` = top_5_communities, `Average daily tested positive rates` = average_daily_tested_positive, `Average daily population proportion in infection statuses` = average_daily_in_infected, `Correlation between infection and positive each day` = cor_random_daily_possitive, `Highest average positive counts` = highest_alert_levels)
} else {
  res <- used_params %>%
    select(`Tracing type` = type, `Initial infections` = init_infects, beta = rate_IS2IE, gamma = org_rate_I2R, R0, k = dispersion, `Probability of hopitalization` = HoP, `Average final infectioned population` = total_infected, `Average days to end` = ave_days2end, `Average final infection proportion of top 5 largest communities` = top_5_communities, `Average daily tested positive rates` = average_daily_tested_positive, `Average daily population proportion in infection statuses` = average_daily_in_infected, `Correlation between infection and positive each day` = cor_random_daily_possitive, `Highest average positive counts` = highest_alert_levels)
}