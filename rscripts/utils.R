# some helper functions in this script
library(tidyverse)

run_simulation_from_python <- function(model_file, parameters, ignore_stdout = TRUE, return_res = FALSE, print_call_string = FALSE) {
  call_string <- str_interp('python3 main.py -pf ${model_file}')
  if (!missing(parameters)) {
    flags <- paste0('-', names(parameters))
    call_string <- str_interp('${call_string} ${paste(flags, format(parameters, scientific = FALSE), collapse = " ")}')
  }

  # calling the python script for running simulations
  if (print_call_string) {
    cat(call_string)
    cat('\n')
  } else {
    system(call_string, ignore.stdout = ignore_stdout)

    if (return_res) {
      return(RJSONIO::fromJSON(parameters$of, simplifyVector = FALSE, simplifyDataFrame = FALSE, flatten = FALSE, simplifyMatrix = FALSE))
    }
  }
}

combine_batch_files <- function(input_jsons, output_dat, remove_tmp = FALSE) {
  filenames <- str_split(input_jsons, "/") %>% map_chr(last)
  jsons <- lapply(input_jsons, RJSONIO::fromJSON) %>% set_names(filenames)
  saveRDS(jsons, file = output_dat)
  if (remove_tmp) {
    file.remove(input_jsons)
  }
}

compute_and_save_exact_degree_distribution <- function(R0, beta, gamma, k, return_res = FALSE) {
  n_max <- 10000
  infection_probability <- beta/(beta + gamma)
  file_name <- str_interp('data/degree_distribution/new_params_R0_${R0}_beta_${beta}_gamma_${gamma}_k_${k}.csv')
  if (!file.exists(file_name)) {
    tryCatch({
    new_ave <- 1/sum(map_dbl(seq(100000), ~dnbinom(.x-1, mu = R0/infection_probability, size = k)/.x))
    map_df(seq(n_max), ~list(degree = .x, prob = dnbinom(.x-1, mu = R0/infection_probability, size = k)*new_ave/.x)) %>%
      fwrite(file = file_name)
    }, error = function(e) browser())
  }
  if (return_res) {
    return(fread(file_name))
  }
}

aggregate_top_5_comminuty_infection <- function(top_5_communities_infection) {
  if (is.list(top_5_communities_infection)) {
    top_5_communities_infection <- map(top_5_communities_infection, ~.x[seq(min(5, length(.x)))]) %>% unlist()
  } else {
    top_5_communities_infection <- top_5_communities_infection[,seq(min(5, ncol(top_5_communities_infection)))]
  }
  mean(top_5_communities_infection[top_5_communities_infection > 0])
}

process_alert_levels <- function(daily_counts) {
  if (length(daily_counts) < 1) return(NA)

  # get daily counts for the past 14 days
  daily_counts <- daily_counts[,2]
  imap_int(daily_counts, function(x, i) as.integer(sum(daily_counts[max(0, i-13):i])))
}

compute_average_ratio <- function(.y, num_tracers, N = 1e5, type = 'random', sim, return_daily = FALSE) {

  if (length(N) == 1) {
    N <- rep(N, length(.y$in_positive_induced_daily))
  }

  daily_in_infected <- map(seq_along(.y$in_positive_induced_daily), function(j) {

    people_in_hospitalized <- (sim$daily_counts$Ho[[j]][,2] %>% head(-1) %>% cumsum)
    if ('in_ctd_daily_start_of_day' %in% names(.y)) {
      people_who_are_in_ctd <- .y$in_ctd_daily_start_of_day[[j]] - people_in_hospitalized
    } else {
      people_who_are_in_ctd <- .y$in_quarantine_daily_start_of_day[[j]] - people_in_hospitalized # when pQ=1
    }

    people_will_not_be_selected_for_test_total_daily_start_of_day <- people_who_are_in_ctd + people_in_hospitalized


    (.y$in_positive_induced_daily[[j]] + people_in_hospitalized)/(N[j])
  })

  daily_tested_positive <- map(seq_along(.y$in_positive_induced_daily), function(j) {
    people_in_hospitalized <- (sim$daily_counts$Ho[[j]][,2] %>% head(-1) %>% cumsum)
    if ('in_ctd_daily_start_of_day' %in% names(.y)) {
      people_will_not_be_selected_for_test_total_daily_start_of_day <- .y$in_ctd_daily_start_of_day[[j]]
    } else {
      people_will_not_be_selected_for_test_total_daily_start_of_day <- .y$in_quarantine_daily_start_of_day[[j]] # when pQ=1
    }

    if (length(.y$tested_positive_daily_end_of_day[[j]]) > 0) {
      if (type == 'random') {
        tested_positive_daily_end_of_day <- .y$tested_positive_daily_end_of_day[[j]][,2]
      } else if (type == 'contact') {
        tested_positive_daily_end_of_day <- .y$tested_positive_daily_end_of_day[[j]][,1]
      } else {
        if (is.matrix(.y$tested_positive_daily_end_of_day[[j]])) {
          tested_positive_daily_end_of_day <- .y$tested_positive_daily_end_of_day[[j]][,1] + .y$tested_positive_daily_end_of_day[[j]][,2]
        } else {
          tested_positive_daily_end_of_day <- .y$tested_positive_daily_end_of_day[[j]]
        }
      }
      if (length(num_tracers) == 1) {
        num_tracers_count <- num_tracers
      } else {
        num_tracers_count <- (num_tracers[[j]][, 1] + num_tracers[[j]][, 2])
      }
      ((tested_positive_daily_end_of_day/num_tracers_count) * (N[j]-people_will_not_be_selected_for_test_total_daily_start_of_day)+people_will_not_be_selected_for_test_total_daily_start_of_day)/N[j]
    } else {
      numeric()
    }
  })

  if (return_daily) return(daily_tested_positive)
  list(mean(discard(map2_dbl(daily_in_infected, daily_tested_positive, ~cor(.x[.x > 0], .y[.x > 0])), is.nan), na.rm = TRUE),
       mean(discard(map2_dbl(daily_in_infected, daily_tested_positive, ~mean(.y[.x > 0]/.x[.x > 0])), is.nan)),
       mean(discard(map2_dbl(daily_in_infected, daily_tested_positive, ~mean(.y[.x > 0])), is.nan)),
       mean(discard(map2_dbl(daily_in_infected, daily_tested_positive, ~mean(.x[.x > 0])), is.nan)))
}

make_label <- function(label, value, decreasing = FALSE) {
  factor(paste(label, value), levels = paste(label, sort(unique(value), decreasing = decreasing)))
}

convert_to_numeric <- function(vector) {
  suppressWarnings({
    converted_vec <- as.numeric(vector)
  })
  if (any(is.na(converted_vec))) {
    return(vector)
  } else {
    return(converted_vec)
  }
}