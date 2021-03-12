library(tidyverse)
library(parallel)
library(jsonlite)
library(data.table)
source('rscripts/utils.R')

seed <- 190
N <- 50000
simIts <- 30
nxIts <- 50
jobs <- 1
batch_num <- c(87)
chosen_disease <- 'SIR_with_Ho'
# chosen_disease <- 'SEIR_with_Ho'

pQ <- 1 # probability of quarantining
pRT <- 1 # probability of random testing
RT_num <- 0 # reserved number of random tests
num_nodes <- 5000 # reserved inital infected nodes for analyzing the node properties


# cache all simulations results
output_path <- str_interp('data/${chosen_disease}_exhaustive_parameter')
dir.create(output_path)

source('rscripts/01_simulation_parameters.R')

varying_params_list <- map(varying_params_list, function(varying_params) {
  varying_params %>%
    mutate(rate_I2Ho = rate_I2R * (HoP), org_rate_I2R = rate_I2R, rate_I2R = rate_I2R * (1 - HoP), num_tracers_perc = round(num_tracers_perc * N)) %>%
    rename(num_tracers = num_tracers_perc) %>%
    mutate(pCT = map_dbl(num_tracers, ~ifelse(.x==0, 1,1 - min(round(0.5 * .x), RT_num) / .x)))
})

mclapply(seq(jobs), function(i) {
    mclapply(split(seq(sum(batch_num)), seq(jobs))[[i]], function(total_index) {
    index_set <- c(0, cumsum(batch_num))
    which_batch <- findInterval(total_index-1, index_set)
    varying_params <- varying_params_list[[which_batch]]

    batch <- split(seq(nrow(varying_params)), seq(batch_num[which_batch]))[[total_index - index_set[which_batch]]]
    print(paste0('running on batch: ', total_index, '; there are ', length(batch), ' parameter sets to test in this batch.'))
    print(batch)
    batch_output_path <- file.path(output_path, paste0('batch-', total_index, '-new-set-selected.rds'))
    if (!file.exists(batch_output_path)) {
        batch_output_files <- mclapply(batch, function(k) {
            output_file <- str_interp('tmp_disease_param_test/${chosen_disease}-${N}-${paste(varying_params[k,], collapse = "-")}.json')
            print(output_file)
            if (!file.exists(output_file)) {
                parameters <- list(
                    num_tracers = varying_params$num_tracers[k],
                    nI = varying_params$init_infects[k],
                    oracle_type = varying_params$type[k],
                    rate_I2R = varying_params$rate_I2R[k],
                    rate_I2Ho = varying_params$rate_I2Ho[k],
                    pCT = varying_params$pCT[k],
                    seed = seed,
                    of = output_file,
                    simIts = simIts,
                    nxIts = nxIts,
                    pQ = pQ,
                    pRT = pRT,
                    N = N,
                    parallel = 1,
                    num_nodes = num_nodes
                )
                beta <- ifelse(chosen_disease == 'SIR_with_Ho', varying_params$rate_IS2II[k], varying_params$rate_IS2IE[k])
                if (!is.na(varying_params$dispersion[k])) {
                    compute_and_save_exact_degree_distribution(varying_params$R0[k], beta, varying_params$rate_I2R[k]+varying_params$rate_I2Ho[k], varying_params$dispersion[k])
                    parameters <- c(parameters, list(nt = "nbinomial",
                                                    np = str_interp('\'{"dist_file": "data/degree_distribution/new_params_R0_${varying_params$R0[k]}_beta_${beta}_gamma_${varying_params$rate_I2R[k] + varying_params$rate_I2Ho[k]}_k_${varying_params$dispersion[k]}.csv"}\'')))
                } else {
                    network_params <- varying_params$R0[k]/(beta / (beta + varying_params$rate_I2R[k]))
                    parameters <- c(parameters, list(nt = "er",
                                                    np = str_interp('\'{"p": ${network_params/(N-1)}}\'')))
                }
                if (chosen_disease == 'SIR_with_Ho') {
                    parameters <- c(parameters, list(rate_IS2II = varying_params$rate_IS2II[k]))
                } else {
                    parameters <- c(parameters, list(rate_IS2IE = varying_params$rate_IS2IE[k],
                                                rate_E2I = varying_params$rate_E2I[k]))
                }
                run_simulation_from_python(model_file = paste0('parameters/abstract/', chosen_disease, '.json'),
                                        parameters = parameters,
                                        ignore_stdout = FALSE)
            }
            return(output_file)
        }, mc.cores = 1) %>% flatten_chr()
        combine_batch_files(batch_output_files, output_dat = batch_output_path, remove_tmp = TRUE)
    }
    }, mc.cores = 20)
}, mc.cores = 10)
