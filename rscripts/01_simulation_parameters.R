params_SIR <- rbind(
    crossing(tibble(num_tracers_perc = c(0.0001, 0.001, 0.01, 0.1)),
             tibble(type = c('random', 'forward', 'backward', 'globaloracle', 'oracleTracer', 'backOracleTracer')),
             tibble(init_infects = c(10, 50, 100)),
             tibble(rate_IS2II = c(0.2, 0.4, 0.6, 1)),
             tibble(rate_I2R = c(0.05, 0.25, 0.5, 1)),
             tibble(R0 = c(1, 1.5, 2, 2.5, 3, 3.5, 10)),
             tibble(dispersion = c(0.1, 0.5, 1, 1.5, NA)), # NA is ER networks
             tibble(HoP = c(0, 0.05, 0.1, 0.3))),
    crossing(tibble(num_tracers_perc = c(0)),
             tibble(type = c('random')),
             tibble(init_infects = c(10, 50, 100)),
             tibble(rate_IS2II = c(0.2, 0.4, 0.6, 1)),
             tibble(rate_I2R = c(0.05, 0.25, 0.5, 1)),
             tibble(R0 = c(1, 1.5, 2, 2.5, 3, 3.5, 10)),
             tibble(dispersion = c(0.1, 0.5, 1, 1.5, NA)), # NA is ER networks
             tibble(HoP = c(0, 0.05, 0.1, 0.3)))
)

params_SEIR <- rbind(
    crossing(tibble(num_tracers_perc = c(0.0001, 0.001, 0.01, 0.1)),
             tibble(type = c('random', 'forward', 'backward', 'globaloracle', 'oracleTracer', 'backOracleTracer')),
             tibble(init_infects = c(10, 50, 100)),
             tibble(rate_IS2IE = c(0.2, 0.4, 0.6, 1)),
             tibble(rate_I2R = c(0.05, 0.25, 0.5, 1)),
             tibble(rate_E2I = c(0.2)),
             tibble(R0 = c(1, 1.5, 2, 2.5, 3, 3.5, 10)),
             tibble(dispersion = c(0.1, 0.5, 1, 1.5, NA)), # NA is ER networks
             tibble(HoP = c(0, 0.05, 0.1, 0.3))),
    crossing(tibble(num_tracers_perc = c(0)),
             tibble(type = c('random')),
             tibble(init_infects = c(10, 50, 100)),
             tibble(rate_IS2IE = c(0.2, 0.4, 0.6, 1)),
             tibble(rate_I2R = c(0.05, 0.25, 0.5, 1)),
             tibble(rate_E2I = c(0.2)),
             tibble(R0 = c(1, 1.5, 2, 2.5, 3, 3.5, 10)),
             tibble(dispersion = c(0.1, 0.5, 1, 1.5, NA)), # NA is ER networks
             tibble(HoP = c(0, 0.05, 0.1, 0.3)))
)


varying_params_list <- list(params_SIR, params_SEIR)