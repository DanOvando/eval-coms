



# setup -------------------------------------------------------------------

library(tidyverse)

library(taxize)

library(countrycode)

library(here)

library(glue)

library(sraplus)

library(vip)

library(ramlegacy)

library(ranger)

library(xgboost)

library(tidymodels)

library(foreach)

library(doParallel)

library(ggtext)

library(FishLife)

library(gt)

library(patchwork)

tune_com <- FALSE

run_loo <- FALSE

functions <- list.files(here::here("R"))

purrr::walk(functions, ~ source(here::here("R", .x)))


##### set options #####


results_name <- "v0.1"

results_path <- here::here("results", results_name)

if (!dir.exists(results_path)) {
  dir.create(results_path, recursive = TRUE)
}

results_description <-
  "testing"

write(results_description,
      file = here::here("results", results_name, "description.txt"))

cores <- 8

get_ram_data <- FALSE #pull RAM data from dropbox

min_years_catch <- 25 # minimum years of catch data to include

crazy_b <- 5 # maximum B/Bmsy to allow

crazy_u <- 10 # maximum U/Umsy to allow

catchability <- 1e-3 # survey catchability

plot_theme <- theme_minimal(base_size = 14)

theme_set(plot_theme)

# get ram data ------------------------------------------------------------

if (get_ram_data) {
  if (file.exists(here("data", "ram.zip")) == FALSE) {
    download.file(
      "https://zenodo.org/record/4824192/files/RAMLDB%20v4.495.zip?download=1",
      destfile = here("data", "ram.zip"),
      mode = "wb"
    )

    unzip(here("data", "ram.zip"), exdir = here("data", "ram")) # unzip = 'unzip' needed for windows
  }

  ram_files <-
    list.files(here("data", "ram", "R Data"), recursive = TRUE)

  ram_files <- ram_files[str_detect(ram_files, ".RData")]

  load(here("data", "ram", "R Data", ram_files[1]))

  # process ram data ------------------------------------------------------------

  stock <- stock %>%
    left_join(area, by = "areaid")
  # catches
  ram_catches <- tcbest.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_tibble() %>%
    gather(stockid, catch, -year)

  # B/Bmsy
  ram_b_v_bmsy <- divbpref.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    tibble() %>%
    gather(stockid, b_v_bmsy, -year)


  # U/Umsy
  ram_u_v_umsy <- divupref.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_tibble() %>%
    gather(stockid, u_v_umsy, -year)

  # Effort
  ram_effort <- effort.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_tibble() %>%
    gather(stockid, effort, -year)

  # biomass

  ram_total_biomass <- tbbest.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_tibble() %>%
    gather(stockid, total_biomass, -year)

  # ssb

  ram_ss_biomass <- ssb.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_tibble() %>%
    gather(stockid, ss_biomass, -year)


  ram_exp_rate <- ram_catches %>%
    left_join(ram_total_biomass, by = c("stockid", "year")) %>%
    mutate(exploitation_rate = catch / total_biomass) %>%
    select(-catch,-total_biomass)

  # put it together

  ram_data <- ram_catches %>%
    left_join(bioparams_values_views, by = "stockid") %>%
    left_join(ram_b_v_bmsy, by = c("stockid", "year")) %>%
    left_join(ram_u_v_umsy, by = c("stockid", "year")) %>%
    left_join(ram_exp_rate, by = c("stockid", "year")) %>%
    left_join(ram_effort, by = c("stockid", "year")) %>%
    left_join(ram_total_biomass, by = c("stockid", "year")) %>%
    left_join(ram_ss_biomass, by = c("stockid", "year")) %>%
    left_join(stock, by = "stockid") %>%
    select(stockid, scientificname, commonname, everything())


  # create new variables

  ram_data <- ram_data %>%
    mutate(tb_v_tb0 = total_biomass / TB0,
           ssb_v_ssb0 = ss_biomass / SSB0)

  # filter data -------------------------------------------------------------

  # for now, only include continuous catch series

  ram_data <- ram_data %>%
    filter(is.na(catch) == FALSE) %>%
    # filter(stockid == "ATBTUNAEATL") %>%
    group_by(stockid) %>%
    mutate(delta_year = year - lag(year)) %>%
    mutate(delta_year = case_when(year == min(year) ~ as.integer(1),
                                  TRUE ~ delta_year)) %>%
    mutate(missing_gaps = any(delta_year > 1)) %>%
    filter(missing_gaps == FALSE) %>%
    mutate(n_years = n_distinct(year)) %>%
    filter(n_years >= min_years_catch) %>%
    filter(all(b_v_bmsy < crazy_b, na.rm = TRUE),
           all(u_v_umsy < crazy_u, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(stockid) %>%
    mutate(
      has_tb0 = !all(is.na(TB0)),
      has_tb = all(!is.na(total_biomass)),
      first_catch_year = year[which(catch > 0)[1]]
    ) %>%
    filter(year >= first_catch_year) %>%
    mutate(
      pchange_effort = lead(u_v_umsy) / (u_v_umsy + 1e-6),
      cs_effort = (u_v_umsy - mean(u_v_umsy)) / sd(u_v_umsy),
      index = total_biomass * catchability,
      approx_cpue = catch / (u_v_umsy / catchability + 1e-3),
      b_rel = dplyr::case_when(
        has_tb0 ~ total_biomass / max(TB0),
        has_tb ~ total_biomass / max(total_biomass),
        TRUE ~ b_v_bmsy / 2.5
      )
    ) %>%
    mutate(approx_cpue = pmin(quantile(approx_cpue, 0.9, na.rm = TRUE), approx_cpue)) %>%
    ungroup()

  ram_b_plot <- ram_data %>%
    ggplot(aes(x = year, y = b_v_bmsy)) +
    geom_bin2d() +
    scale_fill_viridis_c()

  kobe_panels <- ram_data %>%
    filter(year >= 1950) %>%
    mutate(year_block = plyr::round_any(year, 10, floor)) %>%
    ggplot(aes(x = b_v_bmsy, y = u_v_umsy)) +
    geom_bin2d(binwidth = c(0.5, 0.5)) +
    facet_wrap( ~ year_block) +
    scale_fill_viridis_c()

  # kobe_animation <- ram_data %>%
  #   ggplot(aes(x = b_v_bmsy, y = u_v_umsy)) +
  #   geom_bin2d() +
  #   transition_states(factor(year),
  #                     transition_length = 2,
  #                     state_length = 1) +
  #   scale_fill_viridis_c() +
  #   labs(title = "Year: {closest_state}")

  write_rds(ram_data, file = here("data", "ram", "ram-data.rds"))

} else {
  ram_data <- read_rds(file = here("data", "ram", "ram-data.rds"))


}



ram_data <- ram_data %>%
  rename(
    fao_area_code = primary_FAOarea,
    scientific_name = scientificname,
    common_name = commonname,
    capture = catch
  ) %>%
  mutate(
    macroid = paste(scientific_name, fao_area_code, sep = '_'),
    fao_area_code = as.integer(fao_area_code)
  )

ram_stocks <- ram_data %>%
  select(scientific_name, common_name, fao_area_code, macroid) %>%
  unique()


# demonstrating lack of prior updating from catch history shape ------------------------------------

make_catches <-
  function(form = "downhill",
           years = 25,
           mean_catch = 1000) {
    if (form == "downhill") {
      catches <- exp(1+-.1 * 0:(years - 1))

      catches <- catches / mean(catches) * mean_catch

    }

    if (form == "uphill") {
      catches <- exp(1+-.1 * 0:(years - 1))

      catches <- rev(catches / mean(catches) * mean_catch)

    }

    if (form == "peak") {
      catches <- exp(1  + .2 * 0:(years / 2))

      catches <- c(catches, rev(catches))

      catches <- catches[1:years]

      catches <-  catches / mean(catches) * mean_catch

    }


    if (form == "valley") {
      catches <- exp(1+-.2 * 0:(years / 2))

      catches <- c(catches, rev(catches))

      catches <- catches[1:years]

      catches <-  catches / mean(catches) * mean_catch

    }

    if (form == "random") {
      catches <- rlnorm(years)

      catches <- catches / mean(catches) * mean_catch


    }


    if (form == "walk") {
      catches <- rep(mean_catch, years)

      for (i in 2:years) {
        catches[i] <- catches[i - 1] * rlnorm(1, 0, .3)
      }

      catches <- catches / mean(catches) * mean_catch


    }


    return(data.frame(year = 1:years, catch = catches))

  }


catch_histories <-
  tibble(form = c("downhill", "uphill", "peak", "valley", "random", "walk")) %>%
  mutate(catches = map(form, make_catches))


# stop()
fit_example <-
  function(catches,
           initial_state = 1,
           initial_state_cv = .2,
           middle_state = 0.4,
           middle_state_cv = 0.2,
           terminal_state = 0.2,
           terminal_state_cv = 0.2) {
    com_driors <-
      format_driors(
        catch = catches$catch,
        years = catches$year,
        use_heuristics = FALSE,
        initial_state = initial_state,
        initial_state_cv = initial_state_cv,
        terminal_state = terminal_state,
        terminal_state_cv = terminal_state_cv,
        b_ref_type = "k",
        use_catch_prior = FALSE,
        taxa = "Gadus morhua",
        shape_prior_source = "thorson"
      )

    # plot_driors(com_driors)

    com_fit <-
      fit_sraplus(
        driors = com_driors,
        include_fit = TRUE,
        engine = "sir",
        draws = 1e6,
        tune_prior_predictive = FALSE
      )

    # plot_sraplus(com_fit)

    # plot_prior_posterior(com_fit, com_driors)
    #
    k <- com_fit$results %>%
      filter(variable == "k")

    # com_driors$years <- com_driors$years + 1950
    # cmsy_fit <- portedcmsy::funct_cmsy(
    #   catches = com_driors$catch,
    #   catch_years = com_driors$years,
    #   stock = "Gadus morhua" ,
    #   common_name = "cod",
    #   scientific_name = "Gadus morhua",
    #   res = 'medium',
    #   r.low = qlnorm(
    #     0.025,
    #     log(com_driors$growth_rate_prior),
    #     com_driors$growth_rate_prior_cv
    #   ),
    #   r.hi = qlnorm(
    #     0.975,
    #     log(com_driors$growth_rate_prior),
    #     com_driors$growth_rate_prior_cv
    #   ),
    #   start.yr = min(com_driors$years),
    #   stb.low = qlnorm(
    #     0.025,
    #     log(com_driors$initial_state),
    #     com_driors$initial_state_cv
    #   ),
    #   stb.hi = pmin(1,qlnorm(
    #     0.975,
    #     log(com_driors$initial_state),
    #     com_driors$initial_state_cv
    #   )),
    #   endb.low = qlnorm(
    #     0.025,
    #     log(com_driors$terminal_state),
    #     com_driors$terminal_state_cv
    #   ),
    #   endb.hi = pmin(1,qlnorm(
    #     0.975,
    #     log(com_driors$terminal_state),
    #     com_driors$terminal_state_cv
    #   )),
    #   intb.low = qlnorm(0.025, log(middle_state), middle_state_cv),
    #   intb.hi = pmin(1,qlnorm(0.975, log(middle_state), middle_state_cv)),
    #   int.yr = round(median(com_driors$years)),
    #   end.yr = max(com_driors$years),
    #   cores = 12,
    #
    # )

    posterior <- com_fit$fit %>%
      filter(variable == "dep_t") %>%
      mutate(prior = rlnorm(n(),
                            log(terminal_state),
                            terminal_state_cv)) %>%
      mutate(scaled_value = value  / max(value))

    # cmsy_fit %>%
    #   mutate(prior = rlnorm(
    #     n(),
    #     log(terminal_state),
    #     terminal_state_cv
    #   )) %>%
    #   ggplot(aes(prior)) +
    #   geom_density() +
    #   geom_vline(aes(xintercept = median_btv_lastyr))
    #

    return(posterior)

  }


catch_histories <- catch_histories %>%
  mutate(com_fits = map(catches, fit_example))


catch_history_plot <- catch_histories %>%
  select(form, catches) %>%
  unnest(cols = catches) %>%
  ggplot(aes(year, catch)) +
  geom_line() +
  facet_wrap( ~ form)

catches <- catch_histories %>%
  select(form, catches) %>%
  unnest(cols = catches)

catch_status_plot <- catch_histories %>%
  select(-catches) %>%
  unnest(cols = com_fits) %>%
  group_by(form, year) %>%
  mutate(value = mean(value)) %>%
  group_by(form) %>%
  mutate(value = value / max(value)) %>%
  left_join(catches, by  = c("form", "year")) %>%
  select(form, year, value, catch) %>%
  rename(`Mean Estimated B/K` = value, Catch = catch) %>%
  pivot_longer(c(`Mean Estimated B/K`, Catch),
               names_to = "type",
               values_to = "value") %>%
  group_by(form , type) %>%
  mutate(value = value / max(value)) %>%
  ggplot(aes(year, value, linetype = type)) +
  geom_line(alpha = 0.75) +
  facet_wrap( ~ form) +
  # scale_color_discrete(name = '') +
  scale_x_continuous(name = "Year", guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(name = "") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank()) +
  scale_linetype(name = '')



nothing_learned_plot <- catch_histories %>%
  unnest(cols = com_fits) %>%
  group_by(form) %>%
  filter(year == max(year)) %>%
  select(form, value, prior) %>%
  rename(posterior = value) %>%
  pivot_longer(-form, names_to = "type", values_to = "value") %>%
  ggplot(aes(value, fill = type)) +
  geom_density(alpha = 0.75) +
  facet_wrap( ~ form) +
  scale_x_continuous(name = "Final B/K", guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(name = "Density") +
  scale_fill_manual(
    labels = c("Posterior", "Prior"),
    name = '',
    values = c("lightgrey", "darkgrey")
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

catch_and_nothing_plot <- catch_status_plot + nothing_learned_plot +
  plot_annotation(tag_levels = 'A')


# train machine learning model to predict catches -------------------------

# classify stocks by stock history shape

ram_catches <- ram_data %>%
  mutate(catch = capture) %>%
  ungroup() %>%
  select(stockid, year, catch) %>%
  group_by(stockid) %>%
  mutate(stock_year = 1:length(catch),
         n_years = length(catch)) %>%
  mutate(scaled_catch = scale(catch)) %>%
  ungroup() %>%
  filter(n_years > 25,
         stock_year <= 25)

ram_catches %>%
  ggplot(aes(stock_year, scaled_catch, color = stockid)) +
  geom_line(show.legend = FALSE)

ram_catches <- ram_catches %>%
  select(stockid, stock_year, scaled_catch) %>%
  pivot_wider(names_from = stock_year, values_from = scaled_catch) %>%
  ungroup()

nstocks <- nrow(ram_catches)

map_dbl(ram_catches, ~ sum(is.na(.x)))

a = ram_catches %>% select(-stockid) %>% as.matrix()
set.seed(42)
catch_pca <- kernlab::specc(a, centers = 4)

# centers(catch_pca)
# size(catch_pca)
# withinss(catch_pca)

cluster <- as.vector(catch_pca)

ram_catches$cluster <- cluster


ram_catches <- ram_catches  %>%
  pivot_longer(c(-stockid,-cluster),
               names_to = "stock_year",
               values_to = "catch",) %>%
  mutate(stock_year = as.integer(stock_year))

ram_catches %>%
  ggplot(aes(stock_year, catch, group = stockid)) +
  geom_line() +
  facet_wrap( ~ cluster)

cluster_data <- ram_catches %>%
  pivot_wider(names_from = stock_year, values_from = catch) %>%
  ungroup() %>%
  mutate(cluster = as.factor(cluster)) %>%
  janitor::clean_names()

cluster_splits <-
  rsample::initial_split(cluster_data, strata = cluster)


cluster_model <-
  rand_forest(mtry = tune(),
              min_n = tune(),
              trees = 1000) %>%
  set_engine("ranger", num.threads = 8) %>%
  set_mode("classification")


cluster_recipe <-
  recipes::recipe(cluster ~ ., data = training(cluster_splits) %>% select(-stockid)) %>%
  themis::step_upsample(cluster)

cluster_workflow <-
  workflows::workflow() %>%
  workflows::add_model(cluster_model) %>%
  workflows::add_recipe(cluster_recipe)

val_set <- training(cluster_splits) %>% select(-stockid) %>%
  rsample::vfold_cv()

set.seed(345)
cluster_tuning <-
  cluster_workflow %>%
  tune_grid(
    val_set,
    grid = 20,
    control = control_grid(save_pred = TRUE),
    metrics = metric_set(roc_auc)
  )


best_forest <- cluster_tuning %>%
  select_best("roc_auc")

final_workflow <-
  cluster_workflow %>%
  finalize_workflow(best_forest)

cluster_fit <-
  final_workflow %>%
  fit(data = training(cluster_splits) %>% select(-stockid))

cluster_fit <- workflows::pull_workflow_fit(cluster_fit)

training_data <- training(cluster_splits) %>%
  bind_cols(predict(cluster_fit, new_data = .)) %>%
  mutate(split = "training")

testing_data <- testing(cluster_splits) %>%
  bind_cols(predict(cluster_fit, new_data = .)) %>%
  mutate(split = "testing")

cluster_predictions <- training_data %>%
  bind_rows(testing_data) %>%
  rename(predicted_cluster = .pred_class)

cluster_model_performance <- cluster_predictions %>%
  group_by(split, cluster) %>%
  summarise(accuracy = mean(cluster == predicted_cluster))

cluster_model_performance %>%
  ggplot(aes(cluster, accuracy, fill = split)) +
  geom_col(position = "dodge")

cluster_predictions %>%
  group_by(split) %>%
  summarise(accuracy = mean(cluster == predicted_cluster)) %>%
  pivot_wider(names_from = "split", values_from = "accuracy") %>%
  mutate(testing_loss = testing / training - 1)


status_model_data <- ram_data %>%
  mutate(catch = capture) %>%
  filter(stockid %in% unique(cluster_predictions$stockid)) %>%
  left_join(cluster_predictions %>% select(stockid, predicted_cluster),
            by = "stockid") %>%
  left_join(
    sraplus::fao_taxa$fao_species %>% select(scientific_name, isscaap_group),
    by = c("scientific_name" = "scientific_name")
  ) %>%
  group_by(stockid) %>%
  mutate(
    c_div_maxc = catch / max(catch, na.rm = TRUE),
    c_div_meanc = catch / mean(catch, na.rm = TRUE),
    fishery_year = 1:length(catch)
  ) %>%
  mutate(
    c_roll_meanc = RcppRoll::roll_meanr(c_div_meanc, 5),
    c_roll_maxc = catch / cummax(catch),
    c_init_slope = lm(log(catch[1:10] + 1e-3) ~ year[1:10])$coefficients[2]
  ) %>%
  gather(metric, value, b_v_bmsy, u_v_umsy, exploitation_rate) %>%
  select(stockid,
         year,
         contains('c_'),
         metric,
         value,
         predicted_cluster,
         fishery_year) %>%
  mutate(log_value = log(value + 1e-3)) %>%
  unique() %>%
  na.omit() %>%
  ungroup() %>%
  group_by(stockid) %>%
  filter(fishery_year > 20) %>%
  ungroup()

# OK now have dataframe ready to make predictions of stock status based on catch

# add in life history data

taxa <-
  tibble(scientific_name = unique(status_model_data$scientific_name))

get_lh <- function(tax) {
  # tax <- taxa$scientific_name[1]

  loc <- FishLife::Match_species(tax)

  lh <-
    t(as.matrix(FishLife::FishBase_and_RAM$beta_gv[loc$GroupNum, ])) %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    mutate(m_v_k = exp(m) / exp(k),
           linf_v_lmat = exp(loo) / exp(lm))
}

taxa <- taxa %>%
  mutate(lh = map(scientific_name, safely(get_lh)))

found_lh <- map_lgl(map(taxa$lh, "error"), is.null)

taxa <- taxa %>%
  filter(found_lh) %>%
  mutate(lh = map(lh, "result"))

taxa <- taxa %>%
  unnest(cols = lh)

b_data <- status_model_data %>%
  filter(metric == "b_v_bmsy",
         scientific_name %in% unique(taxa$scientific_name)) %>%
  left_join(ram_data %>% select(stockid, primary_country, fao_area_code) %>% unique(),
            by = "stockid") %>%
  left_join(taxa, by = "scientific_name")


b_data %>%
  group_by(stockid) %>%
  # filter(fishery_year == max(fishery_year)) %>%
  ungroup() %>%
  ggplot(aes(ln_fmsy, value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_continuous(name = "Recent B/Bmsy") +
  facet_wrap( ~ predicted_cluster)

# b_data %>%
#   group_by(stockid) %>%
#   filter(fishery_year == max(fishery_year)) %>%
#   ungroup() %>%
#   ggplot(aes(m_v_k, value)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   scale_y_continuous(name = "Recent B/Bmsy") +
#   facet_wrap(~predicted_cluster)

# generate splits
#

training_data <- b_data %>%
  filter(!primary_country %in% c("New Zealand", "Australia", "South Africa"))


testing_data <- b_data %>%
  filter(primary_country %in% c("New Zealand", "Australia", "South Africa"))

analysis_data <- training_data %>%
  rsample::group_vfold_cv(group = "fao_area_code")


analysis_data <- training_data %>%
  rsample::group_vfold_cv(group = "stockid", v = 5)

# prepare model workflow

# com_model <-
#   parsnip::rand_forest(
#     mode = "regression",
#     mtry = tune(),
#     min_n = tune(),
#     trees = 500
#   ) %>%
#   parsnip::set_engine("ranger")
#
# com_workflow <- workflow() %>%
#   add_model(com_model) %>%
#   add_formula(
#      value ~ c_div_maxc + c_div_meanc + c_length + c_roll_meanc + c_roll_maxc + c_init_slope + predicted_cluster + fishery_year + loo + linf_v_lmat + winfinity + tmax + tm + lm + temperature + ln_var + rho +
#       ln_masps + ln_margsd + h + logitbound_h + ln_fmsy_over_m + ln_fmsy + ln_r + ln_g + m_v_k
#   )

tune_grid <-
  parameters(
    min_n(range(2, 10)),
    tree_depth(range(2, 15)),
    learn_rate(range = c(-3,-.25)),
    mtry(),
    loss_reduction(),
    sample_prop(range = c(0.25, 1)),
    trees(range = c(10, 750))
  ) %>%
  dials::finalize(mtry(), x = training_data %>% select(-(1:2)))

xgboost_grid <- grid_latin_hypercube(tune_grid, size = 40)
xgboost_grid$learn_rate %>% hist()
xgboost_model <-
  parsnip::boost_tree(
    mode = "regression",
    mtry = tune(),
    min_n = tune(),
    loss_reduction = tune(),
    sample_size = tune(),
    learn_rate = tune(),
    tree_depth = tune(),
    trees = tune()
  ) %>%
  parsnip::set_engine("xgboost")

xgboost_workflow <- workflow() %>%
  add_model(xgboost_model) %>%
  add_formula(
    value ~ c_div_maxc + c_div_meanc + c_roll_meanc + c_roll_maxc + c_init_slope + predicted_cluster + fishery_year + loo + linf_v_lmat + winfinity  + ln_var + h + ln_fmsy_over_m + m_v_k
  )

# xgboost_workflow <- workflow() %>%
#   add_model(xgboost_model) %>%
#   add_formula(
#     value ~ c_div_maxc + c_div_meanc + c_roll_meanc + c_roll_maxc + c_init_slope + predicted_cluster + fishery_year)

if (tune_com) {
  workers <- 4 # for some reason parallel computation is crashing

  cl <- makePSOCKcluster(workers)
  #
  registerDoParallel(cl)
  #
  getDoParName()
  #
  getDoParWorkers()
  # #
  a <- Sys.time()

  xgboost_tuning <- tune_grid(
    xgboost_workflow,
    resamples = analysis_data,
    grid = xgboost_grid,
    control = control_grid(save_pred = FALSE)
  )

  Sys.time() - a


  write_rds(xgboost_tuning, file = file.path(results_path,"com_tunegrid.rds"))

  # write_rds(com_tunegrid, file = "com_tunegrid.rds")

} else {

  xgboost_tuning <- readr::read_rds(file = file.path(results_path,"com_tunegrid.rds"))

}
autoplot(xgboost_tuning, metric = "rmse") +
  scale_y_continuous(limits = c(NA, 1))

best_vals <- tune::select_best(xgboost_tuning, metric = "rmse")

# best_vals$trees <- 300

com_workflow <- finalize_workflow(xgboost_workflow,
                                  best_vals)

# finalize model

com_model <- com_workflow %>%
  fit(data = training_data)

vip(com_model$fit$fit$fit)

# fit to LOO sample


if (run_loo) {
  loo_test <- training_data %>%
    rsample::group_vfold_cv(group = "stockid") %>%
    mutate(fit = map(splits,  ~ fit(com_workflow, analysis(.x))))

  loo_test <- loo_test %>%
    mutate(pred = map2(fit, splits, ~ cbind(
      assessment(.y), predict(.x, new_data = assessment(.y))
    )))

  write_rds(loo_test %>% select(-fit),
            file = file.path(results_path, "loo_test.rds"))

} else {
  loo_test <- read_rds(file.path(results_path, "loo_test.rds"))

}

loo_results <- loo_test %>%
  select(pred) %>%
  unnest(cols = pred) %>%
  mutate(split = "Leave-One-Out")

# write_rds(com_model, file = "com_model.rds")

# com_vip <- com_model %>%
#   pull_workflow_fit() %>%
#   vip()


com_eval_data <- training_data %>%
  cbind(predict(com_model, new_data =  training_data)) %>%
  mutate(split = "Training") %>%
  rbind(testing_data %>%
          cbind(predict(com_model, new_data =  testing_data)) %>%
          mutate(split = "Testing")) %>%
  rbind(loo_results) %>%
  mutate(split = fct_relevel(split, c("Training", "Leave-One-Out", "Testing")))

com_rsq <- com_eval_data %>%
  group_by(split) %>%
  yardstick::rsq(truth = value, estimate = .pred)

com_rmse <- com_eval_data %>%
  group_by(split) %>%
  yardstick::rmse(truth = value, estimate = .pred)

com_mae <- com_eval_data %>%
  group_by(split) %>%
  yardstick::mae(truth = value, estimate = .pred)

com_eval_data %>%
  group_by(split) %>%
  yardstick::mae(truth = value, estimate = .pred)



com_rf_plot <- com_eval_data %>%
  ungroup() %>%
  group_by(stockid) %>%
  filter(fishery_year > (max(fishery_year) - 10)) %>%
  ggplot(aes(value, .pred)) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  geom_bin2d(alpha = 0.8, aes(fill = after_stat(density))) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  ggtext::geom_richtext(data = com_rsq, aes(
    x = 1,
    y = 3,
    label = paste0("R<sup>2</sup> = ", round(.estimate, 2))
  )) +
  ggtext::geom_richtext(data = com_rmse, aes(
    x = 1,
    y = 2.75,
    label = paste0("RMSE = ", round(.estimate, 2))
  )) +
  facet_wrap( ~ split) +
  scale_x_continuous(name = "RLSADB B/Bmsy", limits = c(0, NA)) +
  scale_y_continuous(name = "Predicted B/Bmsy", limits = c(0, NA)) +
  scale_fill_viridis_c(
    name = "Density",
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      barheight = unit(27, "lines")
    )
  )

com_rf_plot

com_eval_data %>%
  group_by(split) %>%
  yardstick::rsq(truth = value, estimate = .pred)


com_eval_data %>%
  group_by(split) %>%
  yardstick::rmse(truth = value, estimate = .pred)

com_eval_data %>%
  group_by(split) %>%
  yardstick::mae(truth = value, estimate = .pred)


# fit to pre 1990

older_com_model <- com_workflow %>%
  fit(data = training_data %>% filter(year < 1990))


# fit to LOO sample


if (run_loo) {
  older_loo_test <- training_data %>%
    filter(year < 1990) %>%
    rsample::group_vfold_cv(group = "stockid") %>%
    mutate(fit = map(splits,  ~ fit(com_workflow, analysis(.x))))

  older_loo_test <- older_loo_test %>%
    mutate(pred = map2(fit, splits, ~ cbind(
      assessment(.y), predict(.x, new_data = assessment(.y))
    )))

  write_rds(older_loo_test %>% select(-fit),
            file = file.path(results_path, "older_loo_test.rds"))

} else {
  older_loo_test <-
    read_rds(file = file.path(results_path, "older_loo_test.rds"))


}

older_loo_results <- older_loo_test %>%
  select(pred) %>%
  unnest(cols = pred) %>%
  mutate(split = "Leave-One-Out")

# write_rds(com_model, file = "com_model.rds")

# com_vip <- com_model %>%
#   pull_workflow_fit() %>%
#   vip()


older_com_eval_data <- training_data %>%
  filter(year < 1990) %>%
  cbind(predict(older_com_model, new_data =  training_data %>% filter(year < 1990))) %>%
  mutate(split = "Training") %>%
  rbind(testing_data %>%
          filter(year < 1990) %>%
          cbind(predict(
            older_com_model, new_data =  testing_data %>% filter(year < 1990)
          )) %>%
          mutate(split = "Testing")) %>%
  rbind(older_loo_results) %>%
  mutate(split = fct_relevel(split, c("Training", "Leave-One-Out", "Testing")))


older_com_rsq <- older_com_eval_data %>%
  group_by(split) %>%
  yardstick::rsq(truth = value, estimate = .pred)

older_com_rmse <- older_com_eval_data %>%
  group_by(split) %>%
  yardstick::rmse(truth = value, estimate = .pred)


older_com_rf_plot <- older_com_eval_data %>%
  ungroup() %>%
  group_by(stockid) %>%
  ggplot(aes(value, .pred)) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2) +
  geom_bin2d(alpha = 0.8, aes(fill = after_stat(density))) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  ggtext::geom_richtext(data = older_com_rsq, aes(
    x = 1,
    y = 3,
    label = paste0("R<sup>2</sup> = ", round(.estimate, 2))
  )) +
  ggtext::geom_richtext(data = older_com_rmse, aes(
    x = 1,
    y = 2.75,
    label = paste0("RMSE = ", round(.estimate, 2))
  )) +
  facet_wrap( ~ split) +
  scale_x_continuous(name = "RLSADB B/Bmsy", limits = c(0, NA)) +
  scale_y_continuous(name = "Predicted B/Bmsy", limits = c(0, NA)) +
  scale_fill_viridis_c(
    name = "Density",
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      barheight = unit(27, "lines")
    )
  )
older_com_rf_plot

older_com_eval_data %>%
  group_by(split) %>%
  yardstick::rsq(truth = value, estimate = .pred)


older_com_eval_data %>%
  group_by(split) %>%
  yardstick::rmse(truth = value, estimate = .pred)



# save things -------------------------------------------------------------

plots <- ls()[str_detect(ls(), "_plot")]

save(file = file.path(results_path,"eval-cmsy-plots.Rdata"), list = plots)
