make_nothing_learned_plot <- function(plot_col, catch_histories) {
  catch_history_plot <- catch_histories %>%
    select(form, catches) %>%
    unnest(cols = catches) %>%
    ggplot(aes(year, catch)) +
    geom_line() +
    facet_wrap(~ form)

  catches <- catch_histories %>%
    select(form, catches) %>%
    unnest(cols = catches)

  catch_status_plot <- catch_histories %>%
    select(-catches) %>%
    unnest(cols = all_of(plot_col)) %>%
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
    facet_wrap(~ form) +
    # scale_color_discrete(name = '') +
    scale_x_continuous(name = "Year", guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(name = "") +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      strip.background = element_blank()
    ) +
    scale_linetype(name = '')



  nothing_learned_plot <- catch_histories %>%
    unnest(cols = all_of(plot_col)) %>%
    group_by(form) %>%
    filter(year == max(year)) %>%
    select(form, value, prior) %>%
    rename(posterior = value) %>%
    pivot_longer(-form, names_to = "type", values_to = "value") %>%
    ggplot(aes(pmin(value, 1), fill = type)) +
    geom_density(alpha = 0.75) +
    facet_wrap(~ form) +
    scale_x_continuous(name = "Final B/K", guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(name = "Density") +
    scale_fill_manual(
      labels = c("Posterior", "Prior"),
      name = '',
      values = c("lightgrey", "darkgrey")
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "top",
          strip.background = element_blank())

  catch_and_nothing_plot <-
    catch_status_plot + nothing_learned_plot +
    plot_annotation(tag_levels = 'A')

  ggsave(
    filename = file.path(results_path, paste0("fig_2_", plot_col, ".pdf")),
    catch_and_nothing_plot,
    dev = cairo_pdf,
    width = 180,
    height = 100,
    units = "mm",
    dpi = 600
  )

  return(catch_and_nothing_plot)

}
