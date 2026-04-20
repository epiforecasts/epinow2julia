# r_gt_opts_to_julia round-trips a fixed lognormal

    Code
      cat(show_julia(julia_opts))
    Output
      Generation time options:
        - lognormal distribution:
          meanlog:
            1.6
          sdlog:
            0.5

# r_delay_opts_to_julia round-trips a fixed lognormal

    Code
      cat(show_julia(julia_opts))
    Output
      Delay options:
        - lognormal distribution:
          meanlog:
            0.6
          sdlog:
            0.5

# r_rt_opts_to_julia round-trips a weekly random walk

    Code
      cat(show_julia(julia_opts))
    Output
      Rt options:
        prior: lognormal(mean=2.0, sd=0.2)
        use_rt: true
        random walk period: 7
        gp_on: gp_Rt
        future: latest

# r_obs_opts_to_julia round-trips negbin + week effect

    Code
      cat(show_julia(julia_opts))
    Output
      Observation model options:
        family: negbin
        week_effect: true

# r_forecast_opts_to_julia round-trips a 7-day horizon

    Code
      cat(show_julia(julia_opts))
    Output
      EpiNow2.ForecastOpts(7)

# r_inference_opts_to_julia maps sampling args

    Code
      cat(show_julia(julia_opts))
    Output
      Inference options:
        sampler: nuts
        samples: 250, warmup: 100, chains: 2
        adtype: AutoReverseDiff(compile=true)

