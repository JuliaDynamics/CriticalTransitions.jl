#! format: off
pages = [
    "Home" => "index.md",
    "Tutorial" => "examples/tutorial.md",
    "Defining a system" => Any[
        "Stochastic system" => "examples/stochastic-dynamics.md",
        "Nonautonomous system" => "examples/RateSystem.md",],
    "Examples" => Any[
        "Large deviations: Maier-Stein system" => "examples/gMAM_Maierstein.md",
        "Quasipotential: Maier-Stein system" => "examples/quasipotential_maierstein.md",
        "Simple gMAM: Kerr Parametric Oscillator" => "examples/sgMAM_KPO.md",
        "Adaptive step-size control for sgMAM" => "examples/backtracking_KPO.md",
        "Transition path theory: Finite element method" => "examples/transition_path_theory_double_well.md",
        "Multiple shooting method" => "examples/shooting_Maierstein.md",
        "String method: Muller-Brown potential" => "examples/potential_string.md",
    ],
    "Manual" => Any[
        "Choosing an approach" => "man/overview.md",
        "Setting up" => Any[
            "Define your system" => "man/system_construction.md",
            "Stability & basins" => "man/systemanalysis.md",
            "Simulating forward" => "man/simulation.md",
        ],
        "Studying transitions" => Any[
            "Direct sampling" => "man/sampling.md",
            "Large deviation theory" => "man/largedeviations.md",
            "Rates, distributions & the generator" => "man/diffusion_operator.md",
            "Transition path theory" => "man/transition_path_theory.md",
            "Rate-induced transitions" => "man/r-tipping.md",
        ],
        "Utilities" => "man/utils.md",
        "Developer / internals" => "man/dev.md",
        "Bibliography" => "man/bibliography.md",
    ],
    # "Predefined systems" => "man/systems.md",
    "API" => "api.md",
    "References" => "refs.md",
]
#! format: on
