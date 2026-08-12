#! format: off
pages = [
    "Home" => "index.md",
    "Tutorial" => "examples/tutorial.md",
    "Manual" => Any[
        "Overview" => "man/overview.md",
        "Constructing a system" => "man/system_construction.md",
        "Simulation" => "man/simulation.md",
        "Stability analysis" => "man/systemanalysis.md",
        "Stochastic dynamics" => Any[
            "Direct sampling" => "man/sampling.md",
            "Large deviation theory" => "man/largedeviations.md",
            "Rates, distributions & the generator" => "man/diffusion_operator.md",
            "Transition path theory" => "man/transition_path_theory.md",
        ],
        "Nonautonomous dynamics" => Any[
            "Rate-induced tipping" => "man/r-tipping.md",
        ],
        "Utilities" => "man/utils.md",
        "Developer / internals" => "man/dev.md",
    ],
    "Examples" => Any[
        "Large deviations" => Any[
            "Large deviations: Maier-Stein system" => "examples/gMAM_Maierstein.md",
            "Quasipotential: Maier-Stein system" => "examples/quasipotential_maierstein.md",
            "Simple gMAM: Kerr Parametric Oscillator" => "examples/sgMAM_KPO.md",
            "Adaptive step-size control for sgMAM" => "examples/backtracking_KPO.md",
        ],
        "More" => Any[
            "Transition path theory: Finite element method" => "examples/transition_path_theory_double_well.md",
            "Multiple shooting method" => "examples/shooting_Maierstein.md",
            "String method: Muller-Brown potential" => "examples/potential_string.md",
        ],
    ],
    "References" => "refs.md",
]
#! format: on
