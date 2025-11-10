#! format: off
pages = [
    "Home" => "index.md",
    "Getting started" => "quickstart.md",
    "Examples" => Any[
        "Tutorial" => "examples/tutorial.md",
        "Defining a system" => Any[
            "Stochastic system" => "examples/stochastic-dynamics.md",
            "Nonautonomous system" => "examples/RateSystem.md",],
        "System analysis" => Any[  
            "Large deviations: Maier-Stein system" => "examples/gMAM_Maierstein.md",
            "Simple gMAM: Kerr Parametric Oscillator" => "examples/sgMAM_KPO.md",
            "Minimal action method: Optimal Control problem" => "examples/OC_mam.md",
            "Transition path theory: Finite element method" => "examples/transition_path_theory_double_well.md",],
    ],
    "Manual" => Any[
        "Define your system" => "man/system_construction.md",
        "Stability analysis" => "man/systemanalysis.md",
        "Simulating the system" => "man/simulation.md",
        "Sampling transitions" => "man/sampling.md",
        "Large deviation theory" => "man/largedeviations.md",
        "Rate-induced transitions" => "man/r-tipping.md",
        "Transition path theory" => "man/transition_path_theory.md",
        "Utilities" => "man/utils.md",
        "Bibliography" => "man/bibliography.md",
    ],
    # "Predefined systems" => "man/systems.md",
    # "Development stage" => "man/dev.md",
    "References" => "refs.md"
]
#! format: on
