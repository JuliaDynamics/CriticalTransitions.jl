window.BENCHMARK_DATA = {
  "lastUpdate": 1758795939211,
  "repoUrl": "https://github.com/JuliaDynamics/CriticalTransitions.jl",
  "entries": {
    "Benchmark Results": [
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "80ecdfa7e06758b42fdfb45aafc27cc04f96459f",
          "message": "docs: overhaul (#178)\n\nFull revision of docs, first logo draft\n\n* fix install command\n\n* updated CoupledSDEs docs\n\n* added logo\n\n* Added flowchart.drawio\n\n* worked on docs\n\n* fixed CoupledSDEs docstring not showing\n\n* revised quickstart page\n\n* correct \"conecpts\" to \"concepts\"\n\n* worked on docs\n\n* format\n\n* added references page\n\n* revised Manual sections and index page\n\n* format\n\n* fix merge\n\n---------\n\nCo-authored-by: reykboerner <reyk.boerner@reading.ac.uk>\nCo-authored-by: reykboerner <r.borner@uu.nl>",
          "timestamp": "2025-08-06T12:21:56+02:00",
          "tree_id": "49a3c4da2bfcd92669a4e8229083ad7391ccdff0",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/80ecdfa7e06758b42fdfb45aafc27cc04f96459f"
        },
        "date": 1754476057422,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 90124715,
            "unit": "ns",
            "extra": "gctime=11905503\nmemory=220494096\nallocs=410240\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2333895249,
            "unit": "ns",
            "extra": "gctime=375207135\nmemory=4138029184\nallocs=112014889\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 714891421,
            "unit": "ns",
            "extra": "gctime=112004064.5\nmemory=2207174816\nallocs=3377834\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "41898282+github-actions[bot]@users.noreply.github.com",
            "name": "github-actions[bot]",
            "username": "github-actions[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "8a431e9849a4524426774cbd5d50eb578ac58154",
          "message": "CompatHelper: bump compat for IntervalArithmetic to 1, (keep existing compat) (#200)\n\n* CompatHelper: bump compat for IntervalArithmetic to 1, (keep existing compat)\n\n* update reference\n\n* fix typo, bump verison\n\n---------\n\nCo-authored-by: CompatHelper Julia <compathelper_noreply@julialang.org>\nCo-authored-by: reykboerner <r.borner@uu.nl>",
          "timestamp": "2025-09-25T12:20:08+02:00",
          "tree_id": "b33b8b705cd8592e95034398f68b67d2205977c4",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/8a431e9849a4524426774cbd5d50eb578ac58154"
        },
        "date": 1758795936988,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 92238073.5,
            "unit": "ns",
            "extra": "gctime=13385051.5\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2274700691,
            "unit": "ns",
            "extra": "gctime=390285010\nmemory=4138027152\nallocs=112014837\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 745922108.5,
            "unit": "ns",
            "extra": "gctime=116724530\nmemory=2207190160\nallocs=3378334\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}