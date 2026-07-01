window.BENCHMARK_DATA = {
  "lastUpdate": 1782283788176,
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
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=375207135\nmemory=4138029184\nallocs=112014889\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 714891421,
            "unit": "ns",
            "extra": "gctime=112004064.5\nmemory=2207174816\nallocs=3377834\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
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
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=390285010\nmemory=4138027152\nallocs=112014837\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 745922108.5,
            "unit": "ns",
            "extra": "gctime=116724530\nmemory=2207190160\nallocs=3378334\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "3f3d126feeb7f55ade7a533174c7e0eceed10c14",
          "message": "build: benchmark on v1.11 (#206)",
          "timestamp": "2025-10-19T13:41:20+02:00",
          "tree_id": "56466bc78cf5dcdc92831a2483e55401bbcfe8a6",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/3f3d126feeb7f55ade7a533174c7e0eceed10c14"
        },
        "date": 1760874412726,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 88335747,
            "unit": "ns",
            "extra": "gctime=11842418\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=328491072\nmemory=4138027152\nallocs=112014837\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 672580551,
            "unit": "ns",
            "extra": "gctime=94337914\nmemory=2207190624\nallocs=3378374\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "88832655+raphael-roemer@users.noreply.github.com",
            "name": "raphael-roemer",
            "username": "raphael-roemer"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c43082599e29174c5a22540a89bf54181baac445",
          "message": "feat: add RateSystem (#117)\n\n* added first draft for RateSystem\n\n* incorporated RateSystem file\n\n* small typo fix\n\n* undid wrong fix\n\n* Removed RateSystem Lines\n\n* Return RateSystem changes\n\n* Again delete Rate stuff\n\n* Include Rate stuff again\n\n* added new RateSystem Version\n\n* Compressed RateSystemDraft2.jl script, created RateSystemTestsDraft.jl script, created a dev folder, updated src/CriticalTransitions.jl file.\n\n* Corrected src/CriticalTransitions.jl file.\n\n* Moved contents of dev/ file to test/ file and reset src/CriticalTransitions.jl file to be in line with that of the main branch.\n\n* Collected all old and new RateSystem scripts into the test/ratesystem folder.\n\n* yep\n\n* Added PolynomialRoots to dependencies.\n\n* See previous.\n\n* See previous again.\n\n* Exported truscottbrindley-mod-gen-det function.\n\n* Removed CoupledODEs functionality from systems/ folder.\n\n* Removed the system I added.\n\n* Starting a RateSystem Example\n\n* work on example of RateSystem\n\n* work on example of RateSystem\n\n* work on example of RateSystem\n\n* work on example of RateSystem\n\n* added documentation in RateSystem file\n\n* corrected mistake in RateSystem file\n\n* improved documentation of test2\n\n* work on the example for RateSystem\n\n* added codes for the RateSystem example\n\n* correction in MaierStein example\n\n* small corrections in Test and RateSystemDraft\n\n* addition to RateSystem example\n\n* additions to RateSystem example\n\n* additions to RateSystem example\n\n* correction in RateSystem example\n\n* correction in RateSystem.jl\n\n* added RateSystem.md to pages.jl\n\n* moved RateSystem source code to src\n\n* applied formatter, fixed typo, disabled spell check until PR review\n\n* remove CairoMakie dep\n\n* added RateSystem test\n\n* small fix\n\n* applied Formatter\n\n* small edits in docs and added docstring drafts\n\n* deleted test/ratesystem/RateSystem.jl\n\n* deleted test/ratesystem/RateSystemDraft1.jl and test/ratesystem/RateSystemTestDraft1.jl\n\n* deleted test/ratesystem\n\n* deleted examples/RateSystem.jl\n\n* expanded documentaation of RateSystem\n\n* expanded documentation\n\n* added plot of parameter shift in RateSystem documentation\n\n* deleted NLPModelsIpopt from project.toml\n\n* fixed error in Project.toml\n\n* fixed quotation marks in quickstart.md\n\n* removed the ! from apply_ramping\n\n* fixed error in documentation of RateSystem\n\n* resolving documentation issue with RateProtocol-plot\n\n* improved documentation of apply_ramping\n\n* correction in tes of RaateSystem\n\n* corrected typo in RateProtocol docstring\n\n* added plot of lambda to example of RateSystem\n\n* removed unnecessary formatting information and improved beginning of RateSystem example\n\n* correction in docstring\n\n* recovered docs/Project.toml\n\n* updated .toml and applied Formatter\n\n* fixed typo in RateSystem.md\n\n* enable spell checking again\n\n* imporved docstring of apply_ramping\n\n* improved docstring of RateProtocol\n\n* imporved docstring of apply_ramping\n\n* imporved docstring of apply_ramping\n\n* imporved docstring of apply_ramping\n\n* fixed spelling mistake in RateSystem example\n\n* improved RateSytem Example\n\n* changed name of R-Tipping documentation page\n\n* improved RateSytem Example\n\n* fix spelling\n\n* changed t_start,t_end to time_interval\n\n* corrected typo\n\n* corrected typo\n\n* changed how to access the dynamic rule within modified_drift\n\n* changed how to access the dynamic rule within apply_ramping\n\n* included referrenced_sciml_prob in exported functions from DynamicalSystemsBase\n\n* changed back to t_start and t_end to allow compatibility with the R-Tipping functionality branch, but stored the time_interval version in new file r_tipping/RateSystemTimeInterval.jl\n\n* adapted docs as well\n\n* adapted test as well\n\n* correction in test\n\n* changed RateProtocol to RateConfig\n\n* changed RateProtocol to RateConfig in .md file\n\n* changed p_lambda to p_parameter\n\n* corrected typo\n\n* corrected typo\n\n* corrected typo\n\n* changed lambda to p throughout\n\n* corrected typo\n\n* corrected typo\n\n* 1st try implementing dt dp setting\n\n* changed error in test\n\n* exported RateConfig\n\n* added stretching/compression of ramping given by ramp_t_length\n\n* change of reference window of ramping function\n\n* change of reference window of ramping function\n\n* changed rp to rc\n\n* adapted example in docs to explain dp dt structure\n\n* further adapted example in docs to explain dp dt structure\n\n* corrected typo\n\n* improved example of r-tipping\n\n* changed the start and end point conditions to recommendations\n\n* updated docstrings\n\n* Implemented changes to the RateSystem.jl source file (from discussions on 25.08.25).\n\n* 1st and 2ndIteration of RateSystem\n\n* Added the third iteration of RateSystem.jl.\n\n* typo corrected\n\n* typo corrected\n\n* typo corrected\n\n* typo corrected\n\n* typo corrected\n\n* typo corrected\n\n* typo corrected\n\n* added example CriticalTransitions.jl/examples/RateSys3rdIterTest.jl\n\n* Updated example test RateSys3rdIterTest.jl and made corrections to RateSystem.jl.\n\n* corrected typos\n\n* corrected typos\n\n* added benchmarking in 'examples/RateSys3rdIterTest.jl' to test the speed of the RateSystem setup\n\n* corrected error\n\n* made a deepcop of auto_sys to avoid changing auto_sys\n\n* added docs\n\n* added documentation of the example\n\n* adapted docstrings, updated tests\n\n* corrected docstrings\n\n* corrected docstrings\n\n* corrected docstrings\n\n* improved docs\n\n* started changes towards RateSystem interface\n\n* small addition\n\n* corrected typo\n\n* corrected typo\n\n* corrected typo\n\n* added example for new RateSystem version\n\n* corrected errors in RateSystem.jl\n\n* corrected error in RateSystem.jl\n\n* corrected error in RateSystem.jl\n\n* Added documentation in RateSys4thIter.jl\n\n* Renamed RateSys4thIter.jl to RateSys4thIterTest.jl\n\n* added and expanded docstrings in RateSystem.jl\n\n* corrected docstrings in RateSystem.jl\n\n* corrected docstrings in RateSystem.jl\n\n* corrected docstrings in RateSystem.jl\n\n* corrected docstrings in RateSystem.jl\n\n* corrected docstrings in RateSystem.jl\n\n* Added test for new RateSystem setup\n\n* added .md file for current version\n\n* added linkcheck_ignore = [..] to make.jl to ignore this broken link from EXTERNAL docstrings - internal linkcheck is not affected\n\n* trying to fix docs error\n\n* trying to fix docs error\n\n* trying to fix docs error\n\n* trying to fix docs error\n\n* fixed docs\n\n* format and fix spell check\n\n* fix spell check\n\n* add :referrenced_sciml_prob to public import exception\n\n* update benchmarks to use benchamrtools\n\n* removed placeholder reference\n\n* replace 'section_start/end' fields with 'interval' tuple\n\n* work on making RateSystem a proper ContinuousTimeDynamicalSystem subtype\n\n* removed outdated RateSystem files\n\n* format\n\n* remove copied CoupledODEs source code\n\n* build(deps): bump JET to v0.10\n\n* feat: ContinuousTimeDynamicalSystem helper functions\n\n* format\n\n* add obtaining info tests\n\n* make `RateConfig` concrete\n\n* make RateSystemForcing concrete\n\n* fix example\n\n* docs: add DynamicalSystemsBase and Attractors as remotes\n\n* add trajectory method\n\n* fix example to use updated RateConfig reference for parameter function\n\n* format\n\n* simplify structurce and fix performance\n\n* fix spelling\n\n* fix docs\n\n* docs: update example title\n\n* fix docs\n\n* refactor after discussion with @Datseris and @oameye\n\n* remove unused import\n\n* fix get_forcing and add deepcopy's for debugging\n\n* added debug test script\n\n* fix it: use the agreed scenario of NDR becoming the system equations of motion\n\n* deepcopy the params for rate system (debatable)\n\n* clean tests, rename forcing_length to forcing_duration\n\n* finish docstrings\n\n* update docs\n\n* fix benchmarks\n\n* fix example\n\n* format\n\n* small changes in docs\n\n---------\n\nCo-authored-by: Orjan Ameye <orjan.ameye@hotmail.com>\nCo-authored-by: Reyk Börner <reyk.boerner@reading.ac.uk>\nCo-authored-by: Ryan Deeley <ryan.deeley@uni-oldenburg.de>\nCo-authored-by: reykboerner <r.borner@uu.nl>\nCo-authored-by: Datseris <datseris.george@gmail.com>",
          "timestamp": "2025-11-03T21:39:46+01:00",
          "tree_id": "f5b48105fd377406c450336fcc6eca07be3ea1a8",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/c43082599e29174c5a22540a89bf54181baac445"
        },
        "date": 1762202782106,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 90216835,
            "unit": "ns",
            "extra": "gctime=12344444\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=377680889\nmemory=4138027152\nallocs=112014837\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 741955261.5,
            "unit": "ns",
            "extra": "gctime=121679369.5\nmemory=2207520464\nallocs=3397504\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 66855,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 69650,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "r.borner@uu.nl",
            "name": "Reyk Börner",
            "username": "reykboerner"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "971e6889df568378d09bc4e958fb3c7306e1c3a1",
          "message": "bugfix: frozen_system and docs update (#210)\n\n* update install instructions in docs\n\n* fix frozen_system bug\n\n* bump MTK compat in test/Project.toml\n\n* format\n\n* replace mtkcompile with structural_simplify again\n\n* streamline RateSystem example\n\n* refine docs\n\n---------\n\nCo-authored-by: Orjan Ameye <orjan.ameye@hotmail.com>",
          "timestamp": "2025-11-11T00:42:26+01:00",
          "tree_id": "bc8811545edf9bbd6dd1cbf3aae59bf44ff57e56",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/971e6889df568378d09bc4e958fb3c7306e1c3a1"
        },
        "date": 1762818544001,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 88126716,
            "unit": "ns",
            "extra": "gctime=11824751\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=365835574\nmemory=4138027152\nallocs=112014837\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 736896052,
            "unit": "ns",
            "extra": "gctime=118591340.5\nmemory=2207479728\nallocs=3397954\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 60053,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40406,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "71acaed2aa1011a51709e86f6b24fe4a48600f32",
          "message": "build(deps): bump actions/checkout from 5 to 6 (#215)\n\nBumps [actions/checkout](https://github.com/actions/checkout) from 5 to 6.\n- [Release notes](https://github.com/actions/checkout/releases)\n- [Changelog](https://github.com/actions/checkout/blob/main/CHANGELOG.md)\n- [Commits](https://github.com/actions/checkout/compare/v5...v6)\n\n---\nupdated-dependencies:\n- dependency-name: actions/checkout\n  dependency-version: '6'\n  dependency-type: direct:production\n  update-type: version-update:semver-major\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2025-11-24T17:28:22+01:00",
          "tree_id": "e024e3edd9f7cc12bf5975e115ee0e6061de36d2",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/71acaed2aa1011a51709e86f6b24fe4a48600f32"
        },
        "date": 1764002101225,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 112947706,
            "unit": "ns",
            "extra": "gctime=20364923\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=507025414\nmemory=4138027152\nallocs=112014837\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 770254628,
            "unit": "ns",
            "extra": "gctime=137569934\nmemory=2207148416\nallocs=3377234\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 78566,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 54700,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "ffed1b77378c1d65ab186234dcc3d53358022f68",
          "message": "fix(deps): set LinearSolve compat (#221)\n\n* fix(deps): set LinearSolve compat\n\n* Nonlinear solve compat",
          "timestamp": "2026-01-17T11:30:33+01:00",
          "tree_id": "c9a0c55c4f2705f23232cba78b7781b83fd9ee3c",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/ffed1b77378c1d65ab186234dcc3d53358022f68"
        },
        "date": 1768646235355,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 89037247.5,
            "unit": "ns",
            "extra": "gctime=11938786.5\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=338556306\nmemory=4138027152\nallocs=112014837\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 681720791,
            "unit": "ns",
            "extra": "gctime=97396866\nmemory=2207206672\nallocs=3377684\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 65874,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40616,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "59bfdaf3b829a561efeab3583c5c72a98ac5cc0e",
          "message": "fix: add parameter handling in H_x and H_p functions (#218)\n\n* fix: add parameter handling in H_x and H_p functions\n\n* format",
          "timestamp": "2026-01-17T11:59:59+01:00",
          "tree_id": "173248f66ae33ae1ef998d91b398301867553711",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/59bfdaf3b829a561efeab3583c5c72a98ac5cc0e"
        },
        "date": 1768647811880,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 89734682,
            "unit": "ns",
            "extra": "gctime=12156706\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=385821597\nmemory=4138027152\nallocs=112014837\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 739759275.5,
            "unit": "ns",
            "extra": "gctime=119635529.5\nmemory=2207144896\nallocs=3378104\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 83847,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 61224,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "6436d5f981a9cd4285874ee128c2782a293b523a",
          "message": "build(deps): bump JET (#226)\n\n* build(deps): bump JET\n\n* fix JET tests",
          "timestamp": "2026-01-17T18:49:57+01:00",
          "tree_id": "a4264a6357fc3a617e8430feb919aabeefde06e8",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/6436d5f981a9cd4285874ee128c2782a293b523a"
        },
        "date": 1768672407525,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 88641085,
            "unit": "ns",
            "extra": "gctime=11677086\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=320934387\nmemory=4138027152\nallocs=112014837\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 664857005,
            "unit": "ns",
            "extra": "gctime=91212799\nmemory=2207242032\nallocs=3378694\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 58749,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40315,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "198982749+Copilot@users.noreply.github.com",
            "name": "Copilot",
            "username": "Copilot"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "423ae9fd0332381319e7669c766e17a7c9b1be0a",
          "message": "Add warning when Nmax reached before generating N transitions (#223)\n\n* Initial plan\n\n* Add warning when Nmax is reached in transitions function\n\nCo-authored-by: oameye <57623933+oameye@users.noreply.github.com>\n\n* fix: tmax kwargs is respected in transitions\n\n* commit to check the new stats\n\n* add more comparison values\n\n* make tests less strict\n\n* Refactor tests for fitzhugh_nagumo: remove debug output and ensure Nmax warning is properly tested\n\n* Add warning for maximum attempts in fitzhugh_nagumo transitions test\n\n* fix warning test\n\n* restict warning test to release\n\n* fix version string\n\n* Fix condition for Nmax warning test in fitzhugh_nagumo transitions\n\n* fix it once and for all?\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: oameye <57623933+oameye@users.noreply.github.com>\nCo-authored-by: Orjan Ameye <orjan.ameye@hotmail.com>",
          "timestamp": "2026-01-18T09:46:52+01:00",
          "tree_id": "956e9de6def742efe9a2a00768fdb558c76b93ff",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/423ae9fd0332381319e7669c766e17a7c9b1be0a"
        },
        "date": 1768726232536,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 89452577.5,
            "unit": "ns",
            "extra": "gctime=12074430.5\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=347660829\nmemory=4138027152\nallocs=112014837\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 714846631,
            "unit": "ns",
            "extra": "gctime=108644146\nmemory=2207154384\nallocs=3378344\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 58359,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40196,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "255d32ce5761ff0c33d0bb5cce7f3e0d55fd2946",
          "message": "build(deps): update (Non)LinearSolve (#225)\n\n* build(deps): update (Non)LinearSolve\n\n* commit to get new stats",
          "timestamp": "2026-01-18T10:00:54+01:00",
          "tree_id": "9d5ad0fa3f5b8a0bf4cc6a083bfacd0dd60f922e",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/255d32ce5761ff0c33d0bb5cce7f3e0d55fd2946"
        },
        "date": 1768727173642,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 90114364,
            "unit": "ns",
            "extra": "gctime=12678640\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=334604075\nmemory=4138027152\nallocs=112014837\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 684126742,
            "unit": "ns",
            "extra": "gctime=100731148\nmemory=2189214912\nallocs=3278114\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 59851,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40966,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "c693a5b0747bbeeb9e63947a030a1bfc67aab2cd",
          "message": "bump Optimization (#231)\n\n* bump Optimization\n\n* bump MTK to v10\n\n* fix(tests): adjust variable ordering for MTK v10 compatibility\n\n* style(tests): remove unnecessary blank lines in sgMAM.jl",
          "timestamp": "2026-01-18T13:32:25+01:00",
          "tree_id": "4f9e99e3c4c343e27e55f9813da3dc3d50412e90",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/c693a5b0747bbeeb9e63947a030a1bfc67aab2cd"
        },
        "date": 1768739758407,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 88401232,
            "unit": "ns",
            "extra": "gctime=11729852\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=327740426\nmemory=4137128784\nallocs=111994292\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 674861017,
            "unit": "ns",
            "extra": "gctime=93458475\nmemory=2189277904\nallocs=3279124\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 58439,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40295,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19688\nallocs=12\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "8f5adefaf230839c035d52ba9a4658672cb31fdf",
          "message": "refactor: simple geometric minimal action method (#191)\n\n* refactor: simple geometric minimal action method\n\n* format\n\n* rename ExtendedHamiltonianSystem to ExtendedPhaseSpace\n\n* update docstring\n\n* fix benchmarks",
          "timestamp": "2026-01-21T18:01:28+01:00",
          "tree_id": "035f01d996c52aaaf684a1ad6eee32d92178dc4d",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/8f5adefaf230839c035d52ba9a4658672cb31fdf"
        },
        "date": 1769015237819,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 88852228,
            "unit": "ns",
            "extra": "gctime=11525624\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=337803416\nmemory=4137128784\nallocs=111994292\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 676953945,
            "unit": "ns",
            "extra": "gctime=94041394\nmemory=2189234000\nallocs=3277744\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 58389,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40064,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "872ac0be45981a994d5f54308b5cacad646bde85",
          "message": "fix: vendeheydenMethod (#239)\n\n* fix: vendeheydenMethod\n\n* update docs and bump version number\n\n* fix: emove redundant verbosity tests in CoupledSDEs\n\n---------\n\nCo-authored-by: reykboerner <r.borner@uu.nl>",
          "timestamp": "2026-01-21T18:29:47+01:00",
          "tree_id": "38ddced349c6df2cc16c0ae98ea4b44b8823c595",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/872ac0be45981a994d5f54308b5cacad646bde85"
        },
        "date": 1769016809356,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 90139378,
            "unit": "ns",
            "extra": "gctime=11840049\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=340053884\nmemory=4078989744\nallocs=111994292\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 704445979,
            "unit": "ns",
            "extra": "gctime=109548676\nmemory=2189240368\nallocs=3278614\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 58269,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40416,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "9459441d0ea377f754c7b69f5014d8aba37bb237",
          "message": "fix: update sgMAM method to reduce allocations (#242)",
          "timestamp": "2026-01-21T20:05:01+01:00",
          "tree_id": "64664e26ff96fc845e4403bca68e8abf3a68e0db",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/9459441d0ea377f754c7b69f5014d8aba37bb237"
        },
        "date": 1769022563919,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 92048034.5,
            "unit": "ns",
            "extra": "gctime=13134915\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=9740088\nmemory=102759072\nallocs=2006016\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=420042860\nmemory=4078989744\nallocs=111994292\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 736403298,
            "unit": "ns",
            "extra": "gctime=121412482.5\nmemory=2172404944\nallocs=3257866\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 58599,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40385,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "411465f19eb9a9778990fa4224875a8466795aad",
          "message": "improve performance HeymannVandenEijnden (#244)\n\n* feat: add HeymannVandenEijnden method to Maier-Stein benchmark\n\n* improve performance HeymannVandenEijnden\n\n* format",
          "timestamp": "2026-01-21T21:37:54+01:00",
          "tree_id": "1109047fcf02f76855addb946b3bc6baf6124a0d",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/411465f19eb9a9778990fa4224875a8466795aad"
        },
        "date": 1769028111360,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 87920966,
            "unit": "ns",
            "extra": "gctime=11531026\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37759664\nallocs=765104\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=337509360\nmemory=4078989744\nallocs=111994292\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 667299721,
            "unit": "ns",
            "extra": "gctime=96655748\nmemory=2172463504\nallocs=3258026\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 58419,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40225,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "2772d81fab44ecdfa48cdf630db8523e0b5b21b4",
          "message": "docs: references for string method and geometric minimum action method (#248)",
          "timestamp": "2026-01-22T00:17:10+01:00",
          "tree_id": "9c853b99440d99d59502de709e7f5854b5a4ea06",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/2772d81fab44ecdfa48cdf630db8523e0b5b21b4"
        },
        "date": 1769037863108,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 87837226,
            "unit": "ns",
            "extra": "gctime=11694686.5\nmemory=220493552\nallocs=410229\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=5241677\nmemory=37759664\nallocs=765104\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=340578244\nmemory=4078989744\nallocs=111994292\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 680268078,
            "unit": "ns",
            "extra": "gctime=97175603\nmemory=2172404416\nallocs=3257656\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 58159.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40095,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "4c236ddb492213852881ec377c1482b98dc25e07",
          "message": "feat: update string method to sciml integrator (#247)\n\n* feat: update string method to sciml integrator\n\n* format\n\n* update benchmarks\n\n* format\n\n* only bench euler",
          "timestamp": "2026-01-22T09:43:43+01:00",
          "tree_id": "ce9206ee75d030b985e4e579659b70e5a3cea245",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/4c236ddb492213852881ec377c1482b98dc25e07"
        },
        "date": 1769071662858,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 96183257,
            "unit": "ns",
            "extra": "gctime=12217506\nmemory=222466432\nallocs=440738\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37759664\nallocs=765104\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=309413632\nmemory=4078989744\nallocs=111994292\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 655778398.5,
            "unit": "ns",
            "extra": "gctime=87818854\nmemory=2172477664\nallocs=3258116\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 58348,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40235,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "1c26a2096c454768f8a589d1044a0467a60e79fa",
          "message": "test: add LC action test (#250)",
          "timestamp": "2026-01-22T12:13:39+01:00",
          "tree_id": "86c8d5fc4071f5366ae3bbdc1269ad923f7f5ed2",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/1c26a2096c454768f8a589d1044a0467a60e79fa"
        },
        "date": 1769080651125,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 96212608,
            "unit": "ns",
            "extra": "gctime=12315656\nmemory=222466432\nallocs=440738\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37759664\nallocs=765104\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=315645017\nmemory=4078989744\nallocs=111994292\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 657810128,
            "unit": "ns",
            "extra": "gctime=87998298\nmemory=2172431072\nallocs=3257146\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 58449,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40235,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "4a822d7536817bf2d62dc8c1e5c645a3769676f6",
          "message": "feat(perf): optimize linear solve operations in gmam and sgmam (#251)\n\n* feat(perf): optimize linear solve operations in geometric minimum action method and sgMAM\n\n* format",
          "timestamp": "2026-01-22T16:51:55+01:00",
          "tree_id": "692a0230efb66da85d8b16605f9144f2d8bad400",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/4a822d7536817bf2d62dc8c1e5c645a3769676f6"
        },
        "date": 1769097466612,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 102340377.5,
            "unit": "ns",
            "extra": "gctime=13606559\nmemory=222466432\nallocs=440738\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37695664\nallocs=767104\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=343271224\nmemory=4078989744\nallocs=111994292\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 276834640.5,
            "unit": "ns",
            "extra": "gctime=43609258\nmemory=711307312\nallocs=1509598\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 59621,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 42379,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "e8e899d1a8702ef1ba152ceb66cdbb1ceb1428e6",
          "message": "Refactor large deviation (#257)\n\n* refactor: uniformize API and output for action minimizers\n\n* update docs and examples\n\n* rename venheyden method to geometricgradient\n\n* format and update benchmarks\n\n* update docstrings\n\n* fix benchmarks\n\n* format",
          "timestamp": "2026-01-23T16:21:15+01:00",
          "tree_id": "97fcf531b515b6fad7f17aa01064be066052e6eb",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/e8e899d1a8702ef1ba152ceb66cdbb1ceb1428e6"
        },
        "date": 1769182115624,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 109726727,
            "unit": "ns",
            "extra": "gctime=18166375\nmemory=222466480\nallocs=440740\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=7622196\nmemory=37690480\nallocs=766868\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=7106255\nmemory=37690480\nallocs=766868\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 316469501,
            "unit": "ns",
            "extra": "gctime=65319215\nmemory=711216256\nallocs=1509346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 66013,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40505,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "b77c1f46b308a2891c533882309c4b6b9c96b21e",
          "message": "feat: add converged option to GradientDescent (#259)\n\n* feat: add converged option to GradientDescent\n\n* format\n\n* remove duplicate benchmark entry for Maier-Stein (HeymannVandenEijnden)\n\n* fix benchmarks?\n\n* format",
          "timestamp": "2026-01-24T10:11:11+01:00",
          "tree_id": "1ea0fd4abfa27770dc6f72d00c18d32b2ff5db88",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/b77c1f46b308a2891c533882309c4b6b9c96b21e"
        },
        "date": 1769246172352,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 96840025.5,
            "unit": "ns",
            "extra": "gctime=11677569\nmemory=222466480\nallocs=440740\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 34109467,
            "unit": "ns",
            "extra": "gctime=5247746.5\nmemory=37690480\nallocs=766868\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 2242463904,
            "unit": "ns",
            "extra": "gctime=5303982\nmemory=37690480\nallocs=766868\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 265237622,
            "unit": "ns",
            "extra": "gctime=39620328\nmemory=711385360\nallocs=1509814\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 59923,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 40641.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "c03ceb2fc4a76e18275deb8d2ac8dd3255a5b23b",
          "message": "Add verbose option to sgmam (#265)",
          "timestamp": "2026-01-28T12:33:20+01:00",
          "tree_id": "f08779d482159f164d5583cb79c424ac582a9312",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/c03ceb2fc4a76e18275deb8d2ac8dd3255a5b23b"
        },
        "date": 1769600427062,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 115000102,
            "unit": "ns",
            "extra": "gctime=20935728\nmemory=222466480\nallocs=440740\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 37104992,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37690480\nallocs=766868\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1817176098,
            "unit": "ns",
            "extra": "gctime=422587876.5\nmemory=3447340160\nallocs=83292138\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 323129760.5,
            "unit": "ns",
            "extra": "gctime=65777506.5\nmemory=711433056\nallocs=1509946\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 74104.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 54591.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "7ecd256db3d77dabe442f8360578e1073b8f8f3c",
          "message": "build(deps): pin StochasticDiffEq version :/ (#286)",
          "timestamp": "2026-03-20T08:26:47+01:00",
          "tree_id": "cb4a3cc917b81ec4cdcab8e79fdd16f75f8edb9a",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/7ecd256db3d77dabe442f8360578e1073b8f8f3c"
        },
        "date": 1773992036695,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 98063242,
            "unit": "ns",
            "extra": "gctime=12615656\nmemory=222465392\nallocs=440737\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 35932879.5,
            "unit": "ns",
            "extra": "gctime=2737819\nmemory=37690464\nallocs=766868\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1682091950,
            "unit": "ns",
            "extra": "gctime=270879181.5\nmemory=3447341904\nallocs=83292161\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 284196277.5,
            "unit": "ns",
            "extra": "gctime=42611262\nmemory=711346336\nallocs=1509706\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 65662,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 39734,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "cb847a3949d472af6e9da277fc9152b562230d5c",
          "message": "build(deps): bump julia-actions/cache from 2 to 3 (#281)\n\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>\nCo-authored-by: Orjan Ameye <orjan.ameye@hotmail.com>",
          "timestamp": "2026-03-20T08:27:05+01:00",
          "tree_id": "3f585e754cd270e292f06ee3bf3ef5a779b755c8",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/cb847a3949d472af6e9da277fc9152b562230d5c"
        },
        "date": 1773992064965,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 99802651,
            "unit": "ns",
            "extra": "gctime=13228243\nmemory=222465392\nallocs=440737\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 36904958,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37690464\nallocs=766868\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1747983066.5,
            "unit": "ns",
            "extra": "gctime=290740619.5\nmemory=3447341904\nallocs=83292161\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 288322937,
            "unit": "ns",
            "extra": "gctime=42848739\nmemory=711359344\nallocs=1509742\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 63628,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 42408,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "02b0c7ceaa178d5ea262de8021edf99a8c8a6957",
          "message": "refactor and rename: minimum action methods with type based optimizers  (#267)\n\n* refactor: type based optimizers\n\n* format\n\n* stepsize inside optimizer\n\n* fix ambiquity\n\n* update related methods and documentation\n\n* fix benchmarks\n\n* update docstring, rename points -> npoints\n\n* rename minimizer methods\n\n* fix npoints\n\n---------\n\nCo-authored-by: reykboerner <r.borner@uu.nl>",
          "timestamp": "2026-03-28T02:35:59+01:00",
          "tree_id": "6af518941bade72f02e163cf3c5c62fb26d88473",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/02b0c7ceaa178d5ea262de8021edf99a8c8a6957"
        },
        "date": 1774662192759,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 116090312.5,
            "unit": "ns",
            "extra": "gctime=21420311\nmemory=222465392\nallocs=440737\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 36823339,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37690464\nallocs=766868\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1816070678,
            "unit": "ns",
            "extra": "gctime=434382261\nmemory=3447341904\nallocs=83292161\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 332664366,
            "unit": "ns",
            "extra": "gctime=72424184\nmemory=711411376\nallocs=1509886\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 75813,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 54186,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":20,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "fb5ee28003fe06933225c35119be421375c2e534",
          "message": "build(deps): bump julia-actions/setup-julia from 2 to 3 (#299)\n\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2026-04-22T08:13:56+02:00",
          "tree_id": "6230bfe5112eedc99b1db8d6ee86e8859da20f07",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/fb5ee28003fe06933225c35119be421375c2e534"
        },
        "date": 1776838917620,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 36668022,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37690464\nallocs=766868\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1751648700,
            "unit": "ns",
            "extra": "gctime=289726426.5\nmemory=3447341904\nallocs=83292161\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 280735022.5,
            "unit": "ns",
            "extra": "gctime=45201474\nmemory=711376688\nallocs=1509790\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 99885416,
            "unit": "ns",
            "extra": "gctime=12924469\nmemory=222466448\nallocs=440753\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Rate System/trajectory/Hard coded",
            "value": 66024,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":20,\"time_tolerance\":0.05}"
          },
          {
            "name": "Rate System/trajectory/RateSystem",
            "value": 46608,
            "unit": "ns",
            "extra": "gctime=0\nmemory=19448\nallocs=5\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":20,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "b6fc82b02cfc47e83e99306cae7318ca14e20eb6",
          "message": "fix: unbreak CI on main (#310)",
          "timestamp": "2026-05-01T14:36:08+02:00",
          "tree_id": "94b3010ae409113115c3f34fcbe54397268a27da",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/b6fc82b02cfc47e83e99306cae7318ca14e20eb6"
        },
        "date": 1777639391745,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 36692760,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37690464\nallocs=766868\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1768699427,
            "unit": "ns",
            "extra": "gctime=313648785.5\nmemory=3447341904\nallocs=83292161\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 281262395.5,
            "unit": "ns",
            "extra": "gctime=46125992\nmemory=711229264\nallocs=1509382\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 99459426,
            "unit": "ns",
            "extra": "gctime=12970066\nmemory=222466448\nallocs=440753\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "14fb001c4801899fd3e4a88e6e858146cc9bef2d",
          "message": "format: Move to runic.jl (#315)",
          "timestamp": "2026-05-02T18:10:20+02:00",
          "tree_id": "bc37030ea1187e2f8b4726bff17567c4353cfd7d",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/14fb001c4801899fd3e4a88e6e858146cc9bef2d"
        },
        "date": 1777738497133,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 36759760,
            "unit": "ns",
            "extra": "gctime=0\nmemory=37690464\nallocs=766868\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1788924209,
            "unit": "ns",
            "extra": "gctime=295395888.5\nmemory=3447341904\nallocs=83292161\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 276657014,
            "unit": "ns",
            "extra": "gctime=42769275\nmemory=711268288\nallocs=1509490\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 98939685.5,
            "unit": "ns",
            "extra": "gctime=12783591\nmemory=222466448\nallocs=440753\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "874a85b03a542272c5d3d4dcfee352a9e76d4268",
          "message": "feat: backtracking step-size control for GeometricGradient (#268)",
          "timestamp": "2026-05-11T14:55:42+02:00",
          "tree_id": "24dd9023109b97da53405cf0e86db7eaee77f8e7",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/874a85b03a542272c5d3d4dcfee352a9e76d4268"
        },
        "date": 1778504552255,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 3584300,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4672416\nallocs=102639\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1952684817,
            "unit": "ns",
            "extra": "gctime=430138762\nmemory=3447341904\nallocs=83292161\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 262662841,
            "unit": "ns",
            "extra": "gctime=42744563\nmemory=482908928\nallocs=1229792\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 105899222.5,
            "unit": "ns",
            "extra": "gctime=16724810\nmemory=222466448\nallocs=440753\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "198982749+Copilot@users.noreply.github.com",
            "name": "Copilot",
            "username": "Copilot"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "cc18acd3b3ecc34ed3023335dcce50b6849f027c",
          "message": "Remove legacy `StochSystem` alias from public API (#321)\n\n* Initial plan\n\n* Remove legacy StochSystem alias and update references\n\nAgent-Logs-Url: https://github.com/JuliaDynamics/CriticalTransitions.jl/sessions/0d064a18-953e-41c3-9083-e3701c36811e\n\nCo-authored-by: oameye <57623933+oameye@users.noreply.github.com>\n\n* remove legacy test\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: oameye <57623933+oameye@users.noreply.github.com>\nCo-authored-by: reykboerner <r.borner@uu.nl>",
          "timestamp": "2026-05-12T19:19:14+01:00",
          "tree_id": "db29b73dc7285807eda7f464392376164594f101",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/cc18acd3b3ecc34ed3023335dcce50b6849f027c"
        },
        "date": 1778610387794,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 2771950,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4611360\nallocs=101127\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1451949439,
            "unit": "ns",
            "extra": "gctime=336251103\nmemory=3344629600\nallocs=80748578\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 252032770,
            "unit": "ns",
            "extra": "gctime=47483005.5\nmemory=482804864\nallocs=1229504\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 102063879,
            "unit": "ns",
            "extra": "gctime=15071545\nmemory=222466448\nallocs=440753\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "044f6adb86e6c28a7563a2598e2837b2e584cb9c",
          "message": "build(deps): bump StochasticDiffEq (#323)\n\n* build(deps): bump StochasticDiffEq\n\n* fix tests\n\n* refactor: remove unused StochasticSystemsBase alias and related diffusion_matrix assignment",
          "timestamp": "2026-05-13T12:50:23+02:00",
          "tree_id": "a8434d4b2c201bd00a7919bce5c24a614ac0051b",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/044f6adb86e6c28a7563a2598e2837b2e584cb9c"
        },
        "date": 1778669860811,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 2722056.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4611360\nallocs=101127\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1365874019,
            "unit": "ns",
            "extra": "gctime=280660965\nmemory=3344629600\nallocs=80748578\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 231222565,
            "unit": "ns",
            "extra": "gctime=33628354\nmemory=482913264\nallocs=1229804\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 100340349.5,
            "unit": "ns",
            "extra": "gctime=13372589\nmemory=222466448\nallocs=440753\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "ba2c086a209d31dbda5c8012fd05b6440d13eb04",
          "message": "build(deps): bump  OrdinaryDiffEqLowOrderRK, StochasticDiffEq, .. (#330)",
          "timestamp": "2026-05-14T09:46:24+02:00",
          "tree_id": "32c31f3719b642408612469406abf529a88de6c1",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/ba2c086a209d31dbda5c8012fd05b6440d13eb04"
        },
        "date": 1778745196967,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 2724418,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4611360\nallocs=101128\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1358027509,
            "unit": "ns",
            "extra": "gctime=276448464\nmemory=3344629600\nallocs=80748578\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 229693575.5,
            "unit": "ns",
            "extra": "gctime=34493292\nmemory=482917600\nallocs=1229816\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 100385573,
            "unit": "ns",
            "extra": "gctime=13317963\nmemory=222465520\nallocs=440739\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "9322d31f34d019f6ed4ce5d34a20a320b8a29eeb",
          "message": "fix: update DynamicalSystemsBase and remove StochasticSystemsBase (#335)",
          "timestamp": "2026-05-18T11:12:34+02:00",
          "tree_id": "d1d16e16e4513e0daa38bba4fc31e5603d852ba0",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/9322d31f34d019f6ed4ce5d34a20a320b8a29eeb"
        },
        "date": 1779095889755,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 3079754,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4597648\nallocs=100628\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1605094162.5,
            "unit": "ns",
            "extra": "gctime=447230816.5\nmemory=3321373984\nallocs=79900717\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 281522674,
            "unit": "ns",
            "extra": "gctime=53242713\nmemory=482588064\nallocs=1228904\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 108519972,
            "unit": "ns",
            "extra": "gctime=17011554\nmemory=222465536\nallocs=440739\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "bba2bed4bb74e94de6b07ed1b6859fe148b85b98",
          "message": "feat: adaptive geometric gradient (#334)",
          "timestamp": "2026-05-21T11:42:44+02:00",
          "tree_id": "7d0ebb98f91fb413d2607aabd19609f4db35023c",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/bba2bed4bb74e94de6b07ed1b6859fe148b85b98"
        },
        "date": 1779356897713,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 2738106,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4597648\nallocs=100628\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1319344712.5,
            "unit": "ns",
            "extra": "gctime=273622567\nmemory=3321373984\nallocs=79900717\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 226129518,
            "unit": "ns",
            "extra": "gctime=32348945.5\nmemory=482939280\nallocs=1229876\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 99177683,
            "unit": "ns",
            "extra": "gctime=13097675\nmemory=222465536\nallocs=440739\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "150e6a79f78bff4bcad7c9c2a08411e22465da71",
          "message": "fix: Freidlin-Wentzell action conventions and /2 bugs (#329)",
          "timestamp": "2026-05-21T11:46:12+02:00",
          "tree_id": "ed80d461ae735118d7252e83207190c2bf2e95a9",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/150e6a79f78bff4bcad7c9c2a08411e22465da71"
        },
        "date": 1779357114020,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 2931985,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4597648\nallocs=100628\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1535878731,
            "unit": "ns",
            "extra": "gctime=414579538\nmemory=3321373984\nallocs=79900717\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 269391611,
            "unit": "ns",
            "extra": "gctime=44357291\nmemory=482358256\nallocs=1228268\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 106087577.5,
            "unit": "ns",
            "extra": "gctime=16822603.5\nmemory=222465536\nallocs=440739\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "3a9a24376570edb5199aa0841e822eec617be6ca",
          "message": "fix(transitions): forward sys.diffeq and reseed per call (#331)",
          "timestamp": "2026-05-21T13:09:14+02:00",
          "tree_id": "ce3e7019f1ef430ddded5b7c3b0be79b971c5e99",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/3a9a24376570edb5199aa0841e822eec617be6ca"
        },
        "date": 1779361922362,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 2713951,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4597648\nallocs=100628\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 1442248644,
            "unit": "ns",
            "extra": "gctime=322813488\nmemory=3321373984\nallocs=79900717\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 239443047.5,
            "unit": "ns",
            "extra": "gctime=38837261.5\nmemory=482956624\nallocs=1229924\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 99584467,
            "unit": "ns",
            "extra": "gctime=13369931.5\nmemory=222465536\nallocs=440739\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "3e9d5c7822d85bea0f89beff39fef77a791b4b66",
          "message": "feat: multiplicative noise (#338)",
          "timestamp": "2026-05-24T17:30:56+02:00",
          "tree_id": "40af9e9ed7ced930a90d50a6164e16daa3f9502f",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/3e9d5c7822d85bea0f89beff39fef77a791b4b66"
        },
        "date": 1779637137965,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 880813,
            "unit": "ns",
            "extra": "gctime=0\nmemory=243472\nallocs=2148\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 356161497,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5642656\nallocs=59694\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative diagonal",
            "value": 3389630,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2113136\nallocs=66677\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative off-diagonal",
            "value": 6328511.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4121648\nallocs=97464\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 151522657,
            "unit": "ns",
            "extra": "gctime=11221364\nmemory=141004880\nallocs=370309\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative diagonal",
            "value": 4300985,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4513920\nallocs=102051\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative off-diagonal",
            "value": 10298318,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9428744\nallocs=96963\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 43133178.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2590896\nallocs=70282\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "acb1ae6e4c3dd6367b090dd17082c57b469ecd6c",
          "message": "feat: discrete diffusion operator on Cartesian grids (#327)",
          "timestamp": "2026-05-25T11:15:24+02:00",
          "tree_id": "8e50f397236324d8139e7c73f155032e1fcd5f2e",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/acb1ae6e4c3dd6367b090dd17082c57b469ecd6c"
        },
        "date": 1779700811724,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 891579,
            "unit": "ns",
            "extra": "gctime=0\nmemory=243472\nallocs=2148\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 351703657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5642656\nallocs=59694\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative diagonal",
            "value": 3361980,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2113136\nallocs=66677\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative off-diagonal",
            "value": 6382648.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4121648\nallocs=97464\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 151085804,
            "unit": "ns",
            "extra": "gctime=11020813\nmemory=141004880\nallocs=370309\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative diagonal",
            "value": 4384336,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4513920\nallocs=102051\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative off-diagonal",
            "value": 10219046,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9428744\nallocs=96963\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 43526543,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2590896\nallocs=70282\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "206213b726e1f3d7166132fe7f37c3022a2bafec",
          "message": "feat: OLIM quasipotential solver (#346)",
          "timestamp": "2026-06-01T19:12:21+02:00",
          "tree_id": "308f3aab654ef068822c12a887bf328b8d027173",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/206213b726e1f3d7166132fe7f37c3022a2bafec"
        },
        "date": 1780334740917,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 879493,
            "unit": "ns",
            "extra": "gctime=0\nmemory=243600\nallocs=2148\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 350518249,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5642960\nallocs=59694\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative diagonal",
            "value": 3389556,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2113264\nallocs=66677\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative off-diagonal",
            "value": 6352264.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4121776\nallocs=97464\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 2D multiplicative 81x81",
            "value": 19713844696,
            "unit": "ns",
            "extra": "gctime=0\nmemory=603392\nallocs=189\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 3D quadratic 11x11x11",
            "value": 687282826,
            "unit": "ns",
            "extra": "gctime=0\nmemory=181008\nallocs=780\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, Maier-Stein 121x81",
            "value": 11619555990,
            "unit": "ns",
            "extra": "gctime=0\nmemory=895328\nallocs=196\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 158565810.5,
            "unit": "ns",
            "extra": "gctime=14569743.5\nmemory=141005008\nallocs=370309\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative diagonal",
            "value": 4789432.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4514048\nallocs=102051\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative off-diagonal",
            "value": 8106937,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9428872\nallocs=96963\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 42569602,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2591392\nallocs=70282\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "a127c358b4df6c77f4e6d6177e73baab07f0c9cd",
          "message": "feat(tpt): CTMC formulation of TPT on DiffusionGenerator (#343)",
          "timestamp": "2026-06-01T20:44:39+02:00",
          "tree_id": "d7820c075915e54a133c04a33eb3361ce7f81299",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/a127c358b4df6c77f4e6d6177e73baab07f0c9cd"
        },
        "date": 1780339986344,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 695168,
            "unit": "ns",
            "extra": "gctime=0\nmemory=243600\nallocs=2148\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 279197668.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5642960\nallocs=59694\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative diagonal",
            "value": 2686209.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2113264\nallocs=66677\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative off-diagonal",
            "value": 5005096,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4121776\nallocs=97464\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 2D multiplicative 81x81",
            "value": 14741742456,
            "unit": "ns",
            "extra": "gctime=0\nmemory=603392\nallocs=189\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 3D quadratic 11x11x11",
            "value": 528294293,
            "unit": "ns",
            "extra": "gctime=0\nmemory=181008\nallocs=780\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, Maier-Stein 121x81",
            "value": 8606861472.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=895328\nallocs=196\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 122168647,
            "unit": "ns",
            "extra": "gctime=10163733\nmemory=141005008\nallocs=370309\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative diagonal",
            "value": 3567648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4035616\nallocs=90918\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative off-diagonal",
            "value": 6297032.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9428872\nallocs=96963\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 34858568,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2591392\nallocs=70282\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "1a50e967673d28a084d5673174d3bb2c5135a3b1",
          "message": "feat: multiple-shooting BVP for the Freidlin-Wentzell instanton (#347)",
          "timestamp": "2026-06-01T21:56:49+02:00",
          "tree_id": "7a541b70683efa199a795e8571f444e82221730a",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/1a50e967673d28a084d5673174d3bb2c5135a3b1"
        },
        "date": 1780344347322,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 848411.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=243248\nallocs=2145\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 339180774.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5643552\nallocs=59696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative diagonal",
            "value": 3310452,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2112912\nallocs=66674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative off-diagonal",
            "value": 6221954,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4121424\nallocs=97461\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 2D multiplicative 81x81",
            "value": 19069997635.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=603392\nallocs=189\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 3D quadratic 11x11x11",
            "value": 651879303.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=181008\nallocs=780\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, Maier-Stein 121x81",
            "value": 11917714733,
            "unit": "ns",
            "extra": "gctime=0\nmemory=895328\nallocs=196\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 157869631,
            "unit": "ns",
            "extra": "gctime=11164941\nmemory=141006688\nallocs=370311\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative diagonal",
            "value": 4683144,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4514768\nallocs=102053\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative off-diagonal",
            "value": 8714541,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9429928\nallocs=96965\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 42400879,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2593072\nallocs=70284\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "d9c50c1ffd074372ac283c6bef5c98489c5c1212",
          "message": "feat: degenerate noise OLIM (rank-1 diffusion) (#351)",
          "timestamp": "2026-06-03T15:39:37+02:00",
          "tree_id": "587b4e8da2bc7f20f5a8bca6fe5222e623b33e0e",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/d9c50c1ffd074372ac283c6bef5c98489c5c1212"
        },
        "date": 1780494623560,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 869674.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=243248\nallocs=2145\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 350692948,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5643552\nallocs=59696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative diagonal",
            "value": 3380366,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2112912\nallocs=66674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative off-diagonal",
            "value": 6315633,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4121424\nallocs=97461\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 2D multiplicative 81x81",
            "value": 19500851714.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=605472\nallocs=209\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 3D quadratic 11x11x11",
            "value": 612785756,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237176\nallocs=806\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, Maier-Stein 121x81",
            "value": 9943359453.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1211320\nallocs=222\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 150112835,
            "unit": "ns",
            "extra": "gctime=10501574\nmemory=141006688\nallocs=370311\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative diagonal",
            "value": 4258913,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4514768\nallocs=102053\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative off-diagonal",
            "value": 8020044,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9429928\nallocs=96965\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 42954620.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2593072\nallocs=70284\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "b374dd6d1be5761ee0d8d593600d499158fef395",
          "message": "docs(build): answer some TODO's (#344)",
          "timestamp": "2026-06-03T19:51:32+02:00",
          "tree_id": "474e76272b1599dd59c5dc06a914ea3fab8f463f",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/b374dd6d1be5761ee0d8d593600d499158fef395"
        },
        "date": 1780509702459,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 889412,
            "unit": "ns",
            "extra": "gctime=0\nmemory=243248\nallocs=2145\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 358073255.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5643552\nallocs=59696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative diagonal",
            "value": 3431185.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2112912\nallocs=66674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative off-diagonal",
            "value": 6366792,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4121424\nallocs=97461\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 2D multiplicative 81x81",
            "value": 19512475930.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=605472\nallocs=209\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 3D quadratic 11x11x11",
            "value": 609256612,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237176\nallocs=806\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, Maier-Stein 121x81",
            "value": 9880992555.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1211320\nallocs=222\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 155422442,
            "unit": "ns",
            "extra": "gctime=12410892\nmemory=141006688\nallocs=370311\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative diagonal",
            "value": 4445687,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4514768\nallocs=102053\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative off-diagonal",
            "value": 8118746,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9429912\nallocs=96965\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 42981765.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2593072\nallocs=70284\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "88832655+raphael-roemer@users.noreply.github.com",
            "name": "raphael-roemer",
            "username": "raphael-roemer"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e7c49fbf0e41b1429b4c2c69534f3c1947e0d23e",
          "message": "Finalise RateSystem (#316)\n\n* first go at implementing the agreed TODOs (not yet double-checked)\n\n* fixed bugs, should now run as intended\n\n* updated docstrings and implemented last TODO\n\n* tidied up AND added tests\n\n* added CoupledSDEs functionality and corresponding testset\n\n* spell check fixes AND extended SDE testset\n\n* rename forcers, add 1-parameter RateSystem constructor\n\n* fix syntax bug\n\n* added and restructured tests\n\n* fixed bug in new testsets\n\n* included ModelingToolkit again\n\n* added 2-parameter-shifts tests that are exactly like 1-parameter cases, just split up into 2parameters\n\n* included ModelinToolkit again\n\n* add API function for ForcingProfile\n\n* fix parameters function\n\n* move ForcingProfile code to separate file\n\n* update RateSystem\n\n* still WIP\n\n* implement update_parameters! with set_parameters! call\n\n* remove function with try statement\n\n* implemented review comments\n\n* updated RateSystemSpec constructor\n\n* reverted last commit\n\n* Merge main branch into FinaliseRateSystem branch for consistency (#336)\n\n* build(deps): bump crate-ci/typos from 1.45.2 to 1.46.0 (#317)\n\nBumps [crate-ci/typos](https://github.com/crate-ci/typos) from 1.45.2 to 1.46.0.\n- [Release notes](https://github.com/crate-ci/typos/releases)\n- [Changelog](https://github.com/crate-ci/typos/blob/master/CHANGELOG.md)\n- [Commits](https://github.com/crate-ci/typos/compare/v1.45.2...v1.46.0)\n\n---\nupdated-dependencies:\n- dependency-name: crate-ci/typos\n  dependency-version: 1.46.0\n  dependency-type: direct:production\n  update-type: version-update:semver-minor\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>\n\n* feat: backtracking step-size control for GeometricGradient (#268)\n\n* Remove legacy `StochSystem` alias from public API (#321)\n\n* Initial plan\n\n* Remove legacy StochSystem alias and update references\n\nAgent-Logs-Url: https://github.com/JuliaDynamics/CriticalTransitions.jl/sessions/0d064a18-953e-41c3-9083-e3701c36811e\n\nCo-authored-by: oameye <57623933+oameye@users.noreply.github.com>\n\n* remove legacy test\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: oameye <57623933+oameye@users.noreply.github.com>\nCo-authored-by: reykboerner <r.borner@uu.nl>\n\n* build(deps): bump crate-ci/typos from 1.46.0 to 1.46.1 (#322)\n\nBumps [crate-ci/typos](https://github.com/crate-ci/typos) from 1.46.0 to 1.46.1.\n- [Release notes](https://github.com/crate-ci/typos/releases)\n- [Changelog](https://github.com/crate-ci/typos/blob/master/CHANGELOG.md)\n- [Commits](https://github.com/crate-ci/typos/compare/v1.46.0...v1.46.1)\n\n---\nupdated-dependencies:\n- dependency-name: crate-ci/typos\n  dependency-version: 1.46.1\n  dependency-type: direct:production\n  update-type: version-update:semver-patch\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>\n\n* build(deps): bump StochasticDiffEq (#323)\n\n* build(deps): bump StochasticDiffEq\n\n* fix tests\n\n* refactor: remove unused StochasticSystemsBase alias and related diffusion_matrix assignment\n\n* build(deps): bump  OrdinaryDiffEqLowOrderRK, StochasticDiffEq, .. (#330)\n\n* fix: update DynamicalSystemsBase and remove StochasticSystemsBase (#335)\n\n---------\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>\nCo-authored-by: Orjan Ameye <orjan.ameye@hotmail.com>\nCo-authored-by: Copilot <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: oameye <57623933+oameye@users.noreply.github.com>\nCo-authored-by: reykboerner <r.borner@uu.nl>\n\n* updated docstrings in RateSystem.jl\n\n* removed StochasticSystemsBase as this extension was removed from src/sde_utlis.jl\n\n* Restore src/largedeviations/sgmam.jl from origin/main (fix missing include)\n\n* trying to fix tests again\n\n* typofix in referenced_sciml_prob()\n\n* updated RateSystem docs via updating examples/RateSystem.jl and regenerating .md file\n\n* updated docs\n\n* hopefully fixed tests\n\n* exported parameter\n\n* fix wrong parameterization of Specs; fix missing copy\n\n* WIP on updating tutorial and making it work\n\n* finito\n\n* remove obsolete rate system examples\n\n* finish docs for rate system\n\n---------\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: reykboerner <r.borner@uu.nl>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>\nCo-authored-by: Orjan Ameye <orjan.ameye@hotmail.com>\nCo-authored-by: Copilot <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: oameye <57623933+oameye@users.noreply.github.com>\nCo-authored-by: Datseris <datseris.george@gmail.com>",
          "timestamp": "2026-06-03T21:58:46+01:00",
          "tree_id": "f05c14b4988c6eff19a5ec38203b8ef160bb6c4b",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/e7c49fbf0e41b1429b4c2c69534f3c1947e0d23e"
        },
        "date": 1780521160172,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 754836,
            "unit": "ns",
            "extra": "gctime=0\nmemory=243248\nallocs=2145\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 272590760,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5643552\nallocs=59696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative diagonal",
            "value": 3276333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2112912\nallocs=66674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative off-diagonal",
            "value": 6036280,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4121424\nallocs=97461\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 2D multiplicative 81x81",
            "value": 13995509973,
            "unit": "ns",
            "extra": "gctime=0\nmemory=605472\nallocs=209\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 3D quadratic 11x11x11",
            "value": 522289308.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237176\nallocs=806\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, Maier-Stein 121x81",
            "value": 9090628812,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1211320\nallocs=222\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 162642620.5,
            "unit": "ns",
            "extra": "gctime=15596317\nmemory=141006688\nallocs=370311\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative diagonal",
            "value": 3724387,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4036336\nallocs=90920\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative off-diagonal",
            "value": 8412164,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9429928\nallocs=96965\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 43467193.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2593072\nallocs=70284\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "aaabe768eba780f2676251489b91126208627aea",
          "message": "cleanup #316 aftermath (#356)",
          "timestamp": "2026-06-04T11:31:16+02:00",
          "tree_id": "4eb58f0e82b9eac2b2032eba20d53372c31c23a1",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/aaabe768eba780f2676251489b91126208627aea"
        },
        "date": 1780566267979,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 842822.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=243248\nallocs=2145\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 337776545.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5643552\nallocs=59696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative diagonal",
            "value": 3301936.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2112912\nallocs=66674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative off-diagonal",
            "value": 6144501,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4121424\nallocs=97461\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 2D multiplicative 81x81",
            "value": 18792031427,
            "unit": "ns",
            "extra": "gctime=0\nmemory=605472\nallocs=209\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 3D quadratic 11x11x11",
            "value": 586489169,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237176\nallocs=806\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, Maier-Stein 121x81",
            "value": 10439710665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1211320\nallocs=222\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 152314440,
            "unit": "ns",
            "extra": "gctime=8872741\nmemory=141006688\nallocs=370311\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative diagonal",
            "value": 4476158.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4514768\nallocs=102053\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative off-diagonal",
            "value": 8149156,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9429912\nallocs=96965\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 42090839,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2593072\nallocs=70284\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "6761c18aa3813ac09083e578badc661fe716b34f",
          "message": "docs: rewrite the manual (#355)",
          "timestamp": "2026-06-04T12:43:17+02:00",
          "tree_id": "e6bd502895f182ff07f85bd76d5d8e6e18a187b1",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/6761c18aa3813ac09083e578badc661fe716b34f"
        },
        "date": 1780570334193,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 844376,
            "unit": "ns",
            "extra": "gctime=0\nmemory=243248\nallocs=2145\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 337523941.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5643552\nallocs=59696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative diagonal",
            "value": 3293302,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2112912\nallocs=66674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative off-diagonal",
            "value": 6177667.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4121424\nallocs=97461\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 2D multiplicative 81x81",
            "value": 18725380190.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=605472\nallocs=209\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 3D quadratic 11x11x11",
            "value": 587258855,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237176\nallocs=806\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, Maier-Stein 121x81",
            "value": 10363381833,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1211320\nallocs=222\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 155052245.5,
            "unit": "ns",
            "extra": "gctime=10079240\nmemory=141006688\nallocs=370311\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative diagonal",
            "value": 4422544.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4514768\nallocs=102053\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative off-diagonal",
            "value": 8101540.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9429912\nallocs=96965\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 42790199,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2593072\nallocs=70284\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "datseris.george@gmail.com",
            "name": "George Datseris",
            "username": "Datseris"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9c177cee9526b4bf32c04597200b9045d71286db",
          "message": "Rate-return version 2 (#358)\n\n* first go at implementing the agreed TODOs (not yet double-checked)\n\n* fixed bugs, should now run as intended\n\n* updated docstrings and implemented last TODO\n\n* tidied up AND added tests\n\n* added CoupledSDEs functionality and corresponding testset\n\n* spell check fixes AND extended SDE testset\n\n* rename forcers, add 1-parameter RateSystem constructor\n\n* fix syntax bug\n\n* added and restructured tests\n\n* fixed bug in new testsets\n\n* included ModelingToolkit again\n\n* added 2-parameter-shifts tests that are exactly like 1-parameter cases, just split up into 2parameters\n\n* included ModelinToolkit again\n\n* add API function for ForcingProfile\n\n* fix parameters function\n\n* move ForcingProfile code to separate file\n\n* update RateSystem\n\n* still WIP\n\n* implement update_parameters! with set_parameters! call\n\n* remove function with try statement\n\n* implemented review comments\n\n* updated RateSystemSpec constructor\n\n* reverted last commit\n\n* Merge main branch into FinaliseRateSystem branch for consistency (#336)\n\n* build(deps): bump crate-ci/typos from 1.45.2 to 1.46.0 (#317)\n\nBumps [crate-ci/typos](https://github.com/crate-ci/typos) from 1.45.2 to 1.46.0.\n- [Release notes](https://github.com/crate-ci/typos/releases)\n- [Changelog](https://github.com/crate-ci/typos/blob/master/CHANGELOG.md)\n- [Commits](https://github.com/crate-ci/typos/compare/v1.45.2...v1.46.0)\n\n---\nupdated-dependencies:\n- dependency-name: crate-ci/typos\n  dependency-version: 1.46.0\n  dependency-type: direct:production\n  update-type: version-update:semver-minor\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>\n\n* feat: backtracking step-size control for GeometricGradient (#268)\n\n* Remove legacy `StochSystem` alias from public API (#321)\n\n* Initial plan\n\n* Remove legacy StochSystem alias and update references\n\nAgent-Logs-Url: https://github.com/JuliaDynamics/CriticalTransitions.jl/sessions/0d064a18-953e-41c3-9083-e3701c36811e\n\nCo-authored-by: oameye <57623933+oameye@users.noreply.github.com>\n\n* remove legacy test\n\n---------\n\nCo-authored-by: copilot-swe-agent[bot] <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: oameye <57623933+oameye@users.noreply.github.com>\nCo-authored-by: reykboerner <r.borner@uu.nl>\n\n* build(deps): bump crate-ci/typos from 1.46.0 to 1.46.1 (#322)\n\nBumps [crate-ci/typos](https://github.com/crate-ci/typos) from 1.46.0 to 1.46.1.\n- [Release notes](https://github.com/crate-ci/typos/releases)\n- [Changelog](https://github.com/crate-ci/typos/blob/master/CHANGELOG.md)\n- [Commits](https://github.com/crate-ci/typos/compare/v1.46.0...v1.46.1)\n\n---\nupdated-dependencies:\n- dependency-name: crate-ci/typos\n  dependency-version: 1.46.1\n  dependency-type: direct:production\n  update-type: version-update:semver-patch\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>\n\n* build(deps): bump StochasticDiffEq (#323)\n\n* build(deps): bump StochasticDiffEq\n\n* fix tests\n\n* refactor: remove unused StochasticSystemsBase alias and related diffusion_matrix assignment\n\n* build(deps): bump  OrdinaryDiffEqLowOrderRK, StochasticDiffEq, .. (#330)\n\n* fix: update DynamicalSystemsBase and remove StochasticSystemsBase (#335)\n\n---------\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>\nCo-authored-by: Orjan Ameye <orjan.ameye@hotmail.com>\nCo-authored-by: Copilot <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: oameye <57623933+oameye@users.noreply.github.com>\nCo-authored-by: reykboerner <r.borner@uu.nl>\n\n* updated docstrings in RateSystem.jl\n\n* removed StochasticSystemsBase as this extension was removed from src/sde_utlis.jl\n\n* Restore src/largedeviations/sgmam.jl from origin/main (fix missing include)\n\n* trying to fix tests again\n\n* typofix in referenced_sciml_prob()\n\n* updated RateSystem docs via updating examples/RateSystem.jl and regenerating .md file\n\n* updated docs\n\n* hopefully fixed tests\n\n* exported parameter\n\n* fix wrong parameterization of Specs; fix missing copy\n\n* WIP on updating tutorial and making it work\n\n* finito\n\n* remove obsolete rate system examples\n\n* finish docs for rate system\n\n* add reversable tipping parameters\n\n* add reverse example in tutorial\n\n* copy original parameters to rate system\n\n* add forcing reverse function\n\n* fix incorrect forcing of time and parameter types\n\n* update tutorial to have reyks code\n\n* add rtipping phase diagrams (not working)\n\n* revert everything to start and just add `reverse = true`.\n\n* fix tests by commenting out DistMesh2D\n\n* reexport Attractors, Remove unecessary StaticArrays\n\n* almost finish rate phase diagram\n\n* remove unecessary using static arrays\n\n* finish docstring of phase diagrams\n\n* remove sparse arrays more\n\n* finish the source code and tutorial but FUCK it doesn't work again\n\n* use correct steady state\n\n* change the parameter value to bistable\n\n* revert unnecsaary changes from merge commit\n\n* format\n\n* more clean up\n\n* fix explicit imports\n\n* fix docs\n\n* fix LLM review\n\n* add current_time to DynamicalSystemsBase imports\n\n---------\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: raphael-roemer <raphael-roemer@t-online.de>\nCo-authored-by: reykboerner <r.borner@uu.nl>\nCo-authored-by: raphael-roemer <88832655+raphael-roemer@users.noreply.github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>\nCo-authored-by: Orjan Ameye <orjan.ameye@hotmail.com>\nCo-authored-by: Copilot <198982749+Copilot@users.noreply.github.com>\nCo-authored-by: oameye <57623933+oameye@users.noreply.github.com>",
          "timestamp": "2026-06-12T09:42:20+02:00",
          "tree_id": "a04988cbf61143a8a6ed1fcf3ee3ee745ba75eed",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/9c177cee9526b4bf32c04597200b9045d71286db"
        },
        "date": 1781250890562,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 698335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=243264\nallocs=2145\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 284896299,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5643552\nallocs=59696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative diagonal",
            "value": 2690505,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2112928\nallocs=66674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative off-diagonal",
            "value": 5001658,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4121440\nallocs=97461\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 2D multiplicative 81x81",
            "value": 14504090677,
            "unit": "ns",
            "extra": "gctime=0\nmemory=605472\nallocs=209\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 3D quadratic 11x11x11",
            "value": 459327723.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237176\nallocs=806\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, Maier-Stein 121x81",
            "value": 7410026296,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1211320\nallocs=222\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 122596868.5,
            "unit": "ns",
            "extra": "gctime=9504608.5\nmemory=141006704\nallocs=370311\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative diagonal",
            "value": 3037225,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4036352\nallocs=90920\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative off-diagonal",
            "value": 6382271,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9429928\nallocs=96965\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 35443511.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2593072\nallocs=70284\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a279b5f8f0cb921139e412c0ce177f0063f62160",
          "message": "build(deps): bump actions/checkout from 6 to 7 (#369)\n\nBumps [actions/checkout](https://github.com/actions/checkout) from 6 to 7.\n- [Release notes](https://github.com/actions/checkout/releases)\n- [Changelog](https://github.com/actions/checkout/blob/main/CHANGELOG.md)\n- [Commits](https://github.com/actions/checkout/compare/v6...v7)\n\n---\nupdated-dependencies:\n- dependency-name: actions/checkout\n  dependency-version: '7'\n  dependency-type: direct:production\n  update-type: version-update:semver-major\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2026-06-24T08:35:27+02:00",
          "tree_id": "84c3e1f1d0a42a2a98a7a248ba1739711ce6e70d",
          "url": "https://github.com/JuliaDynamics/CriticalTransitions.jl/commit/a279b5f8f0cb921139e412c0ce177f0063f62160"
        },
        "date": 1782283782459,
        "tool": "julia",
        "benches": [
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (HeymannVandenEijnden)",
            "value": 874639.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=243264\nallocs=2145\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/Maier-Stein (Optimisers.Adam; AutoFiniteDiff)",
            "value": 337114506,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5643552\nallocs=59696\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative diagonal",
            "value": 3288676,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2112928\nallocs=66674\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Geometric minimal action/multiplicative off-diagonal",
            "value": 6148215,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4121440\nallocs=97461\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 2D multiplicative 81x81",
            "value": 18810259842,
            "unit": "ns",
            "extra": "gctime=0\nmemory=605472\nallocs=209\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, 3D quadratic 11x11x11",
            "value": 587861215,
            "unit": "ns",
            "extra": "gctime=0\nmemory=237176\nallocs=806\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Quasipotential/OLIM, Maier-Stein 121x81",
            "value": 10451245040,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1211320\nallocs=222\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":30,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 151545231,
            "unit": "ns",
            "extra": "gctime=7864106\nmemory=141006704\nallocs=370311\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative diagonal",
            "value": 4228807.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4514784\nallocs=102053\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/multiplicative off-diagonal",
            "value": 7916555,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9429944\nallocs=96965\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Large deviation/String method/Kerr parametric resonator",
            "value": 42954755,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2593120\nallocs=70285\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":10,\"time_tolerance\":0.05}"
          }
        ]
      }
    ]
  }
}