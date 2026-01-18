window.BENCHMARK_DATA = {
  "lastUpdate": 1768727177274,
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
            "value": 2224884285,
            "unit": "ns",
            "extra": "gctime=328491072\nmemory=4138027152\nallocs=112014837\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Large deviation/Simple geometric minimal action/KPO",
            "value": 672580551,
            "unit": "ns",
            "extra": "gctime=94337914\nmemory=2207190624\nallocs=3378374\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
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
            "value": 2238500481,
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
            "value": 2251822088,
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
            "value": 2365841917,
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
            "value": 2268668879,
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
            "value": 2327762273,
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
            "value": 2250953877,
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
            "value": 2262480032,
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
            "value": 2231509041,
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
          }
        ]
      }
    ]
  }
}