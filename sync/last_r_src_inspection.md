# R src sync inspection

- R repo: OVVO-Financial/NNS
- R commit: `2241afd55be348a628bddc94b53819380a435227`
- R version: `13.0`
- Previous src hash: `03ac013e5cdaa606300eea2b5c0ac16aa23a61362b555d5522e42eccd069d581`
- Incoming src hash: `f3e1a50da9f80c94be44ec1046f805b8b4fa00f2d67a83944cd3d4c835995425`
- Changed files:
  - `DESCRIPTION`
  - `NNS_13.0.tar.gz`
  - `NNS_13.0.zip`
  - `R/NNS_meboot.R`
  - `doc/NNSvignette_05_Sampling.R`
  - `doc/NNSvignette_05_Sampling.Rmd`
  - `doc/NNSvignette_05_Sampling.html`
  - `man/NNS.meboot.Rd`
  - `src/NNS.dll`
  - `vignettes/NNSvignette_05_Sampling.R`
  - `vignettes/NNSvignette_05_Sampling.Rmd`
  - `vignettes/NNSvignette_05_Sampling.html`

## Result

`src/**` changed. Review and port the affected R C++ files into NNS-core.

This inspection does not auto-port arbitrary C++ changes. After the core source has actually been updated and C++ build and tests pass, update `sync/nns_r_source.json` with the new `r_commit`, `r_version`, `r_src_tree_hash`, and `core_commit`.
