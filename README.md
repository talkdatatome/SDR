# SDR
Formulas for expected values of Successive Difference Replication under Simple Random Sampling without Replacement.

This package supplements the work found in the citation below. While the package itself provides a small vignette validating the equivalency of the three versions of the variance estimation, the separate SDRpkg.RLog file provides further examples to support the assertions made in the paper. On GitHub, it is found in the inst/ folder. Once installed locally as a package, it can be found via the R command `system.file("SDRpkg.RLog", package="SDR")`. Notation can be found in the SDR-package help via `?SDR` in R after loading the package.

## Citation

If you use the SDR package, please cite it as:

> Slud, E., & Trudell, T. (2026). SDR Variance Estimates in Small Domains. arXiv preprint arXiv:2608.17353. https://doi.org/10.48550/arXiv.2608.17353
```bibtex
@misc{slud2026sdrvarianceestimatessmall,
      title={SDR Variance Estimates in Small Domains}, 
      author={Eric Slud and Tim Trudell},
      year={2026},
      eprint={2608.17353},
      archivePrefix={arXiv},
      primaryClass={stat.ME},
      url={https://arxiv.org/abs/2608.17353}, 
}
```


