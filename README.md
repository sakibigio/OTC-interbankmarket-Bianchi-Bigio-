# Portfolio Theory with Settlement Frictions - Replication Code

This repository contains MATLAB code for replicating the numerical results and figures in **"Portfolio Theory with Settlement Frictions"** by Javier Bianchi and Saki Bigio.

## Paper Information

- **Authors**: Javier Bianchi & Saki Bigio
- **Journal**: Journal of Economic Theory (JET)
- **Manuscript ID**: JET-D-23-00798

## Overview

This project analyzes portfolio optimization with settlement frictions in an interbank market model. The code solves for liquidity yield functions under different matching technologies and generates comparative statics across market tightness, matching efficiency, and bargaining power parameters.

## Repository Structure

```
interbankmarket-code/
├── README.md                  # This file
├── LICENSE                    # License information
├── code/                      # All source code
│   ├── interbankmarket_example.m  # Main script to run all analyses
│   ├── solvers/              # Numerical solution algorithms
│   ├── analytics/            # Analytical solution functions
│   ├── utilities/            # Helper functions
│   └── plotting/             # Plotting utilities
├── output/                    # Generated figures (gitignored)
│   └── figures/              
├── paper/                     # Paper and reference documents
│   └── references/           
└── .gitignore                # Git ignore patterns
```

## Requirements

### Software
- **MATLAB** R2020b or later (tested on R2024a)
- Required Toolboxes: None (base MATLAB only)

### System
- Memory: 4 GB RAM minimum
- Disk Space: 500 MB for code and figures

## Quick Start

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/interbankmarket-code.git
   cd interbankmarket-code
   ```

2. **Open MATLAB** and navigate to the code directory:
   ```matlab
   cd code
   ```

3. **Run the main script**:
   ```matlab
   interbankmarket_example
   ```

This will generate all figures in the `output/figures/` directory.

## Code Organization

### Main Script
- **`interbankmarket_example.m`**: Master script that runs all analyses and generates figures

### Solvers (`code/solvers/`)
Numerical solution methods for the interbank market equilibrium:
- `interbanksolve_discrete.m` - Discrete-time solver with finite trading rounds
- `interbanksolve_continuous.m` - Continuous-time limit solver
- `interbanksolve_analyticleontief.m` - Analytical solver for Leontief matching
- `interbanksolve_analyticcobbdouglas.m` - Analytical solver for Cobb-Douglas matching
- `interbanksolve_analyticharmonic.m` - Analytical solver for harmonic mean matching

### Analytics (`code/analytics/`)
Closed-form solutions for specific matching functions:
- `analytic_leontief.m` - Leontief (min) matching function solutions
- `analytic_cobbdouglas.m` - Cobb-Douglas matching function solutions  
- `analytic_harmonic.m` - Harmonic mean matching function solutions
- `analyticleontieff.m` - Alternative Leontief implementation

### Utilities (`code/utilities/`)
Helper functions and parameter management:
- `build_matching.m` - Constructs CES matching function with specified elasticity
- `interbank_resetmatchtech.m` - Resets matching technology parameters
- `Plot_Migration.m` - Auto-sync figures to Overleaf (optional)

### Plotting (`code/plotting/`)
Visualization utilities for figures:
- `interbankplot_addcorridortheta.m` - Add corridor bounds (theta plots)
- `interbankplot_addcorridorround.m` - Add corridor bounds (round plots)
- `interbankplot_addcorridorlambda.m` - Add corridor bounds (lambda plots)
- `interbankplot_addstaticrate.m` - Add static rate reference line

## Key Parameters

The model is calibrated with the following baseline parameters (defined in `interbankmarket_example.m`):

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Discount window rate | r^w | 0.0035 | Annual rate at which banks can borrow from central bank |
| Excess reserves rate | r^m | 0.0025 | Annual rate paid on excess reserves |
| Matching efficiency | λ | 1.0 | Speed of matching in OTC market |
| Bargaining power | η | 0.5 | Borrower's bargaining power in Nash bargaining |
| Trading rounds | N | 10,000 | Number of trading rounds in discrete-time model |

## Figures Generated

The script generates the following figures (saved as PDF in `output/figures/`):

### Matching Function Comparisons
- `F_growthrate.pdf` - Growth rate comparison across matching functions
- `F_chip_comp.pdf` - Chi-plus liquidity yield comparison
- `F_chim_comp.pdf` - Chi-minus liquidity yield comparison
- `F_trajectories.pdf` - Theta evolution trajectories

### Leontief Matching (L-series)
- `F_l_Chiplus_theta.pdf` - Chi-plus as function of theta
- `F_l_Chiminus_theta.pdf` - Chi-minus as function of theta
- `F_l_InterbankRate_theta.pdf` - OTC rate as function of theta
- `F_l_Surplus_theta.pdf` - Surplus as function of theta
- `F_l_Chiplus_tau.pdf` - Chi-plus evolution over trading rounds
- `F_l_Chiminus_tau.pdf` - Chi-minus evolution over trading rounds
- `F_l_InterbankRate_tau.pdf` - OTC rate evolution over trading rounds
- `F_l_prices_lambda.pdf` - Prices as function of matching efficiency
- `F_l_Q_lambda.pdf` - Rate dispersion as function of lambda
- `F_l_Vol_lambda.pdf` - Trading volume as function of lambda

### Cobb-Douglas Matching (CD-series)
- Similar set of figures for Cobb-Douglas matching function

### Harmonic Mean Matching (H-series)
- `F_h_InterbankRate_tau.pdf` - OTC rates over trading rounds
- `F_h_Chiplus_tau.pdf` - Chi-plus over trading rounds
- `F_h_Chiminus_tau.pdf` - Chi-minus over trading rounds
- `F_h_Surplus_tau.pdf` - Surplus over trading rounds
- `F_h_gammaplus_tau.pdf` - Surplus-side matching rates
- `F_h_gammaminus_tau.pdf` - Deficit-side matching rates

### Special Cases
- `F_Walrasian_prices_theta.pdf` - Walrasian limit (high efficiency)

## Customization

### Changing Parameters

To modify baseline parameters, edit the relevant section in `interbankmarket_example.m`:

```matlab
%% Parameters
lambda_o = 1;      % Matching efficiency
eta_o = 0.5;       % Bargaining power
theta_bench = 0.75; % Benchmark market tightness
N = 10000;         % Number of rounds
```

### Selecting Figures to Generate

Control figure generation with flags:

```matlab
plotit = 0;   % Generate extra diagnostic figures (0=no, 1=yes)
printit = 1;  % Save figures to PDF (0=no, 1=yes)
```

### Output Directory

By default, figures are saved to `output/figures/`. To change the output directory, modify:

```matlab
foldername = 'your/custom/path/';  % Must end with '/'
```

## Matching Functions

The code supports several CES matching functions parameterized by elasticity ρ:

| Function | ρ | p = (ρ-1)/ρ | Description |
|----------|---|--------------|-------------|
| Leontief | 0 | -∞ | Perfect complementarity (min function) |
| Harmonic | -1 | -1 | Harmonic mean matching |
| Cobb-Douglas | 1 | 0 | Geometric mean matching |
| Linear | ∞ | 1 | Perfect substitutability |

The matching function is: G(a,b) = λ(1/2)^(ρ/(ρ-1))(a^(1-1/ρ) + b^(1-1/ρ))^(ρ/(ρ-1))

## Mathematical Notation

| Code Variable | Mathematical Symbol | Description |
|---------------|---------------------|-------------|
| `theta_0` | θ₀ | Initial market tightness (deficit/surplus ratio) |
| `CHIp` | χ⁺ | Liquidity yield for surplus positions |
| `CHIm` | χ⁻ | Liquidity yield for deficit positions |
| `delta_r` | r^w - r^m | Spread between discount window and reserve rates |
| `lambda` | λ | Matching efficiency parameter |
| `eta` | η | Borrower's bargaining power |
| `r_f` | r^f | OTC market equilibrium rate |
| `gammas` | γˢ | Matching probability for surplus side |
| `gammad` | γᵈ | Matching probability for deficit side |

## Computational Notes

- **Precision**: All computations use MATLAB's double precision (≈15-17 decimal digits)
- **Runtime**: Full script execution takes approximately 2-5 minutes on modern hardware
- **Convergence**: Analytical solutions are exact; numerical solutions converge as N → ∞

## Troubleshooting

### Common Issues

1. **"Out of memory" error**
   - Reduce `N` (number of rounds) from 10,000 to 1,000
   - Close other applications to free RAM

2. **Figures not saving**
   - Check that output directory exists and is writable
   - Verify `printit = 1` in the script

3. **"Function not found" errors**
   - Ensure MATLAB's current directory includes all subdirectories
   - Run `addpath(genpath('code'))` to add all subdirectories

## Citation

If you use this code, please cite:

```bibtex
@article{bianchi2025portfolio,
  title={Portfolio Theory with Settlement Frictions},
  author={Bianchi, Javier and Bigio, Saki},
  journal={Journal of Economic Theory},
  year={2025},
  note={Manuscript ID: JET-D-23-00798}
}
```

## License

This code is released under the MIT License. See `LICENSE` file for details.

## Contact

For questions about the code or paper:
- **Javier Bianchi**: [email/website]
- **Saki Bigio**: [email/website]

## Acknowledgments

This research was supported by [funding information if applicable].

## Version History

- **v1.0.0** (2025-04-29): Initial release
- Complete implementation of all matching functions
- All figures for JET submission

## Contributing

This is research code for a published paper. For bug reports or suggestions, please open an issue on GitHub.

---

**Last Updated**: October 2025
