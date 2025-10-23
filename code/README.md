# Code Documentation

This directory contains all MATLAB code for "Portfolio Theory with Settlement Frictions".

## Directory Structure

```
code/
├── interbankmarket_example.m     # MAIN SCRIPT - Run this file
├── solvers/                       # Equilibrium solution algorithms
├── analytics/                     # Closed-form analytical solutions
├── utilities/                     # Helper functions
└── plotting/                      # Figure formatting utilities
```

## Running the Code

### Basic Usage

```matlab
% Navigate to code directory
cd /path/to/interbankmarket-code/code

% Run main script
interbankmarket_example
```

### Advanced Usage

```matlab
% Add all subdirectories to path
addpath(genpath(pwd))

% Run with custom parameters
lambda_o = 2.0;  % Higher matching efficiency
eta_o = 0.75;    % Higher borrower bargaining power
interbankmarket_example
```

## Module Descriptions

### Main Script: `interbankmarket_example.m`

The master file that:
1. Sets global parameters and figure formatting
2. Compares different matching functions
3. Generates comparative statics plots
4. Exports figures to PDF

**Key Sections**:
- Lines 57-108: Parameter definitions
- Lines 110-138: Matching function growth rate comparison
- Lines 140-183: Liquidity yield comparisons across matching functions
- Lines 185-224: Theta evolution trajectories
- Lines 226-400: Leontief matching analysis
- Lines 401-800: Cobb-Douglas matching analysis  
- Lines 801-1100: Harmonic mean matching analysis
- Lines 1426-1465: Walrasian limit analysis
- Lines 1467-1588: Additional harmonic mean plots

### Solvers

#### `interbanksolve_discrete.m`
Solves the model with discrete trading rounds using backward induction.

**Inputs**:
- `delta_r`: Interest rate spread (r^w - r^m)
- `theta_o`: Initial market tightness
- `N`: Number of trading rounds
- `matchtech`: Structure with matching parameters

**Outputs**:
- `xi_d`, `xi_s`: Initial liquidity yields (deficit/surplus)
- `r_i`: OTC rates vector
- `gammad`, `gammas`: Matching probability vectors
- `theta_t`: Market tightness evolution
- `Vd`, `Vs`: Value functions

**Algorithm**:
1. Forward iteration: Compute theta evolution and matching probabilities
2. Backward iteration: Solve value functions via dynamic programming
3. Returns equilibrium at t=0

#### `interbanksolve_continuous.m`
Solves the continuous-time limit as Δ → 0.

Similar structure to discrete solver but uses:
- Differential equations for theta evolution
- Continuous-time matching rates
- Integration over [0,1] interval

#### `interbanksolve_analyticleontief.m`
Uses closed-form solutions for Leontief matching (perfect complements).

**Key Formula**:
For θ₀ < 1: θₜ = (1 + ((1-θ₀)/θ₀)exp(λt))⁻¹

For θ₀ > 1: θₜ = 1 + (θ₀-1)exp(λt)

#### `interbanksolve_analyticcobbdouglas.m`
Uses closed-form solutions for Cobb-Douglas matching.

**Key Formula**:
θₜ = ((α·exp(-λt) - 1)/(α·exp(-λt) + 1))²
where α = (1 + √θ₀)/(1 - √θ₀)

**Special Feature**: Handles finite stopping time T* when market clears

#### `interbanksolve_analyticharmonic.m`
Uses closed-form solutions for harmonic mean matching (p = -1).

**Key Feature**: Handles both stable (θ₀ < 1) and unstable (θ₀ > 1) branches of solution

### Analytics

These functions compute liquidity yields and OTC rates for specific matching technologies.

#### `analytic_leontief.m`
Returns yields and rates for Leontief matching at time t.

**Mathematical Implementation**:
```matlab
% For theta_0 < 1
theta_t = (1 + ((1-theta_0)/theta_0)*exp(lambda*t))^(-1);
CHIp_t = (delta_r)*(theta_1/theta_t)^eta * ...
         ((theta_t - theta_t^eta*theta_1^(1-eta))/(1-theta_1));
```

#### `analytic_cobbdouglas.m`  
Handles all three regimes:
- θ₀ < 1: Surplus market
- θ₀ = 1: Balanced market (special limit)
- θ₀ > 1: Deficit market

**Stopping Time**: Computes T* = (1/λ)log|α| where market exhausts

#### `analytic_harmonic.m`
Most complex analytical solution due to two-branch structure.

**Solution Branches**:
- Stable branch (θ₀ < 1): θₜ = 2a/(1 + 2a + √(1+4a))
- Unstable branch (θ₀ > 1): θₜ = (1 + 2a + √(1+4a))/(2a)

where a = A·exp(-λt) and A = θ₀/(1-θ₀)²

### Utilities

#### `build_matching.m`
Constructs CES matching function with elasticity ρ.

**Function Signature**:
```matlab
function [G] = build_matching(rho, lambda)
```

**Returns**:
Anonymous function G(a,b) implementing:
- ρ = 1: Cobb-Douglas (geometric mean)
- ρ = 0: Leontief (minimum)
- ρ ∈ (0,1): Intermediate CES

**Formula**:
```matlab
G = @(a,b) lambda*(1/2)^(rho/(rho-1))* ...
           (a.^(1-1/rho) + b.^(1-1/rho)).^(rho/(rho-1))
```

#### `interbank_resetmatchtech.m`
Resets matching technology structure to baseline values.

**Usage Pattern**:
```matlab
interbank_resetmatchtech;  % Sets matchtech.lambda and matchtech.eta
```

#### `Plot_Migration.m`
Automatically copies generated figures to Overleaf directory.

**Configuration**:
```matlab
src = 'Figures/';
dst = '/path/to/Overleaf/Figures/';
```

### Plotting Utilities

These add formatting elements to figures (corridors, reference lines, labels).

#### Common Pattern

```matlab
% After creating a plot
plot(x, y);
interbankplot_addcorridortheta;  % Adds r^w-r^m corridor and formatting
```

All plotting utilities:
- Add horizontal lines at 0 and δr = r^w - r^m
- Add vertical line at θ = 1 (balanced market)
- Format axes with appropriate tick marks
- Add LaTeX-formatted labels

## Data Flow

```
Main Script
    │
    ├─> Set Parameters (λ, η, θ₀, N)
    │
    ├─> Build Matching Function (build_matching.m)
    │
    ├─> Solve Equilibrium
    │   ├─> Analytical (analytic_*.m)
    │   └─> Numerical (interbanksolve_*.m)
    │
    ├─> Generate Figure
    │   ├─> Plot data
    │   └─> Add formatting (interbankplot_*.m)
    │
    └─> Export PDF (exportgraphics)
```

## Numerical Considerations

### Precision and Accuracy

- All calculations in double precision (~10⁻¹⁵ relative error)
- Analytical solutions are exact (no discretization error)
- Numerical solvers: O(1/N) convergence rate

### Stability

**Potential Issues**:
1. **Extreme theta values**: For θ → 0 or θ → ∞
   - Solution: Use log-space for theta grid
   - Implemented: `ln_THETA = linspace(-5, 5, N_theta)`

2. **Numerical overflow**: In exponential terms
   - Solution: Check for exp overflow before computation
   - Implemented: Tests like `if theta_0 < 1` branch logic

3. **Division by zero**: When θ₁ → 1
   - Solution: L'Hôpital's rule applied
   - Implemented: Special case `if theta_0 == 1`

### Performance Optimization

**Current Performance** (N=10,000):
- Full script: ~3 minutes
- Single solve: ~0.1 seconds
- Memory: <500 MB

**Bottlenecks**:
1. Figure generation (45% of runtime)
2. Discrete solver backward iteration (30%)
3. Grid computations for theta plots (25%)

**Optimization Tips**:
```matlab
% Reduce figure count
printit = 0;  % Skip PDF export

% Reduce resolution
N_theta = 50;  % Instead of 101

% Reduce rounds
N = 1000;  % Instead of 10,000
```

## Matching Function Theory

### CES Specification

General form: G(a,b) = λ((aᵖ + bᵖ)/2)^(1/p)

**Special Cases**:
- p → -∞: min(a,b) [Leontief]
- p = -1: 2ab/(a+b) [Harmonic]
- p = 0: √(ab) [Cobb-Douglas]
- p → 1: (a+b)/2 [Linear]

**Elasticity of Substitution**:
σ = 1/(1-p) = ρ

### Matching Rates

Given tightness θ = S⁻/S⁺:

**Surplus side**: γˢ(θ) = G(θ,1)
**Deficit side**: γᵈ(θ) = G(1/θ,1)

**Properties**:
- γˢ increasing in θ
- γᵈ decreasing in θ  
- γˢ(1) = γᵈ(1) = λ/2

## Testing and Validation

### Unit Tests

To verify analytical solutions match numerical:

```matlab
% Set parameters
theta_0 = 0.5;
matchtech.lambda = 1;
matchtech.eta = 0.5;
N = 10000;

% Analytical
[CHIm_a, CHIp_a] = analytic_leontief(0, theta_0, delta_r, matchtech);

% Numerical  
matchtech.rho = 0;
[CHIm_n, CHIp_n] = interbanksolve_continuous(delta_r, theta_0, N, matchtech);

% Compare
assert(abs(CHIm_a - CHIm_n) < 1e-3);
assert(abs(CHIp_a - CHIp_n) < 1e-3);
```

### Consistency Checks

Built into main script:
1. Symmetry: CHI⁺(θ) = CHI⁻(1/θ) when η = 0.5
2. Bounds: 0 ≤ CHI± ≤ δr
3. Convergence: Analytical vs numerical agreement
4. Balance: CHI⁺(1) = CHI⁻(1) when θ₀ = 1

## Extension Guide

### Adding New Matching Function

1. Create `analytic_newfunction.m` in `analytics/`
2. Implement with signature:
   ```matlab
   function [CHIm_t, CHIp_t, r_f_t, gammad_t, gammas_t, theta_t] = ...
       analytic_newfunction(t, theta_0, delta_r, matchtech)
   ```

3. Create solver `interbanksolve_analyticnewfunction.m`
4. Add to main script comparison sections

### Adding New Figure

1. Create plotting block in main script
2. Set up grid and initialize arrays
3. Loop over parameters calling solvers
4. Plot with formatting:
   ```matlab
   figure('Name', 'Descriptive Title')
   plot(x, y)
   interbankplot_addcorridor[type]
   % Additional formatting
   if printit==1
       exportgraphics(gca, [foldername 'F_newplot.pdf'])
   end
   ```

## Common Pitfalls

1. **Forgetting to reset matchtech**: Always call `interbank_resetmatchtech` before new analysis section

2. **Wrong path for figures**: Ensure `foldername` exists or script will fail silently

3. **Mixing analytical/numerical solvers**: They use different time grids (continuous vs discrete)

4. **Not scaling rates**: Parameters are annual, but plots often show BPS (basis points)

## References

For mathematical details, see paper:
- Section 2: Model setup
- Section 3: Analytical solutions (Leontief, Cobb-Douglas)
- Section 4: Numerical methods
- Appendix A: Proofs of analytical formulas

---

**Maintained by**: Saki Bigio  
**Last Updated**: April 2025
