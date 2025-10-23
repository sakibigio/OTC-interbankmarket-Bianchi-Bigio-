# Quick Start Guide

Get up and running with the code in 5 minutes.

## Prerequisites

- MATLAB R2020b or later
- No additional toolboxes required

## Installation

```bash
# Clone repository
git clone https://github.com/yourusername/interbankmarket-code.git
cd interbankmarket-code
```

## Run Everything

```matlab
% In MATLAB
cd code
interbankmarket_example
```

Wait ~3 minutes. All figures will be saved to `output/figures/`.

## Generate Specific Figures

### Just Leontief Matching

```matlab
% Edit interbankmarket_example.m
% Comment out lines for other matching functions
% Keep only Leontief sections (lines ~226-400)
```

### Just One Figure

```matlab
% Example: Generate only F_l_Chiplus_theta.pdf

% Set parameters
lambda_o = 1;
eta_o = 0.5;
delta_r = 12;  % In BPS
N = 10000;

% Initialize
matchtech.lambda = lambda_o;
matchtech.eta = eta_o;
N_theta = 101;
ln_THETA = linspace(log(0.125), log(1/0.125), N_theta);
CHIp = zeros(N_theta, 1);

% Solve
for jj = 1:N_theta
    theta_0 = exp(ln_THETA(jj));
    [~, CHIp_0] = analytic_leontief(0, theta_0, delta_r, matchtech);
    CHIp(jj) = CHIp_0;
end

% Plot
figure
plot(ln_THETA, CHIp, 'LineWidth', 2)
xlabel('ln(\theta_0)')
ylabel('BPS')
grid on
```

## Modify Parameters

In `interbankmarket_example.m`, lines 76-78:

```matlab
lambda_o = 1;      % Try: 0.5, 2, 5
eta_o = 0.5;       % Try: 0.25, 0.75
theta_bench = 0.75; % Try: 0.5, 1.5
```

Then rerun.

## Common Tasks

### Change Output Directory

```matlab
foldername = '/your/path/here/';
```

### Skip Figure Export (Faster)

```matlab
printit = 0;  % Line 34
```

### Reduce Computation Time

```matlab
N = 1000;        % Instead of 10000 (line 62)
N_theta = 50;    % Instead of 101 (line 96)
```

### Compare to Paper

Paper figures are in:
- Figure 1 → `F_l_Chiplus_theta.pdf`
- Figure 2 → `F_cd_Chiplus_theta.pdf`
- etc.

## Troubleshooting

**Problem**: "Function not found"
```matlab
% Solution
addpath(genpath('code'))
```

**Problem**: "Out of memory"
```matlab
% Solution
N = 1000;  % Reduce rounds
```

**Problem**: Figures look different
```matlab
% Check parameters match paper
lambda_o = 1
eta_o = 0.5
delta_r = 12  % Should be 12 BPS
```

## Next Steps

- Read `README.md` for full documentation
- Check `code/README.md` for technical details
- See paper for theory

## Getting Help

1. Check `README.md`
2. Read code comments
3. Open issue on GitHub

---

Ready to customize? See `README.md` for detailed documentation.
