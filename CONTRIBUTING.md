# Contributing to Portfolio Theory with Settlement Frictions

Thank you for your interest in this project! This repository contains replication code for an academic paper. 

## Reporting Issues

If you find bugs or have suggestions:

### Bug Reports

Please include:
1. **MATLAB version** (run `version` in MATLAB)
2. **Operating system** (Windows/Mac/Linux and version)
3. **Exact error message** (copy full text)
4. **Steps to reproduce** (minimal example)
5. **Expected vs actual behavior**

### Feature Requests

For new features or enhancements:
1. Check existing issues first
2. Explain the use case
3. Describe proposed solution
4. Note any paper references if applicable

## Code Style

If contributing code improvements:

### MATLAB Style Guide

```matlab
% Use descriptive variable names
theta_t  % Good
x        % Bad (unless mathematically standard)

% Comment complex operations
% Compute theta evolution using closed-form solution
theta_t = (1 + ((1-theta_0)/theta_0)*exp(lambda*t))^(-1);

% Use consistent spacing
a = b + c;        % Good
a=b+c;           % Bad

% Function headers
function [output1, output2] = myfunction(input1, input2)
%MYFUNCTION Brief description
%   Detailed description
%
%   Inputs:
%       input1 - description
%       input2 - description
%
%   Outputs:
%       output1 - description
%       output2 - description
```

### Documentation

- Add comments for non-obvious code
- Update README.md if adding features
- Include references to paper sections
- Document all function inputs/outputs

## Testing

Before submitting:

```matlab
% Run main script
cd code
interbankmarket_example

% Verify output
% Check that figures are generated
% Compare results with paper figures
```

## Pull Request Process

1. **Fork** the repository
2. **Create branch** (`git checkout -b feature/description`)
3. **Make changes** following style guide
4. **Test thoroughly**
5. **Commit** with clear messages
6. **Push** to your fork
7. **Open Pull Request** with:
   - Clear description of changes
   - Reference to related issues
   - Test results/verification

### Commit Messages

```
Short (50 chars or less) summary

More detailed explanatory text, if necessary. Wrap at 72 characters.
Explain the problem this commit solves. Focus on why you are making 
this change as opposed to how.

References:
- Fixes #123
- Related to paper Section 4.2
```

## Academic Citation

If you use or extend this code:

```bibtex
@article{bianchi2025portfolio,
  title={Portfolio Theory with Settlement Frictions},
  author={Bianchi, Javier and Bigio, Saki},
  journal={Journal of Economic Theory},
  year={2025}
}
```

## Questions?

- **Paper questions**: Contact authors directly
- **Code questions**: Open an issue on GitHub
- **General discussion**: Use GitHub Discussions

## Code of Conduct

- Be respectful and professional
- Focus on constructive feedback
- Welcome newcomers and help them learn
- Academic integrity: properly cite sources

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for helping improve this code!
