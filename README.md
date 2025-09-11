# integralARB

[![R-CMD-check](https://github.com/tranbaokhue/integralARB/workflows/R-CMD-check/badge.svg)](https://github.com/tranbaokhue/integralARB/actions)

Rigorous numerical integration for mathematical functions using ball arithmetic with mathematically guaranteed error bounds. Based on the FLINT/ARB library for arbitrary-precision computation.

## Table of Contents

- [Features](#features) • [System Requirements](#system-requirements) • [Installation](#r-package-installation) • [Quick Start](#quick-start)
- [Supported Expressions](#supported-mathematical-expressions) • [Examples](#examples) • [Precision](#precision-and-error-bounds)
- [Troubleshooting](#troubleshooting) • [License](#license) • [Citation](#citation) • [References](#references)

## Features

- **Mathematically rigorous results** with guaranteed error bounds
- **High-precision integration** (64-1024 bits of precision)
- **Expression parser** supporting common mathematical functions
- **Flexible integration limits** using mathematical expressions (e.g., `pi/2`, `ln(5)/3`, `sqrt(2)`)
- **Cross-platform compatibility** (Linux, macOS, Windows with Rtools)

## System Requirements

This package requires the FLINT library and its dependencies. Follow the installation instructions for your operating system below.

### Linux (Ubuntu/Debian)

```bash
# Update package list
sudo apt-get update

# Install FLINT and dependencies
sudo apt-get install libflint-dev libgmp-dev libmpfr-dev

# For older Ubuntu versions, you may need:
sudo apt-get install build-essential
```

### Linux (CentOS/RHEL/Fedora)

```bash
# CentOS/RHEL with EPEL
sudo yum install epel-release
sudo yum install flint-devel gmp-devel mpfr-devel

# Fedora
sudo dnf install flint-devel gmp-devel mpfr-devel
```

### macOS

```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install FLINT and dependencies
brew install flint gmp mpfr

# Ensure Xcode command line tools are installed
xcode-select --install
```

### Windows

Windows installation requires Rtools and MSYS2:

1. **Install Rtools**: Download and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) compatible with your R version.

2. **Install FLINT via MSYS2**:
   ```bash
   # Open MSYS2 terminal (comes with Rtools)
   pacman -S mingw-w64-x86_64-flint
   pacman -S mingw-w64-x86_64-gmp
   pacman -S mingw-w64-x86_64-mpfr
   ```

3. **Alternative - vcpkg**:
   ```bash
   # If you prefer vcpkg
   vcpkg install flint:x64-windows
   vcpkg install gmp:x64-windows
   vcpkg install mpfr:x64-windows
   ```

## R Package Installation

Once system dependencies are installed:

```r
# Install from GitHub
if (!require(devtools)) install.packages("devtools")
devtools::install_github("tranbaokhue/integralARB")
```

## Examples

```r
# Load the package
library(integralARB)

result <- integrate_rigorous("2*sin(2*x) + 2*exp(x*3)", "-pi/6", "pi/3")
print(result)

Rigorous Integration Result
===========================
Expression:  2*sin(2*x) + 2*exp(x*3) 
Limits:     [ -pi/6 ,  pi/3 ]
Precision:   128  bits

Value:
   16.2885420376190047314547538320757124056766981590 
Rigorous error bound:
  ± 2.72714689390770701848746876635560997253930409171e-36 
Guaranteed accuracy: ≥ 35  decimal digits
Evaluations: 10000  (estimated)
```

## Supported Mathematical Expressions

### Functions
- Trigonometric: `sin(x)`, `cos(x)`, `atan(x)`
- Exponential/Logarithmic: `exp(x)`, `log(x)`, `ln(x)`
- Hyperbolic: `sinh(x)`, `cosh(x)`
- Other: `sqrt(x)`

### Operations
- Basic arithmetic: `+`, `-`, `*`, `/`, `^`
- Parentheses for grouping: `()`

### Constants
- Mathematical constants: `pi`, `e`
- Complex expressions: `pi/2`, `ln(5)/3`, `sqrt(2)/4`

### Chain Rule Support
First derivatives are supported for simple compositions such as: `sin(2*x)`, `cos(3*x)`, `exp(-x)`, `log(2*x)`, `atan(5*x)`, `sqrt(3*x)`.

## Troubleshooting

### Common Issues

1. **FLINT not found**: Ensure FLINT is installed and in your system's library path. Windows support requires manual FLINT installation. See detailed Windows instructions above.
2. **Compilation errors**: Check that development tools are installed (Xcode on macOS, build-essential on Linux, Rtools on Windows).
3. **Windows DLL issues**: Make sure MSYS2 paths are in your system PATH.
4. **Max depth**: Max depth reached. Make sure you increase the Maximum subdivision depth with increased precision.
5. **Failed to parse**: The integrand needs to follow the format for it to be parsed correctly and send to ARB precisely. Moreover, the expression you want might not be supported yet, so please reach out to Khue to inquire about addiung your specific function as needed.

### Getting Help

- Check supported expressions: `list_supported_expressions()`
- Report bugs: [GitHub Issues](https://github.com/tranbaokhue/integralARB/issues)
- For FLINT installation issues, consult the [FLINT documentation](https://flintlib.org/)

## License

This package is licensed under LGPL-3. It uses the FLINT library, which is also under LGPL-2.1+.

## Citation

Tran, B.K. (2024). integralARB: Rigorous Numerical Integration using ARB Library. 
R package version 0.1.0. https://github.com/tranbaokhue/integralARB

## References

[1] F. Johansson, “Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval Arithmetic,” IEEE Transactions on Computers, vol. 66, no. 8, pp. 1281–1292, Aug. 2017, doi: 10.1109/TC.2017.2690633.

[2] FLINT Development Team. FLINT: Fast Library for Number Theory. https://flintlib.org/
