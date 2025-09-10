# integralARB

[![R-CMD-check](https://github.com/tranbaokhue/integralARB/workflows/R-CMD-check/badge.svg)](https://github.com/tranbaokhue/integralARB/actions)

Rigorous numerical integration for mathematical functions using ball arithmetic with mathematically guaranteed error bounds. Based on the FLINT/ARB library for arbitrary-precision computation.

## Table of Contents

- [Features](#features)
- [System Requirements](#system-requirements)
  - [Linux (Ubuntu/Debian)](#linux-ubuntudebian)
  - [Linux (CentOS/RHEL/Fedora)](#linux-centosrhelfedora)
  - [macOS](#macos)
  - [Windows](#windows)
- [R Package Installation](#r-package-installation)
- [Quick Start](#quick-start)
- [Supported Mathematical Expressions](#supported-mathematical-expressions)
  - [Functions](#functions)
  - [Operations](#operations)
  - [Constants](#constants)
  - [Chain Rule Support](#chain-rule-support)
- [Examples](#examples)
- [Precision and Error Bounds](#precision-and-error-bounds)
- [Troubleshooting](#troubleshooting)
  - [Common Issues](#common-issues)
  - [Getting Help](#getting-help)
- [License](#license)
- [Citation](#citation)
- [References](#references)

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

# Load the package
library(integralARB)
```

## Quick Start

```r
library(integralARB)

# Basic integration: ∫₀^(π/2) cos(x) dx = 1
result <- integrate_rigorous("cos(x)", "0", "pi/2")
print(result)

# High precision integration
result <- integrate_rigorous("sin(x)", "0", "pi", precision = 256)
print(result)

# Complex expressions as limits: ∫₁^(ln(5)) 1/x dx = ln(ln(5))
result <- integrate_rigorous("1/x", "1", "ln(5)")
print(result)
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
First derivatives are supported for simple compositions:
- `sin(2*x)`, `cos(3*x)`, `exp(-x)`
- `log(2*x)`, `atan(5*x)`, `sqrt(3*x)`

## Examples

```r
# Trigonometric integrals
integrate_rigorous("sin(2*x)", "0", "pi")        # = 0
integrate_rigorous("cos(x)", "0", "pi/2")        # = 1

# Polynomial integrals  
integrate_rigorous("x^2", "0", "1")              # = 1/3
integrate_rigorous("x^3", "0", "2")              # = 4

# Exponential integrals
integrate_rigorous("exp(-x)", "0", "1")          # = 1 - 1/e
integrate_rigorous("exp(2*x)", "0", "ln(2)/2")   # = 1/2

# Complex limits
integrate_rigorous("1/(1+x^2)", "0", "sqrt(3)")  # = π/3
integrate_rigorous("sqrt(x)", "0", "pi^2/4")     # = (2/3)*(π/2)^(3/2)
```

## Precision and Error Bounds

The package provides mathematically rigorous error bounds:

```r
# Specify working precision (64-1024 bits)
result <- integrate_rigorous("cos(x)", "0", "pi/2", precision = 128)

# Access high-precision string representation
cat("Value:", result$value_str, "\n")
cat("Error bound: ±", result$error_bound, "\n")

# Rigorous accuracy guarantee
if (result$status == 0) {
  error_num <- as.numeric(result$error_bound)
  accuracy_digits <- -log10(error_num)
  cat("Guaranteed accuracy: >=", floor(accuracy_digits), "decimal digits\n")
}
```

## Troubleshooting

### Common Issues

1. **FLINT not found**: Ensure FLINT is installed and in your system's library path
2. **Compilation errors**: Check that development tools are installed (Xcode on macOS, build-essential on Linux, Rtools on Windows)
3. **Windows DLL issues**: Make sure MSYS2 paths are in your system PATH

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

- Johansson, F. (2017). Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval Arithmetic. IEEE Transactions on Computers.
- FLINT Development Team. FLINT: Fast Library for Number Theory. https://flintlib.org/
