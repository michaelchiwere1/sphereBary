# Introduction
Here we have the software necessary to perform barycentric interpolation on the sphere. It uses the formulas derived for the sphere on four standard grid types: equally spaced (EQ), shifted equally spaced (SEQ), Gaussian-Legendre (GL), and HEALPix grids.

# Software Dependencies
Our software is written in MATLAB, and apart from having a MATLAB installation, one has to install Chebfun software in his MATLAB environment. If one does not have a Chebfun installation can get one [here](https://www.chebfun.org/).

# Contents of the software
## 1. Home folder
Inside the home folder, we have sample test files that demonstrate how to use this software to interpolate. Each test file was named depending on the kind of interpolation grid that it utilizes.
## 2. +tensor
In this subfolder, we have functions that perform barycentric interpolation on  tensor product grids, namely:

- Equally spaced grid: we use function ```sphereBaryInterpEQ.m ```.
- Shifted equally spaced grid: we use function ```sphereBaryInterpSEQ.m```.
- Gauss-Legendre grid: we use function ```sphereBaryInterpGL.m```. The interpolation weights are computed using function; ```sphereBaryWeights.m```.
## 3. +HEALPix
This subfolder contains the code necessary to perform interpolation on a HEALPix grid.  The file ``` healBaryInterp.m``` implements the barycentric interpolation algorithm for interpolating data on the sphere using the HEALPix grid. Here we use the same trigonometric weights computed using function ```sphereBaryWeights.m``` available in subfolder __+tensor__.





