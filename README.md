# Introduction 
Here we have software necessary to compute barycentric interpolation necessary to perform barycentric interpolation on the sphere. It use the formulas derived for the sphere on four types of grid: equally spaced (EQG), shifted equally spaced, Gaussian-Legegendre and HEALPix.encies

# Software Dependencies
Our software is written in MATLAB, and apart from having a MATLAB one has to install Chebfun sotware in his MATLAB environment. If one does not have a Chebfun installation can get one [here](https://www.chebfun.org/).

# Contents of the sotware
## 1. Home folder
Inside the home folder we have sample test files on how to use this software to interpolate depending on the kind of interpolation grid one is using.
## 2. +tensor
In this folder we have functions that perform barycentric interpolation on a tensor product grids, namely:

- Equally spaced grid: we use function ```sphereBaryInterpEQG.m ```.
- Shifted equally spaced grid: we use function ```sphereBaryInterpSEQG.m```.
- Gauss-Legendre grid: we use function ```sphereBaryInterpGL.m```. The interpolation weights are computed using function; ```sphereBaryWeights.m```.
## 3. +HEALPix
This folder contains colder necessary to perform interpolation on a HEALPix grid. The file ```FloaterHormannWght.m``` is the code that implements an algorithm for computing Floater-Hormann weights with cosine, and ```FloaterHormannWghtNcos.m``` evaluates Floater-Hormann weights with cosine. The file ``` healBaryInterp.m``` implements the barycentric interpolation algorithm for interpolating on the sphere using HEALPix grid.
