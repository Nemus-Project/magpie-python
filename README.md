# m<img src="./img/magpie.svg" style="height:1ch;"/>gpie 


## About

Magpie is an open source framework for working with thin plates with generalized elastic boundary conditions. The framework is designed around two key use cases

1. Exploration of plates for understanding how material properties affect mode shapes and modal frequencies.
2. Derivation of material properties from experimental data.

## Structure

All source code is found in the `src/` directory of the repository. In `src/` you can find a directory for each supported language as well as a `data/` directory which contains any datasets shared across implementations.

### Magpie function

The magpie function accepts:

- plate material properties: density, elasticity (Young's modulus), poisson ratio
- plate constraints: dimensions, boundary conditions

The user can stipulate an accuracy coefficient and how many modes they wish to calculate.

The output from `magpie` is:

- `Q`: A list of eigen vectors, one for each mode. This can be used to visualize the mode shape.
- `Om`: The angular frequency for the corresponding eigen vector.
- `N`:  Number of grid points in the $x$ and $y$ dimensions for the plate. A smaller `h` value will result in a larger number of grid points
- `biharm`: The biharmonic used for deriving eigen vectors and modal frequencies. The biharmonic can also be used for a finite difference difference time domain scheme of the plate.

## References
