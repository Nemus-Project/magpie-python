# m<img src="./img/magpie.svg" style="height:1ch;"/>gpie 


## About

Finite Difference implementation of a stiff plate with adjustable boundary conditions.

Elastic boundary conditions adjustable both in torsional $R$ stiffness and transversal stiffness $K$

## Functionality

- Inputs
  - Material
    - from a data set of coefficients from literature  
      - named
      - or set range manually
  - Continuous scale for elastic boundary conditions
- Outputs
  - Image (or series)
    - downloadable in common formats (jpg, png)
      - modal patterns
      - in 2D
        - with a colour map
        - or with simple node lines
      - in 3D 
        - ThreeJS interface
        - adjustable camera
    - Cached Results
      - slide across scale of elastic conditions and Young's
  - data as text / common data exchange formats (csv / tsv / JSON)
    - input paremeters
      - across a scale Young's moduli
    - derived coefficients
      - Modal frequencies
  - Audio
    - Impulse response from a FDTD simulation of the plate
    - Pure tones of eigen frequencies
  - Plots
    - plot frequency of IR against Eigen Frequencies
    - sonify tones of eigen frequencies

## Goals

- create an interface for user to explore parameters
  - Parameters
    - For each edge
      - 2 springs
        - rotational
        - transversal
    - Grid size
    - Material
    - Geometric Parameters
    - Young's Modulus
    - density
- through the browser

## Dependencies

### MATLAB

- fastsmooth [(Fast smoothing function Version 2.0.0.0 (52.4 KB) by Tom O'Haver)](https://it.mathworks.com/matlabcentral/fileexchange/19998-fast-smoothing-function)
- Symbolic Math Toolbox

## References
