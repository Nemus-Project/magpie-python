# Paper abstract

## Outline

Max: 200 words

- Description of context
- what do you want to talk about
  - Where you from
  - what you doing
  - method
  - outcomes

- MATLAB source for deriving mode patterns of a thin plate.
  - under abritrary / parametric boundary conditions
- Design an interactive tool for exploring parameters
  -  modal analysis of elastic rectangularly restrained plates
  -  Easily generate numerical mode patterns for comparison against modes obtained against experimental measurements

## Main Abstract

Online their exist tools for deriving modes for physical acoustic system but rarely are they extensive and open source. The Eigen Plate project differs in it's aim to not only be accessible but also transparent and free open source software. The Eigen Plate is software developed as part of the Numerical Restoration of Historical Musical Instruments [NEMUS] project. Part of the NEMUS project remit is developing projects for the digital simulation and conservation of historic musical instruments. Eigen Plate is derived from an in-house, MATLAB source, tool developed for deriving mode shapes and mode frequencies for a parametric isotropic thin plate. Parameters for the plate model are width, length, thickness, density, young's modulus, poisson ratio and elastically constrained boundary conditions. Interactive elements provide a more intuitive approach to parameter exploration.
<!-- summarize width/length/thickness etc... into something like "physical properties", then I guess you could be more specific into the paper -->
Using current, open-source, non-proprietary development platforms Eigen Plate allows for easy access an exploration of parameters of a thin plate mode and their effect on mode shapes and frequencies. 

Obtaining mode shapes for plates can provide useful way for instrument makers to further understand the material with which they are working. Eigen Plate aims to also make conceptualising concepts such as mode patterns to a wider audience and has potential as a research tool for curators of musical instrument collections as well as creative tool for artists.


## Transcript

### The context

Online tools for deriving modes for physical acoustic systems, but the EigenPlate differs in it's aim to be entirely trasparent and a free open source software.

### The plan

Part of the Numerical Restoration of Historical Musical Instruments [NEMUS] project, aimed at develop projects for the digital simulation and conservation of historic musical instruments.
the plan is folks, taht I want to take the biharmEigs script that michele has writeen as a source, as a base for generating an open source tool for eploring the results of the script. The point is to make it easily accessible so that you don't relly need to understand what is going on underneath the hood, but rather you can just give a tool a few parameters and  from there be presented a some mode patterns.  Also, want to make easier to explire, providing a few interactive elements that provid ea more intuitive means to provide parameters to the tool.Interesting thing is the accessibility, not the tool.

### Why is it rlevenat

- Pracycally: 
  - knowing  the mode shapes from an instrument makers perspective can be useful to further understand the material they are working with. Tol to understand the modal behaviour of materiual.
  - as a research tool for curators of musical instrumet collections
  - For create purposes (hunt down literature, possibly something Gadi Sassoon did)
- Making modal analysis more accessible and easier to conceptualise
- and make it open source for easy expansion
- 


#### Questions

- What is `biharmEigs.md`
  - it is a code that give you the mode shapes and mode frequencies for a gievn isotropic thin plate (with width, length, thickness, density, young's modulus, poisson ratio) with parametric elastically constrained boundary conditions.
