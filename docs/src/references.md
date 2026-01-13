# References

## Primary Sources

### Rothermel (1972)

**Rothermel, R.C.** 1972. *A mathematical model for predicting fire spread in wildland fuels.* USDA Forest Service Research Paper INT-115. Intermountain Forest and Range Experiment Station, Ogden, Utah. 40 p.

**Key Contributions:**
- Original formulation of the fire spread model
- Derivation of heat source and heat sink equations
- Reaction intensity calculations
- Optimum packing ratio theory
- Moisture damping coefficient formulation

**Citation:**
```bibtex
@techreport{rothermel1972,
  author = {Rothermel, Richard C.},
  title = {A mathematical model for predicting fire spread in wildland fuels},
  institution = {USDA Forest Service, Intermountain Forest and Range Experiment Station},
  year = {1972},
  number = {INT-115},
  type = {Research Paper},
  address = {Ogden, Utah},
  pages = {40}
}
```

**Availability:** [USFS Treesearch](https://www.fs.usda.gov/treesearch/pubs/32533)

---

### Albini (1976a)

**Albini, F.A.** 1976. *Estimating wildfire behavior and effects.* USDA Forest Service General Technical Report INT-30. Intermountain Forest and Range Experiment Station, Ogden, Utah. 92 p.

**Key Contributions:**
- Wind factor formulation and empirical constants
- Slope factor derivation
- Propagating flux ratio equations
- Effective heating number
- Extensions for heterogeneous fuel beds

**Citation:**
```bibtex
@techreport{albini1976,
  author = {Albini, Frank A.},
  title = {Estimating wildfire behavior and effects},
  institution = {USDA Forest Service, Intermountain Forest and Range Experiment Station},
  year = {1976},
  number = {INT-30},
  type = {General Technical Report},
  address = {Ogden, Utah},
  pages = {92}
}
```

**Availability:** [USFS Treesearch](https://www.fs.usda.gov/treesearch/pubs/29574)

---

## Comprehensive Review

### Andrews (2018)

**Andrews, P.L.** 2018. *The Rothermel surface fire spread model and associated developments: A comprehensive explanation.* USDA Forest Service General Technical Report RMRS-GTR-371. Rocky Mountain Research Station, Fort Collins, Colorado. 121 p.

**Key Contributions:**
- Modern comprehensive explanation of all equations
- Historical context and model development
- Detailed derivations and assumptions
- Applications and limitations
- Connections to BEHAVE and other implementations

**Citation:**
```bibtex
@techreport{andrews2018,
  author = {Andrews, Patricia L.},
  title = {The Rothermel surface fire spread model and associated developments: A comprehensive explanation},
  institution = {USDA Forest Service, Rocky Mountain Research Station},
  year = {2018},
  number = {RMRS-GTR-371},
  type = {General Technical Report},
  address = {Fort Collins, Colorado},
  pages = {121},
  doi = {10.2737/RMRS-GTR-371}
}
```

**Availability:** [USFS Treesearch](https://www.fs.usda.gov/treesearch/pubs/55928)

**Note:** This is the most comprehensive modern reference and was the primary source for this implementation.

---

## Fuel Models

### Anderson (1982)

**Anderson, H.E.** 1982. *Aids to determining fuel models for estimating fire behavior.* USDA Forest Service General Technical Report INT-122. Intermountain Forest and Range Experiment Station, Ogden, Utah. 22 p.

**Key Contributions:**
- 13 standard NFFL fuel models
- Fuel parameter values for each model
- Visual guides for fuel model selection
- Application guidelines

**Citation:**
```bibtex
@techreport{anderson1982,
  author = {Anderson, Hal E.},
  title = {Aids to determining fuel models for estimating fire behavior},
  institution = {USDA Forest Service, Intermountain Forest and Range Experiment Station},
  year = {1982},
  number = {INT-122},
  type = {General Technical Report},
  address = {Ogden, Utah},
  pages = {22}
}
```

**Availability:** [USFS Treesearch](https://www.fs.usda.gov/treesearch/pubs/6447)

---

### Scott and Burgan (2005)

**Scott, J.H.; Burgan, R.E.** 2005. *Standard fire behavior fuel models: a comprehensive set for use with Rothermel's surface fire spread model.* USDA Forest Service General Technical Report RMRS-GTR-153. Rocky Mountain Research Station, Fort Collins, Colorado. 72 p.

**Key Contributions:**
- Expanded to 40 standard fuel models
- Improved fuel model descriptions
- Better regional applicability
- Integration with LANDFIRE

**Citation:**
```bibtex
@techreport{scott2005,
  author = {Scott, Joe H. and Burgan, Robert E.},
  title = {Standard fire behavior fuel models: a comprehensive set for use with Rothermel's surface fire spread model},
  institution = {USDA Forest Service, Rocky Mountain Research Station},
  year = {2005},
  number = {RMRS-GTR-153},
  type = {General Technical Report},
  address = {Fort Collins, Colorado},
  pages = {72}
}
```

**Availability:** [USFS Treesearch](https://www.fs.usda.gov/treesearch/pubs/9521)

---

## Model Validation and Testing

### Rothermel and Wilson (1980)

**Rothermel, R.C.; Wilson, R.A.** 1980. *Modeling moisture content of fine dead wildland fuels: input to the BEHAVE fire prediction system.* USDA Forest Service Research Paper INT-226. Intermountain Forest and Range Experiment Station, Ogden, Utah. 61 p.

**Key Contributions:**
- Fuel moisture prediction
- Timelag fuel moisture models
- Integration with fire behavior prediction

**Citation:**
```bibtex
@techreport{rothermel1980,
  author = {Rothermel, Richard C. and Wilson, Ralph A.},
  title = {Modeling moisture content of fine dead wildland fuels: input to the BEHAVE fire prediction system},
  institution = {USDA Forest Service, Intermountain Forest and Range Experiment Station},
  year = {1980},
  number = {INT-226},
  type = {Research Paper},
  address = {Ogden, Utah},
  pages = {61}
}
```

---

## Software Implementations

### BEHAVE

**Andrews, P.L.; Bevins, C.D.; Seli, R.C.** 2008. *BehavePlus fire modeling system, version 4.0: User's Guide.* USDA Forest Service General Technical Report RMRS-GTR-106WWW Revised. Rocky Mountain Research Station, Fort Collins, Colorado. 116 p.

**Description:** Widely used operational fire behavior prediction system implementing Rothermel equations.

**Website:** [BehavePlus](https://www.frames.gov/behaveplus/home)

---

### FARSITE

**Finney, M.A.** 1998. *FARSITE: Fire Area Simulator—model development and evaluation.* USDA Forest Service Research Paper RMRS-RP-4. Rocky Mountain Research Station, Ogden, Utah. 47 p.

**Description:** Spatially explicit fire growth simulation system using Rothermel for spread rate.

**Website:** [FARSITE](https://www.firelab.org/project/farsite)

---

## Julia Ecosystem

### ModelingToolkit.jl

**Rackauckas, C.; Ma, Y.; Dixit, V.; Guo, X.; et al.** 2021. *ModelingToolkit: A composable graph transformation system for equation-based modeling.*

**Description:** Symbolic modeling framework used as the foundation for WildlandFire.jl

**Website:** [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl)

**Documentation:** [ModelingToolkit Docs](https://docs.sciml.ai/ModelingToolkit/stable/)

**Citation:**
```bibtex
@software{modelingtoolkit2021,
  author = {Rackauckas, Christopher and Ma, Yingbo and Dixit, Vaibhav and Guo, Xing and others},
  title = {ModelingToolkit: A composable graph transformation system for equation-based modeling},
  year = {2021},
  url = {https://github.com/SciML/ModelingToolkit.jl}
}
```

---

### SciML Ecosystem

**Rackauckas, C.; Nie, Q.** 2017. *DifferentialEquations.jl – A performant and feature-rich ecosystem for solving differential equations in Julia.* Journal of Open Research Software, 5(1), p.15. DOI: 10.5334/jors.151

**Description:** Comprehensive differential equation solving ecosystem

**Website:** [SciML](https://sciml.ai/)

**Citation:**
```bibtex
@article{rackauckas2017,
  title = {DifferentialEquations.jl – A performant and feature-rich ecosystem for solving differential equations in Julia},
  author = {Rackauckas, Christopher and Nie, Qing},
  journal = {Journal of Open Research Software},
  volume = {5},
  number = {1},
  pages = {15},
  year = {2017},
  doi = {10.5334/jors.151}
}
```

---

## Related Fire Behavior Research

### Frandsen (1971)

**Frandsen, W.H.** 1971. *Fire spread through porous fuels from the conservation of energy.* Combustion and Flame 16:9-16.

**Description:** Theoretical foundation for energy balance in fire spread

---

### Byram (1959)

**Byram, G.M.** 1959. *Combustion of forest fuels.* In: Davis, K.P., ed. Forest fire: control and use. New York: McGraw-Hill: 61-89.

**Description:** Foundational work on fire intensity and flame length

---

### Thomas (1971)

**Thomas, P.H.** 1971. *Rates of spread of some wind-driven fires.* Forestry 44(2):155-175.

**Description:** Theoretical analysis of wind-driven fire spread

---

## Review Articles

### Sullivan (2009)

**Sullivan, A.L.** 2009. *Wildland surface fire spread modelling, 1990-2007. 1: Physical and quasi-physical models.* International Journal of Wildland Fire 18:349-368. DOI: 10.1071/WF06143

**Description:** Comprehensive review of fire spread models including Rothermel

**Citation:**
```bibtex
@article{sullivan2009,
  title = {Wildland surface fire spread modelling, 1990-2007. 1: Physical and quasi-physical models},
  author = {Sullivan, Andrew L.},
  journal = {International Journal of Wildland Fire},
  volume = {18},
  pages = {349--368},
  year = {2009},
  doi = {10.1071/WF06143}
}
```

---

## Online Resources

### USDA Forest Service Fire Research

**Website:** [https://www.fs.usda.gov/research/fire/](https://www.fs.usda.gov/research/fire/)

**Description:** Central hub for USFS fire research including publications and tools

---

### National Wildfire Coordinating Group (NWCG)

**Website:** [https://www.nwcg.gov/](https://www.nwcg.gov/)

**Description:** Operational fire behavior information and training materials

---

### Treesearch

**Website:** [https://www.fs.usda.gov/treesearch/](https://www.fs.usda.gov/treesearch/)

**Description:** USFS research publications database

---

### FRAMES (Fire Research and Management Exchange System)

**Website:** [https://www.frames.gov/](https://www.frames.gov/)

**Description:** Portal for fire science research, tools, and data

---

## Citing WildlandFire.jl

If you use WildlandFire.jl in your research, please cite:

```bibtex
@software{wildlandfire_jl,
  author = {{EarthSciML authors and contributors}},
  title = {WildlandFire.jl: Wildland Fire Behavior Modeling in Julia},
  year = {2026},
  url = {https://github.com/EarthSciML/WildlandFire.jl},
  version = {0.1.0}
}
```

And please cite the original Rothermel work:

```bibtex
@techreport{rothermel1972,
  author = {Rothermel, Richard C.},
  title = {A mathematical model for predicting fire spread in wildland fuels},
  institution = {USDA Forest Service},
  year = {1972},
  number = {INT-115},
  type = {Research Paper}
}
```

For a comprehensive modern reference, also cite:

```bibtex
@techreport{andrews2018,
  author = {Andrews, Patricia L.},
  title = {The Rothermel surface fire spread model and associated developments: A comprehensive explanation},
  institution = {USDA Forest Service},
  year = {2018},
  number = {RMRS-GTR-371},
  type = {General Technical Report},
  doi = {10.2737/RMRS-GTR-371}
}
```

---

## Historical Context

The Rothermel fire spread model was developed in the early 1970s at the USDA Forest Service Intermountain Fire Sciences Laboratory in Missoula, Montana. It represented a major advancement in quantitative fire behavior prediction and has become the foundation for operational fire management systems worldwide.

The model's success stems from:
1. **Physical basis**: Grounded in conservation of energy principles
2. **Practical inputs**: Uses measurable fuel properties
3. **Operational validation**: Extensively tested in real fire conditions
4. **Software implementation**: Integrated into widely-used tools like BEHAVE and FARSITE
5. **Continued development**: Ongoing refinement and extension by the fire research community

WildlandFire.jl brings this established fire behavior model to the modern Julia scientific computing ecosystem, enabling:
- Integration with contemporary modeling workflows
- Automatic differentiation for optimization and sensitivity analysis
- Seamless coupling with other Julia-based Earth system models
- High-performance computing applications
- Reproducible research with open-source tools

---

## Contact and Community

**Package Repository:** [github.com/EarthSciML/WildlandFire.jl](https://github.com/EarthSciML/WildlandFire.jl)

**Documentation:** [wildlandfire.earthsci.dev](https://wildlandfire.earthsci.dev)

**Issues and Questions:** [GitHub Issues](https://github.com/EarthSciML/WildlandFire.jl/issues)

**EarthSciML Project:** [earthsci.dev](https://earthsci.dev)

---

## Acknowledgments

This implementation is based on the foundational work of Richard C. Rothermel and Frank A. Albini on wildland fire behavior modeling. We acknowledge the USDA Forest Service for decades of fire research that has enabled quantitative fire behavior prediction.

The implementation leverages the SciML ecosystem, particularly ModelingToolkit.jl, developed by Christopher Rackauckas and contributors.
