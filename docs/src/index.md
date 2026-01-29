# WildlandFire.jl

```@meta
CurrentModule = WildlandFire
```

WildlandFire.jl provides wildland fire modeling components for the EarthSciML ecosystem, including implementations of the National Fire Danger Rating System (NFDRS).

## Overview

This package implements equation systems for wildland fire danger assessment and behavior modeling using [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/). The components can be used standalone or composed with other EarthSciML modules for integrated Earth system modeling.

## Features

- **Fuel Moisture Models**: Dead fuel moisture (1-hr, 10-hr, 100-hr, 1000-hr timelag classes) and live fuel moisture (herbaceous and woody)
- **Fire Behavior Indices**: Spread Component, Energy Release Component, Burning Index
- **Fire Occurrence Indices**: Ignition Component, Human-Caused Fire Occurrence Index, Fire Load Index
- **Fuel Model Database**: All 20 NFDRS fuel models (A-U, excluding M) with vegetation-specific parameters

## Quick Start

```julia
using WildlandFire
using ModelingToolkit

# Create an equilibrium moisture content model
emc = EquilibriumMoistureContent()

# Get fuel model parameters
fuel_model_a = get_fuel_model(:A)  # Western grasses (annual)

# Create a burning index model
bi = BurningIndex()
```

## Contents

```@contents
Pages = ["nfdrs.md"]
Depth = 2
```

## Index

```@index
```
