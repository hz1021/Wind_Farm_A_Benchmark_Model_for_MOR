# Wind Farm: A Benchmark Model for Model Order Reduction

This project provides the model of a benchmark wind farm model as is detailed in the MATHMOD 2025 conference paper "A Benchmark Model for Model Order Reduction: Large-scale Wind Farms" and the numerical examples therein.

## Prerequisites

The entire wind farm model is modelled using Matlab/Simulink, and can be run on MATLAB R2024a or higher versions.

## Getting Started

This provides both nonlinear (original) model and linear (linearized) model of the wind farm benchmark model.

### Linear Model

The linear model is represented in the state-space form and the system matrices are stored in `Lin_WFModel_20WTG.mat`.

### Nonlinear Model

The nonlinear model is built up using Simulink, see `Non_WindFarm.slx`, and the parameter values used in the Simulink model are given by `Non_ParameterInitialize.m`.


