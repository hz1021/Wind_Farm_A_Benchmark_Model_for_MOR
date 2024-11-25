# Wind Farm: A Benchmark Model for Model Order Reduction

This project provides the model of a benchmark wind farm model as is detailed in the MATHMOD 2025 conference paper "A Benchmark Model for Model Order Reduction: Large-scale Wind Farms" and the numerical examples therein.


## Prerequisites

The entire wind farm model is modelled using Matlab/Simulink, and can be run on MATLAB R2024a or higher versions.


## Getting Started

This provides both nonlinear (original) model and linear (linearized) model of the wind farm benchmark model.

### Nonlinear model

The nonlinear model is built up using Simulink, see `Non_WindFarm.slx`, and the parameter values used in the Simulink model are given by `Non_ParameterInitialize.m`. For the model of a single wind turbine, please check in `Non_WindFarm.slx` following the path `Non_WindFarm/FOM/Wind Farm/Wind Turbine Strings/String1/WT1/WFG`.

### Linear model

The linear model is represented in the state-space form and the system matrices are stored in `Lin_WFModel_20WTG.mat`.


## Running Numerical Examples

This section provides instructions on how to obtain the numerical example results shown in the conference paper.

### MOR for nonlinear systems

Herein, we show step by step how to use the data-driven nonlinear model reduction method proposed in [Scarciotti and Astolfi (2017)](https://www.sciencedirect.com/science/article/pii/S0005109817300249) to realize model reduction on the original nonlinear wind farm model. 

After downloading the entire project and saving it in a single folder, we open the folder in MATLAB. Firstly, we add the folder to the search path for future MATLAB sessions. Run in the **Command Window**

```
current_folder = pwd
addpath(current_folder)  
```
Then, we initialize the parameters that will be used for the Simulink model. Run

```
Non_ParameterInitialize
```

After finishing parameter initialization, we open and simulate the wind farm model to collect the time snapshots (the input and output data which will be used to realize model reduction afterwards, see [Scarciotti and Astolfi (2017)](https://www.sciencedirect.com/science/article/pii/S0005109817300249)). Run

```
Non_WindFarm
sim Non_WindFarm
```

**Be patient!** This data collection process would take **~5hrs** and that is why people in power engineering would like to do model reduction. When it is done, we save the collected data

```
save Wt.mat Wt -v7.3
save Yt.mat Yt
```

With the collected data, we can now apply the aforementioned model reduction method to the wind farm model. Run

```
Non_ModelReduction
```

Finally, we validate the obtained reduced order model (ROM) by comparing its behaviour with the behaviour of the original full order model (FOM). Run

```
Non_ROMTest
sim Non_ROMTest
out = ans
Non_Plots
```

This final step would take another ~5hrs as the original full model is simulated at the same time. To simulate the ROM only we can comment out the `FOM` module in `Non_ROMTest.slx`, which should take within **1sec**!

### MOR for linear systems

Herein, we show step by step how to use the two-sided moment matching method proposed in [Ionescu (2015)](https://ieeexplore.ieee.org/abstract/document/7336499) to realize model reduction on the linearized wind farm model. Run

```
Lin_WindFarm
Lin_Plots
```

It is expected to obtain the figure which is the same as that in the conference paper, which is **way simpler** than the nonlinear case!


## Authors

* **Hanqing Zhang** and **Zilong Gong** - [MAC-X Lab](https://giordanoscarciotti.com/mac-x-lab/)


## License

This project is licensed under the MIT License - see `LICENSE` file for details.

## Acknowledgments

The authors would like to thank [Dr. Giordano Scarciotti](https://profiles.imperial.ac.uk/g.scarciotti) and [Dr. Adrià Junyent-Ferré](https://profiles.imperial.ac.uk/adria.junyent-ferre) for their dedicated help with building up and reducing the benchmark wind farm model given in this project.


## Appendix: Table of The Model Parameters

In this last section, we give the table of the model parameters which are mentioned in the conference for complementary needs.

![image](https://github.com/user-attachments/assets/6cf0e5e2-4bcc-498b-b011-b591ff58b049)


