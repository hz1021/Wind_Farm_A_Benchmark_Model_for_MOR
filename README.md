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

### MOR for nonlinear systems

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc



