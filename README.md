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

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

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
