# Chasing Local Realism with Qiskit

### Repository description

This repository holds the source code implementation for [The Chasing Local Realism with Qiskit project](https://qiskitfallfest.hypeinnovation.com/servlet/hype/IMT?documentTableId=396317333666442236&userAction=Browse&templateName=&documentId=6a036544ae39543c84ea2ffe63841209), a project realized in the context of the [Quantum Hackathon CIC-IPN 2022 (Mexico)](https://qiskitfallfest.hypeinnovation.com/servlet/hype/IMT?documentTableId=396317333666442202&userAction=Browse&templateName=&documentId=a239a36c6092232735d7fc1e7e52aa03).

### Motivation and intention

Recently, the 2022 Nobel Prize in Physics was awarded to a team that proved the existence of quantum entanglement. It was based on this achievement that we decided to study in more detail the Bell's inequalities.


Since there is a family of inequalities of this type the Qiskit bell_checker module provides the necessary framework to compose quantum programs that evaluate any kind of bell inequalitie for a two level system. 

Our goal is to create a module called ***bell_checker***, to allow users to both create the quantum circuits necesary to characterize a bell inequality and also evaluate if it violate or not the desired bell inequality. This would be a value resource for students that want to understand and play with some bell inequalites as well for researchers that want to evaluate a lot of different scenarios in a scalable program.

## Example notebooks

### Can I see some examples of this package in action?

Some jupyter notebooks with examples using `bell_checker` are contained in the folder `Example_notebooks/`. Remember to have `bell_checker` installed on your system before running the notebooks!

### Can I download the example notebooks to experiment locally?

Yes! The folder `Example_notebooks/` can be downloaded as a compressed zip file from [this link](https://gitlab.com/fequinteros/bell_checker-/archive/master/bell_checker-master.zip?path=Example_notebooks).


## How to get this package?

### Installation

This package can be installed via pip by running

```sh
pip install git+https://github.com/fequinteros/bell_checker
```
    
another way is running

```sh
python -m pip install git+https://github.com/fequinteros/bell_checker
```
    
### Uninstallation

The `bell_checker` package can be uninstalled as

```sh
pip uninstall bell_checker
```

also works,

```sh
python -m pip uninstall bell_checker
```

# License acknowledgments

This project makes use of several well known python libraries; [qiskit](https://qiskit.org/), [numpy](https://numpy.org/), [scipy](https://www.scipy.org/) and [matplotlib](https://matplotlib.org/). We acknowledge their developers and express full compliment to their respective licenses.
