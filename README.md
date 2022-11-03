# Chasing Local Realism with Qiskit


![gato](https://github.com/fequinteros/bell_checker/blob/3cbaa1d6c64576f28bbd1aee671659ff5d70175b/logo-project.jpg)


### Repository description

This repository holds the source code implementation for [The Chasing Local Realism with Qiskit project](https://qiskitfallfest.hypeinnovation.com/servlet/hype/IMT?documentTableId=396317333666442236&userAction=Browse&templateName=&documentId=6a036544ae39543c84ea2ffe63841209), a project realized in the context of the [Quantum Hackathon CIC-IPN 2022 (Mexico)](https://qiskitfallfest.hypeinnovation.com/servlet/hype/IMT?documentTableId=396317333666442202&userAction=Browse&templateName=&documentId=a239a36c6092232735d7fc1e7e52aa03).

### Motivation and intention

Recently, the 2022 Nobel Prize in Physics was awarded to a team that proved the existence of quantum entanglement. It was based on this achievement that we decided to study in more detail the Bell's inequalities.


Since there is a family of inequalities of this type the Qiskit bell_checker package provides the necessary framework to compose quantum programs that evaluate any kind of bell inequalitie for a two level system. 

Our goal is to create a package called ***bell_checker***, to allow users to both create the quantum circuits necesary to characterize a bell inequality and also evaluate if it violate or not the desired bell inequality. The bell_checker package will be a value resource for researchers that want to evaluate different bell inequalities and differet scenarios for it in in a scalable program as well for students that want to understand and play with the bell inequalites.

## Example notebooks

### Where I can learn to use this package?

Some jupyter notebooks with examples using `bell_checker` are contained in the folder `Example_notebooks/`. Remember to have `bell_checker` installed on your system before running the notebooks on your computer! Down in the README are the commands to install the `bell_checker` package.


## How to get this package?

### Installation

The `bell_checker` can be installed via pip by running

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

This project makes use of several well known python libraries; [qiskit](https://qiskit.org/), [numpy](https://numpy.org/), [scipy](https://www.scipy.org/), [tqdm](https://tqdm.github.io/) and [matplotlib](https://matplotlib.org/). We acknowledge their developers and express full compliment to their respective licenses.
