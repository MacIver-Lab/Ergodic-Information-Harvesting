# Code and data to reproduce results from "Ergodic information harvesting as a behavioral strategy for complex environments"

## Ergodic Information Harvesting (EIH) Tutorial
To ease the effort of understanding how EIH works, we have included an interactive tutorial in Jupyter as well as a video published for a prior publication ("Ergodic Exploration of Distributed Information", Miller et al., IEEE Trans. Robotics, 2016) paper. You can find the link to them below:
- Interactive Jupyter notebook tutorial [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/MacIver-Lab/Ergodic-Information-Harvesting/master?filepath=Tutorial%2FErgodic_Information_Harvesting_Tutorial.ipynb)
- Video of how EIH works in a two-dimensional tracking case using Fisher Information instead of entropy. [Click here to watch the video on YouTube](https://youtu.be/QZ9fGYmJ0G0)

## Steps to reproduce the results shown in the EIH paper
All of the simulation code is written with Python 3.6 using [Jupyter Notebook](http://jupyter.org/) or simply IPython. All of the figure plotting files are written in MATLAB (2017). The simulation can be run on:
- A local computer, which is very easy to setup but limited by the number of accessible CPU cores.
- Cloud computing virtual servers throug any popular Infrastructure as a Service (IaaS) provider, *e.g.* [Amazon Elastic Compute Cloud](https://aws.amazon.com/ec2/) or [Google Cloud Compute Engine](https://cloud.google.com/compute/). CLoud computing is easy to setup and provides a way to scale up the total number of running threads (*e.g.* Google Cloud Compute Engine allows up to 96 CPU threads per instance).

### Setting up the runtime environment
The results generation system uses Jupyter as the environment for interaction. Jupyter is very easy to setup if you have not used it before.

#### Setup Anaconda 
[Anaconda](https://www.anaconda.com/download/) is the required Python environment. It runs on MacOS, Unix, and Windows environments and provide easy package management for setting up the runtime environments for the simulation. To download, go to [https://www.anaconda.com/download/](https://www.anaconda.com/download/) and install Anaconda.

Once installed, open the Anaconda Prompt or bash and run the following code to install all of the dependencies.
The required packages are `scipy` and `numpy`. Note that a specific version has been specified to avoid issues due to any future update to these packages that may cause a backward incompatibility with our code. Most of the other dependencies, such as `jupyter` and `IPython`, has been satisfied by the Anaconda bundle already and are thus not included in the code below.
```bash
conda install scipy=1.0.1 numpy=1.14.2
```
Alternatively, if you prefer `pip`, you could use the code below instead.
```bash
pip install scipy==1.0.1 numpy==1.14.2
```

#### Launch Jupyter
Once the environment is setup, you can launch the simulation through a `jupyter notebook` client. Here's how to do it.
- Navigate to the simulation code folder `/SimulationCode/`
- Launch a new Jupyter notebook instance
  - If using Windows, open a command prompt under this code folder and run `jupyter notebook`
  - If using Linux or MacOS, simply run this command in a terminal: `jupyter notebook` 
- A new browser tab should pop up. Find and open `Ergodic-Information-Harvesting-Simulation.ipynb` under the file list. Note that jupyter notebook will use the current directory of your command prompt as the working directory, if you don't see the code folder, that means your command window is on a different folder and you need to navigate to `/SimulationCode/` before running `jupyter notebook`
- Once opened, you should see the code for simulation and you are now good to go to reproduce the results of our study.

### How to Reproduce Figure Results
To reproduce figure results, first follow the simulation section below to run EIH simulation trials. Then proceed to the figure plotting code to reproduce the figure results. Note that previously simulated data has been included in the step 2 folder so step 1 is optional.

#### Step 1 - Local EIH Simulation (optional)
##### Code Structure
The simulation code files are all stored under `./SimulationCode/` and are organized in a centralized fashion. `Ergodic-Information-Harvesting-Simulation.ipynb` is the only notebook file you need to run and you can use it to reproduce the raw simulation data used for all of the figures in the paper and Supplementary Information. Each dataset is organized in a per-figure fashion and the parameter supporting the simulation is stored under `/SimulationCode/FigParameters/` folder as `json` format.

##### Example Procedure of Simulating Data for figure 1
To simulate the raw data for a given figure, figure 1 for instance, you just need to specify `targetFigure = 'fig1'` in the `Ergodic-Information-Harvesting-Simulation.ipynb` which is further documented inside the notebook. Then run the entire notebook by selecting on the top menu `Cell -> Run All`. Depends on the size of the simulation and the number of thread you allow the program to use, it could take some time to finish. 

Once completed (the program will display `All done! EOF at timestamp = ***` at the bottom), you can find the simulated data under `/SimulationCode/SimData/fig1/` (`/fig1` because we are simulating for figure 1 here). The raw data will be in `*.mat` MATLAB data format that can be loaded by MATLAB.

#### Step 2 - Reproduce Figure Results
##### Code Structure
The figure code files are under `./Production-Figure-Code/`. Similarly, it is centralized in a single MATLAB file `makeFigurePanels.m`. Code for each figure panels are included under `./Production-Figure-Code/FigureCode/` alongside with the previously simulated data if step 1 is skipped.

##### Example Procedure of Reproducing figure 1
To reproduce figure 1 for example, use the following procedures:
- Launch `./Production-Figure-Code/makeFigurePanels.m` using MATLAB. Note that the code has been tested with MATLAB `R2017a` and `R2018a` so any version in between should work. Previous versions after `R2016a` should also work but not tested.
- Specify input parameters
  - Set `targetFig = 'fig1'` to select figure 1 as the target
  - Set `USE_PREV_DATASET = 1` to use previously simulated dataset. Alternatively, use `USE_PREV_DATASET = 0` if simulation data for figure 1 has been created through step 1
- Run the MATLAB code

You should see new MATLAB figure panel pops up contains figure 1 panels and PDF(s) will be saved under `./FigureOutput/fig1/`.
