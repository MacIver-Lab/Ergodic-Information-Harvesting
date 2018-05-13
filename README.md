# Code and data to reproduce results from "Ergodic information harvesting as a behavioral strategy for complex environments" by Chen Chen, Todd D. Murphey, and Malcolm A. MacIver, Northwestern University, Evanston IL, USA

## Ergodic Information Harvesting (EIH) Tutorial
To ease the effort of understanding how EIH works, we have included an interactive tutorial in Jupyter as well as a video published for a prior publication (["Ergodic Exploration of Distributed Information"](https://nxr.northwestern.edu/sites/default/files/publications/Mill16a_ergodic_control_distributed_info.pdf), Miller et al., IEEE Trans. Robotics, 2016) paper. 
- Interactive Jupyter notebook tutorial [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/MacIver-Lab/Ergodic-Information-Harvesting/master?filepath=Tutorial%2FErgodic_Information_Harvesting_Tutorial.ipynb)
- Video of how EIH works (Two differences from study: 1. Object to be found is not moving; 2. Uses Fisher Information instead of entropy.) [Click here to watch the video on YouTube](https://youtu.be/QZ9fGYmJ0G0)

## Steps to reproduce the results shown in the EIH paper
All of the simulation code is written with Python 3.6 using [Jupyter Notebook](http://jupyter.org/). All of the figure plotting files are written in MATLAB (2017). The code can be run on:
- A local computer, which is very easy to set up but limited by the number of accessible CPU cores.
- Cloud computing virtual servers throug any popular Infrastructure as a Service (IaaS) provider, *e.g.* [Amazon Elastic Compute Cloud](https://aws.amazon.com/ec2/) or [Google Cloud Compute Engine](https://cloud.google.com/compute/). Cloud computing is easy to setup and provides a way to scale up the total number of running threads (*e.g.* Google Cloud Compute Engine allows up to 96 CPU threads per instance).

### Setting up the runtime environment
The results generation system uses Jupyter as the environment for interaction. Jupyter is very easy to set up if you have not used it before.

#### Install Anaconda 
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
- A new browser tab should pop up. Find and open `Ergodic-Information-Harvesting-Simulation.ipynb` under the file list. Note that Jupyter notebook will use the current directory of your command prompt as the working directory: if you don't see the code folder, that means your command window is in a different folder and you need to navigate to `/SimulationCode/` before running `jupyter notebook`
- Once opened, you should see the code for simulation and you are now good to go to reproduce the results of our study.
- **NOTE**: Depending on your operating system, you **may need to prevent your system from going to sleep**. This is necessary with MacOS. With MacOS, there are two ways to do this: 1. Open a terminal, and type `caffeinate` and hit return. Your system will be prevented from sleeping until you hit Control-C. 2. Go to System Preferences, Energy Saver panel, and click the box that says "Prevent computer from sleeping automatically when the display is off".

### How to Reproduce Figure Results
To reproduce the figure results, first follow the simulation section below to run EIH simulation trials. Then proceed to the figure plotting code to reproduce the figure results. As it can take several days to run the hundreds of simulations needed for the computational results we show in the paper, these data have been included and can be viewed by setting a flag in the figure generation code in Step 2 below. This also provides a baseline of comparison for results that you generate.

#### Step 1 - Local EIH Simulation (optional)
##### Code Structure
The simulation code files are all stored under `./SimulationCode/` and are organized in a centralized fashion. `Ergodic-Information-Harvesting-Simulation.ipynb` is the only notebook file you need to run and you can use it to reproduce the raw simulation data used for all of the figures in the paper and Supplementary Information. Each dataset is organized in a per-figure fashion and the parameter supporting the simulation is stored under `/SimulationCode/FigParameters/` folder in `json` format.

##### Example Procedure of Simulating Data for figure 1
To simulate the raw data for a given figure, figure 1 for instance, you just need to specify `targetFigure = 'fig1'` in the `Ergodic-Information-Harvesting-Simulation.ipynb` which is further documented inside the notebook. Then run the entire notebook by selecting on the top menu `Cell -> Run All`. Depending on the size of the simulation and the number of threads you allow the program to use, it could take some time to finish. 

Once completed (the program will display `All done! EOF at timestamp = ***` at the bottom), you can find the simulated data under `/SimulationCode/SimData/fig1/` (`/fig1` because we are simulating for figure 1 here). The raw data will be in `*.mat` MATLAB data format that can be loaded by MATLAB.

#### Step 2 - Reproduce Figure Results
##### Code Structure
The figure code files are under `./Production-Figure-Code/`. Similarly, it is centralized in a single MATLAB file `makeFigurePanels.m`. Code for each figure panels are included under `./Production-Figure-Code/FigureCode/` along with the previously simulated data if step 1 is skipped.

##### Example Procedure of Reproducing figure 1
To reproduce figure 1 for example, use the following procedures:
- Launch `./Production-Figure-Code/makeFigurePanels.m` using MATLAB. Note that the code has been tested with MATLAB `R2017a` and `R2018a`.
- Specify input parameters
  - Set `targetFig = 'fig1'` to select figure 1 as the target
  - Set `USE_PUBLISHED_DATASET = 1` to show the data we published. Alternatively, use `USE_PUBLISHED_DATASET = 0` if the simulation data for figure 1 has been created through step 1
- Run the MATLAB code

You should see a new MATLAB figure containing Figure 1 panels. PDF(s) will be saved under `./FigureOutput/fig1/`.
