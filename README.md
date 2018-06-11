# Code and data to reproduce results from "Ergodic information harvesting as a behavioral strategy for complex environments" by Chen Chen, Todd D. Murphey, and Malcolm A. MacIver, Northwestern University, Evanston IL, USA

# Ergodic Information Harvesting (EIH) Video & Tutorial
For a prior publication (["Ergodic Exploration of Distributed Information"](https://nxr.northwestern.edu/sites/default/files/publications/Mill16a_ergodic_control_distributed_info.pdf), Miller et al., IEEE Trans. Robotics, 2016) we made a video that goes through key steps of the EIH algorithm. For the ergodic harvesting study, we've also developed an interactive tutorial to cover some of the concepts.
- Video of how EIH works (Two differences from study: 1. Object to be found is not moving; 2. Uses Fisher Information instead of entropy.) [Click here to watch the video on YouTube](https://youtu.be/QZ9fGYmJ0G0)
- Interactive Jupyter notebook tutorial [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/MacIver-Lab/Ergodic-Information-Harvesting/master?filepath=Tutorial%2FErgodic_Information_Harvesting_Tutorial.ipynb)


# Steps to reproduce the results shown in the EIH paper
All of the simulation code is written with Python 3.6 using [Jupyter Notebook](http://jupyter.org/). All of the figure plotting files are written in MATLAB (R2017a). The code can be run on:
- A local computer, which is very easy to set up but limited by the number of accessible CPU cores.
- Cloud computing virtual servers throug any popular Infrastructure as a Service (IaaS) provider, *e.g.* [Amazon Elastic Compute Cloud](https://aws.amazon.com/ec2/) or [Google Cloud Compute Engine](https://cloud.google.com/compute/). Cloud computing is easy to setup and provides a way to scale up the total number of running threads (*e.g.* Google Cloud Compute Engine allows up to 96 CPU threads per instance). A [cloud computing setup guide](https://github.com/MacIver-Lab/Ergodic-Information-Harvesting/blob/master/CloudComputingSetup.md) covers the steps needed to setup EIH simulation on a Linux instance running on the cloud.

## Setting up the runtime environment
The results generation system uses Jupyter as the environment for interaction. Jupyter is very easy to set up if you have not used it before. Below is a quick start for those who does not yet have Jupyter environment setup locally.

### Install Anaconda 
[Anaconda](https://www.anaconda.com/download/) is the required Python environment. It runs on MacOS, Unix, and Windows environments and provide easy package management for setting up the runtime environments for the simulation. To download, go to [https://www.anaconda.com/download/](https://www.anaconda.com/download/) and install Anaconda.

Once installed, open the Anaconda Prompt or bash and run the following code to install all of the dependencies.
The required packages are `scipy` and `numpy`. Note that a specific version has been specified to avoid issues due to any future update to these packages that may cause a backward incompatibility with our code. Most of the other dependencies, such as `jupyter` and `IPython`, has been satisfied by the Anaconda bundle already and are thus not included in the code below.
```bash
conda install scipy=1.0.1 numpy=1.14.3
```
Alternatively, if you prefer `pip`, you could use the code below instead.
```bash
pip install scipy==1.0.1 numpy==1.14.3
```

### Launch Jupyter
Once the environment is setup, you can launch the simulation through a `jupyter notebook` client. Here's how to do it.
- Navigate to the simulation code folder `/SimulationCode/`
- Launch a new Jupyter notebook instance
  - If using Windows, open a command prompt under this code folder and run `jupyter notebook`
  - If using Linux or MacOS, simply run this command in a terminal: `jupyter notebook` 
- A new browser tab should pop up. Find and open `Ergodic-Information-Harvesting-Simulation.ipynb` under the file list. Note that Jupyter notebook will use the current directory of your command prompt as the working directory: if you don't see the code folder, that means your command window is in a different folder and you need to navigate to `/SimulationCode/` before running `jupyter notebook`
- Once opened, you should see the code for simulation and you are now good to go to reproduce the results of our study.
- **NOTE**: Depending on your operating system, you **may need to prevent your system from going to sleep**. This is necessary with MacOS. With MacOS: Open a terminal, and type `caffeinate` and hit return. Your system will be prevented from sleeping until you hit Control-C.

## How to Reproduce Figure Results
There are two stages required to reproduce the published figure results. First, follow the simulation section below to run EIH simulation trials to reproduce the data required for figures. Then proceed to the figure plotting code to reproduce the figure results. 

### Benchmark Running Time
**The total run time to run both the simulation and figure plotting code is 45.6 hours for a 2015 MacOS desktop system (iMac 2015, Intel i7 Quad Core with 4.4GHz turboboost, running with `nThreads = 8`). Most of this time is for Supplementary Figures 1 and 4, which have close to 200 simulations. Therefore, in addition to the full simulation code provided and detailed in Step 1, we have also included our published dataset which can be viewed by setting a flag in the figure generation code in Step 2 below, which allows immediate reproduction of published figures. This also provides a baseline of comparison for results that you generate.**

|   Figure  | Total Running Time (Hours) |
|:---------:|:--------------------------:|
|   `fig1`  |             1.1            |
|   `fig2`  |             3.3            |
|   `fig3`  |             1.1            |
| `sm-fig1` |             29             |
| `sm-fig4` |             10             |
| `sm-fig5` |             1.1            |

### Additional Note for Linux and MacOS Users
#### Prevent System from Sleeping During Simulation
As noted above, sleeping causes the Jupyter notebook to lose connection and stops the simulation. To prevent MacOS from sleeping, use `caffeinate` at a Terminal window before submitting any local simulation.

#### Be Sure to Use the Right Version of Python
For `sm-fig4`, we use MATLAB to initiate batch python simulations which might requires additional care. For some systems you might see errors while reproducing `sm-fig4` in MATLAB. This is due to python version mismatch which is detailed in the paragraph below (for those who want to know why this happens):

> Since Linux and Unix has built-in python in addition to Anaconda, please be sure to add Anaconda to the system's path and use the python that comes with Anaconda instead of the system's default python. For `sm-fig4`, the python simulation is handled through MATLAB using the `system()` command to invoke bash commands. We have found for some systems, calling `system('python -V');` in MATLAB will invoke the default system python under `/usr/lib/python*` instead of Anaconda's python, which is typically located under `$HOME/anaconda/bin`.

The solution is simple. Find the path where you installed Anaconda. If you have added Anaconda to the system's PATH you can find it easily with `which conda`, which will be Anaconda's install path. If the previous command did not run, that means you either have not installed Anaconda or it has not yet been added to your PATH. In the former case, please install Anaconda and be sure to say `yes` when asked whether or not to add Anaconda to system's PATH. In the latter case, you need to first find where Anaconda's python was installed (the default path is `$HOME/anaconda/bin`), then use `export PATH="$HOME/anaconda/bin:$PATH"` to add to system's path.

### Step 1 - Local EIH Simulation (optional)
#### Code Structure
The simulation code files are all stored under `./SimulationCode/` and are organized in a centralized fashion. `Ergodic-Information-Harvesting-Simulation.ipynb` is the only notebook file you need to run and you can use it to reproduce the raw simulation data used for all but one of the figures in the paper and Supplementary Information (all but `sm-fig4`, which uses a combination of MATLAB and python). Each dataset is organized in a per-figure fashion and the parameter supporting the simulation is stored under `/SimulationCode/FigParameters/` folder in `json` format.

All of the data used in published figures will be simulated through this step (except for `sm-fig4`: see Step 2).

#### Example Procedure of Simulating Data for figure 1
To simulate the raw data for a given figure, figure 1 for instance, you just need to:
- specify `targetFigure = ['fig1']` in the `Ergodic-Information-Harvesting-Simulation.ipynb` which is further documented inside the notebook. To run multiple figures, simply append the list with the desired figures, *e.g.* `
targetFiguretargetF  = ['fig1', 'fig2', 'fig3', 'sm-fig1']`.
- specify `nThreads` parameter based on the number of CPU threads available on your machine. Use `nThreads = cpu_count()` will ensure the best performance but will slow down your machine significantly on some figures which utilize all of the threads. Use `nThreads = cpu_count() - 1` to leave 1 thread for your normal workflow.
- run the entire notebook by selecting on the top menu `Cell -> Run All`. 

Depending on the size of the simulation and the number of threads you allow the program to use, it could take some time to finish. 

Once completed (the program will display `All done! EOF at ...` at the bottom), you can find the simulated data under `/SimulationCode/SimData/fig1/` (`/fig1` because we are simulating for figure 1 here). The raw data will be in `*.mat` MATLAB data format that can be loaded by MATLAB.

### Step 2 - Reproduce Figure Results
#### Code Structure
The figure code files are under `./Production-Figure-Code/`. Similarly, it is centralized in a single MATLAB file `makeFigurePanels.m`. Code for each figure panels are included under `./Production-Figure-Code/FigureCode/` along with the previously simulated data if step 1 is skipped.

Note that running `sm-fig4` with `USE_PUBLISHED_DATASET = 0` (which uses locally reproduced dataset instead of the previously simulated) will start a batch simulation before making figure panels. The simulations for this figure took around 10 hours on the benchmark machine.

#### Example Procedure of Reproducing figure 1
To reproduce figure 1 for example, use the following procedures:
- Launch `./Production-Figure-Code/makeFigurePanels.m` using MATLAB. Note that the code has been tested with MATLAB `R2017a` and `R2018a`.
- Specify input parameters
  - Set `targetFig = 'fig1'` to select figure 1 as the target
  - Set `USE_PUBLISHED_DATASET = 1` to show the data we published. Alternatively, use `USE_PUBLISHED_DATASET = 0` if the simulation data for figure 1 has been created through step 1
- Run the MATLAB code

You should see a new MATLAB figure containing Figure 1 panels. PDF(s) will be saved under `./FigureOutput/fig1/`.
