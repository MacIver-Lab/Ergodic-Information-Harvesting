# Code and data to reproduce results from "Sense organ control in moths to moles is a gamble on information through motion" by Chen Chen, Todd D. Murphey, and Malcolm A. MacIver, Northwestern University, Evanston IL, USA

# Ergodic Information Harvesting (EIH) Video & Tutorial
For a prior publication (["Ergodic Exploration of Distributed Information"](https://nxr.northwestern.edu/sites/default/files/publications/Mill16a_ergodic_control_distributed_info.pdf), Miller et al., IEEE Trans. Robotics, 2016) we made a video that goes through key steps of the EIH algorithm. For the ergodic harvesting study, we've also developed an interactive tutorial to cover some of the concepts.
- Video of how EIH works (Two differences from ["Ergodic Exploration of Distributed Information"](https://nxr.northwestern.edu/sites/default/files/publications/Mill16a_ergodic_control_distributed_info.pdf): 1. Object to be found is not moving; 2. Uses Fisher Information instead of entropy.) [Click here to watch the video on YouTube](https://youtu.be/eF6J-YmPdIA)
- Interactive Jupyter notebook tutorial [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/MacIver-Lab/Ergodic-Information-Harvesting/master?filepath=Tutorial%2FErgodic_Information_Harvesting_Tutorial.ipynb)

# Steps to reproduce the results shown in the EIH paper
All of the simulation code is written with [Python 3](https://www.python.org/). All of the figure plotting files are written in MATLAB (R2017a+). The code can be run on:
- A local computer, which is very easy to set up but the performance is ultimately limited by the number of locally accessible CPU cores
- Cloud computing virtual servers throug any Infrastructure as a Service (IaaS) provider, *e.g.* [Amazon Elastic Compute Cloud](https://aws.amazon.com/ec2/), [Google Cloud Compute Engine](https://cloud.google.com/compute/), or academic [HPCC (High Performance Computing Cluster)](https://en.wikipedia.org/wiki/HPCC) systems. Cloud computing is easy to setup and provides a way to scale up the total number of running threads (*e.g.* Google Cloud Compute Engine allows up to 96 CPU threads per instance). Our code's runtime envionrment and dependencies are fully containerized through [Singularity](https://www.sylabs.io/singularity/) to minimize the effort for envionrment set and scales easily for better performance.

Note that this repository has included all the published data, including all the simulations, to reproduce all figures in the paper and supplemental materials. For lazy result reproduction and verification, simply jump to step 3 to directly reproduce the figures using published dataset.

## Detailed Steps to Reproduce
To accomodata the possible [dependencies hell](https://en.wikipedia.org/wiki/Dependency_hell) and minimize the effort of setting up runtime environment, we prebuilt [container image](https://en.wikipedia.org/wiki/Container_(virtualization)) to be used for all the simulation code execution in Python. Note that this is only for reproducing simulations, a local installation of MATLAB (not provided in the container) is still required to reproduce figures.

### Download or Clone This Repository
To start, you will need to download the most recent version of the EIH repository. This can be done by either cloning it as a `git` repository, which can be done by executing:
```bash
git clone --depth=1 https://github.com/MacIver-Lab/Ergodic-Information-Harvesting
```
or, simply [download this repository as an archive from GitHub](https://github.com/MacIver-Lab/Ergodic-Information-Harvesting/archive/master.zip).

### Install Singularity and Pull the EIH Container Image
[Singularity](https://www.sylabs.io/singularity/) is a software for building, managing, and executing container images. It is required to run the simulations through our prebuilt container image. To install Singularity, follow the [official installation guide for Windows/MacOS/Linux](https://www.sylabs.io/guides/2.6/user-guide/installation.html).

Once Singularity is installed, open a command line tool at the EIH directory `./Ergodic-Information-Harvesting/` and pull the prebuilt EIH container image from cloud by running the following command in the command line:
```bash
singularity --name EIH.img pull shub://MacIver-Lab/Ergodic-Information-Harvesting
```

### Invoke Shell in the EIH Container Image
The container image is basically a fully self-contained Linux OS image with Python 3 dependencies setup for EIH simulation. We will invoke the command line tool inside of the EIH container image to interact with the resources inside to start our simulations.

First, invoke the shell inside of the image by running:
```bash
singularity shell -B ./:/EIH ./EIH.img
```

We used [Cython](https://cython.org/) to accelerate the simulation which requires compiling some of the code before running the simulation. Compile the accelerated code by calling the following command:
```bash
. ./BuildCython.sh
```

### Start Reproducing Simulation
You are all set for the environment setup. You can start reproducing all the simulation results by running the main simulation code:
```bash
cd SimulationCode/
Python3 RunAllSims.py
```


## TODO - Organize
**NOTE**: Depending on your operating system, you **may need to prevent your system from going to sleep**. This is necessary with MacOS. With MacOS: Open a terminal, and type `caffeinate` and hit return. Your system will be prevented from sleeping until you hit Control-C.

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
As noted above, sleeping causes the Jupyter notebook to lose connection and stops the simulation. Also, 'sm-fig4' is a ten hour run in MATLAB to regenerate simulation results. To prevent MacOS from sleeping in these instances, use `caffeinate` at a Terminal window before starting the Jupyter notebook or before regenerating 'sm-fig4' within MATLAB.


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
