# Code and data to reproduce results from "Sense organ control in moths to moles is a gamble on information through motion" by Chen Chen, Todd D. Murphey, and Malcolm A. MacIver, Northwestern University, Evanston IL, USA


[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2511)
[![nbviewer](https://camo.githubusercontent.com/bfeb5472ee3df9b7c63ea3b260dc0c679be90b97/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f72656e6465722d6e627669657765722d6f72616e67652e7376673f636f6c6f72423d66333736323626636f6c6f72413d346434643464)](https://nbviewer.jupyter.org/github/MacIver-Lab/Ergodic-Information-Harvesting/blob/master/Tutorial/Ergodic_Information_Harvesting_Tutorial.ipynb)
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/MacIver-Lab/Ergodic-Information-Harvesting/master?filepath=Tutorial%2FErgodic_Information_Harvesting_Tutorial.ipynb)
 
# Ergodic Information Harvesting (EIH) Video & Tutorial
- Video explainer of key concepts underlying ergodic information harvesting, and an example application of the EIH algorithm to controlling an underwater electrolocation robot [Click here to watch the video on YouTube](https://youtu.be/eF6J-YmPdIA)
- Interactive Jupyter notebook tutorial, click to view online: [![nbviewer](https://camo.githubusercontent.com/bfeb5472ee3df9b7c63ea3b260dc0c679be90b97/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f72656e6465722d6e627669657765722d6f72616e67652e7376673f636f6c6f72423d66333736323626636f6c6f72413d346434643464)](https://nbviewer.jupyter.org/github/MacIver-Lab/Ergodic-Information-Harvesting/blob/master/Tutorial/Ergodic_Information_Harvesting_Tutorial.ipynb)
  
  or use  [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/MacIver-Lab/Ergodic-Information-Harvesting/master?filepath=Tutorial%2FErgodic_Information_Harvesting_Tutorial.ipynb) to run interactively through online Jupyter Notebook

# Steps to reproduce the results shown in the EIH paper
All of the simulation code is written with [Python 3](https://www.python.org/). All of the figure plotting files are written in MATLAB (R2017a+). The code can be run on:
- A local computer, which is very easy to set up but the performance is ultimately limited by the number of locally accessible CPU cores
- Cloud computing virtual servers throug any Infrastructure as a Service (IaaS) provider, *e.g.* [Amazon Elastic Compute Cloud](https://aws.amazon.com/ec2/), [Google Cloud Compute Engine](https://cloud.google.com/compute/), or academic [HPCC (High Performance Computing Cluster)](https://en.wikipedia.org/wiki/HPCC) systems. Cloud computing is easy to setup and provides a way to scale up the total number of running threads (*e.g.* Google Cloud Compute Engine allows up to 96 CPU threads per instance). Our code's runtime envionrment and dependencies are fully containerized through [Singularity](https://www.sylabs.io/singularity/) to minimize the effort for envionrment set and scales easily for better performance.

Note that this repository has included all the published data, including all the simulations, to reproduce all figures in the paper and supplemental materials. For lazy result reproduction and verification, simply jump to step 5 to directly reproduce the figures using published dataset.

## Detailed Steps to Reproduce
To avoid possible [dependencies hell](https://en.wikipedia.org/wiki/Dependency_hell) and minimize the effort of setting up the runtime environment we used for our results, we prebuilt a [container image](https://en.wikipedia.org/wiki/Container_(virtualization)) to be used for executing all the simulation code in Python. Note that this is only for reproducing simulations: for generation of the figures from the simulations, a local installation of MATLAB (not provided in the container) is still required.

### 1. Download or Clone This Repository
To start, you will need to download the most recent version of the EIH repository. This can be done by either cloning it as a `git` repository, as follows:
```bash
git clone --depth=1 https://github.com/MacIver-Lab/Ergodic-Information-Harvesting
```
or, you can simply [download this repository as an archive from GitHub](https://github.com/MacIver-Lab/Ergodic-Information-Harvesting/archive/master.zip).

### 2. Install Singularity and Pull the EIH Container Image
[Singularity](https://www.sylabs.io/singularity/) is software for building, managing, and executing container images. It is required to run the simulations through our prebuilt container image. To install Singularity, follow the [official installation guide for Windows/MacOS/Linux](https://www.sylabs.io/guides/2.6/user-guide/installation.html).

Once Singularity is installed, open a command line tool at the EIH directory `./Ergodic-Information-Harvesting/` and pull the prebuilt EIH container image from the cloud by running the following command on the command line:
```bash
singularity --name EIH.img pull shub://MacIver-Lab/Ergodic-Information-Harvesting
```

### 3. Invoke Shell in the EIH Container Image
The container image is a fully self-contained Linux OS image with Python 3 dependencies setup for generating the EIH simulations developed for the study. We will invoke the command line tool inside of the EIH container image to interact with the resources inside the container and start the simulations.

First, invoke the shell inside of the image by running:
```bash
singularity shell -B ./:/EIH ./EIH.img
```

We used [Cython](https://cython.org/) to accelerate the simulation which requires compiling some of the code before running the simulation. Compile the accelerated code by calling the following command (this only needs to be done once):
```bash
cd /EIH
. ./BuildCython.sh
```

### 4. Start Reproducing Simulation
You are all set for the environment setup. You can start reproducing all the simulation results by running the main simulation code:
```bash
cd /EIH/SimulationCode/
python3 RunAllSims.py
```
By default, `RunAllSims.py` will check the number of available CPU threads and automally run parallel simulation jobs with the maximum number of threads possible. Nonetheless, the number of threads can be manually specified by passing the desired parallel thread count argument to it, for example
```bash
python3 RunAllSims.py 20
```
will run 20 threads in parallel.

**NOTE**: The simulation will take a long time to finish. Depending on your operating system, you may need to **prevent your system from going to sleep**. This is necessary with MacOS. With MacOS: Open a terminal, and type `caffeinate` and hit return. Your system will be prevented from sleeping until you hit Control-C.

Once all the simulation jobs are done, exit the Singularity shell environment by calling the `exit` command. 

### 5. Reproduce Figure Results
The figure generation code is written in MATLAB so MATLAB R2017a or a more recent version is required. To start, open the `makeFigurePanels.m` code in MATLAB under the `Production-Figure-Code` folder. To reproduce figure 2, for example, use the following procedure:
- Launch `Ergodic-Information-Harvesting/Production-Figure-Code/makeFigurePanels.m` using MATLAB. Note that the code has been tested with MATLAB `R2017a` and `R2018a`.
- Specify input parameters
  - Set `targetFig = 'fig2'` to select figure 2 as the target
  - Set `USE_PUBLISHED_DATASET = 1` to use the published dataset included in the repository. Alternatively, if the local simulation jobs are completed, use of `USE_PUBLISHED_DATASET = 0` will force the code to use reproduced data located at `Ergodic-Information-Harvesting/SimulationCode/SimData/`
- Run the MATLAB code

You should see a new MATLAB figure containing Figure 2 panels. PDF(s) will be saved under `Ergodic-Information-Harvesting/Production-Figure-Code/FigureOutput/fig2/`.

### Benchmark Running Time
**The total run time to run both the simulation and figure plotting code is TBA hours for a 2015 MacOS desktop system (iMac 2015, Intel i7 Quad Core with 4.4GHz turboboost, running with `nThreads = 8`).**

### Additional Note for Linux and MacOS Users
#### Prevent System from Sleeping During Simulation
To prevent MacOS from sleeping in these instances, use `caffeinate` at a Terminal window running simulation jobs.
