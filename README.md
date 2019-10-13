# Code and data to reproduce results from "Sense organ control in moths to moles is a gamble on information through motion" by Chen Chen, Todd D. Murphey, and Malcolm A. MacIver, Northwestern University, Evanston IL, USA


[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2511)
[![Docker status](https://img.shields.io/docker/cloud/build/maciverlabnu/ergodic-information-harvesting.svg)](https://hub.docker.com/r/maciverlabnu/ergodic-information-harvesting)
[![Docker pulls](https://img.shields.io/docker/pulls/maciverlabnu/ergodic-information-harvesting.svg)](https://hub.docker.com/r/maciverlabnu/ergodic-information-harvesting)
[![nbviewer](https://camo.githubusercontent.com/bfeb5472ee3df9b7c63ea3b260dc0c679be90b97/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f72656e6465722d6e627669657765722d6f72616e67652e7376673f636f6c6f72423d66333736323626636f6c6f72413d346434643464)](https://nbviewer.jupyter.org/github/MacIver-Lab/Ergodic-Information-Harvesting/blob/master/Tutorial/Ergodic_Information_Harvesting_Tutorial.ipynb)

# Ergodic Information Harvesting (EIH) Video & Tutorial
- Video explainer of key concepts underlying ergodic information harvesting ([link to YouTube](https://youtu.be/vOAzSa6Pna4)), and an example application of the EIH algorithm to control an underwater electrolocation robot ([link to YouTube](https://youtu.be/kHKvaxuLnq4)).
- Jupyter notebook tutorial on how the ergodic information harvesting algorithm computes the expected information density, click to view online: [![nbviewer](https://camo.githubusercontent.com/bfeb5472ee3df9b7c63ea3b260dc0c679be90b97/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f72656e6465722d6e627669657765722d6f72616e67652e7376673f636f6c6f72423d66333736323626636f6c6f72413d346434643464)](https://nbviewer.jupyter.org/github/MacIver-Lab/Ergodic-Information-Harvesting/blob/master/Tutorial/Ergodic_Information_Harvesting_Tutorial.ipynb)

# Steps to reproduce the results shown in the EIH paper
All of the simulation code is written with [Python 3](https://www.python.org/). All of the figure plotting files are written in MATLAB (R2017a+). The code can be run on:
- A local computer, which is very easy to set up but the performance is ultimately limited by the number of locally accessible CPU cores
- Cloud computing virtual servers through any Infrastructure as a Service (IaaS) provider, *e.g.* [Amazon Elastic Compute Cloud](https://aws.amazon.com/ec2/), [Google Cloud Compute Engine](https://cloud.google.com/compute/), or academic [HPCC (High Performance Computing Cluster)](https://en.wikipedia.org/wiki/HPCC) systems. Cloud computing is easy to setup and provides a way to scale up the total number of running threads (*e.g.* Google Cloud Compute Engine allows up to 96 CPU threads per instance). Our code's runtime environment and dependencies are fully containerized through [Singularity](https://www.sylabs.io/singularity/) to minimize the effort needed to re-run our simulations and data analyses and to permit easy scaling to reduce runtime.

Note that this repository has included all of the published data, including all the simulations, to reproduce all figures in the paper and supplemental materials. To reproduce our figures from the published data, rather than re-run all the simulations from scratch, simply jump to step 5 below.

## Detailed Steps
To avoid possible [dependency hell](https://en.wikipedia.org/wiki/Dependency_hell) and minimize the effort of setting up the runtime environment we used for our results, we prebuilt a [container image](https://en.wikipedia.org/wiki/Container_(virtualization)) to be used for executing all the simulation code in Python. Here is an article explaining the utility of containers for reproducibility of research: [Singularity: Scientific containers for mobility of compute](https://doi.org/10.1371/journal.pone.0177459). Note that this is only for reproducing simulations: for generation of the figures from the simulations, a local installation of MATLAB (not provided in the container) is still required.

### 1. Clone This Repository
To start, you will need to clone the most recent version of the EIH repository. Two tools are required:
- `git` - is used to pull all the non-data files from this repository. Go to [git's official release page](https://git-scm.com/downloads) to download and install git. Then use the following command to clone this repository:
 ```bash
 git clone --depth=1 https://github.com/MacIver-Lab/Ergodic-Information-Harvesting
 ```
- `git-lfs` - is used to pull all the published data. Go to [git-lfs's official release page](https://git-lfs.github.com/) to download and install. Then run the following command **inside the root directory of the cloned EIH repo `./Ergodic-Information-Harvesting/`**
```bash
cd Ergodic-Information-Harvesting
git lfs install
git lfs pull
```
If you succeeded, you should see files being downloaded by `git-lfs`.

### 2. Install Singularity or Docker and Pull the EIH Container Image
The following steps are different depending on the operating system environment. Please follow the instructions that best suits your system.

#### Linux and Linux based HPCC/Cloud Environments
[Singularity](https://www.sylabs.io/singularity/) is software for building, managing, and executing container images on Linux. It is required to run the simulations through our prebuilt container image on Linux and is easier to setup. For HPCC (High Performance Computing Cluster) environments, Singularity is usually used ahead of Docker for security reasons. To install Singularity, follow the [official installation guide for Linux](https://www.sylabs.io/guides/2.6/user-guide/installation.html).

Once Singularity is installed, open a command line tool at the EIH directory `./Ergodic-Information-Harvesting/` and pull the prebuilt EIH container image from the cloud by running the following command on the command line:
```bash
singularity pull --name EIH.img shub://MacIver-Lab/Ergodic-Information-Harvesting
```

#### Windows and MacOS
[Docker](https://en.wikipedia.org/wiki/Docker_%28software%29) is used on Windows and MacOS for running the EIH container image since Singularity only has native support for Linux OS. To install Docker, follow the official installation guide for [Windows](https://hub.docker.com/editions/community/docker-ce-desktop-windows) or [MacOS](https://hub.docker.com/editions/community/docker-ce-desktop-mac).

Once installed, first go to the Docker settings and update the desired number of CPUs and memory available for Docker (and EIH simulation) to use. The number of CPU thread available during the simulation is limited by this configuration.

![](https://docs.docker.com/docker-for-windows/images/settings-advanced.png)

Next, open a command line tool at any directory and pull the prebuilt EIH container image from the cloud by running the following command on the command line:
```bash
docker pull maciverlabnu/ergodic-information-harvesting
```

### 3. Invoke Shell in the EIH Container Image
The container image is a fully self-contained Linux OS image with Python 3 dependencies setup for generating the EIH simulations developed for the study. We will invoke the command line tool inside of the EIH container image to interact with the resources inside the container and start the simulations.

- When using Singularity on Linux OS, invoke the shell inside of the image by running:
 ```bash
 singularity shell -B ./:/EIH ./EIH.img
 ```
- When using Docker on Windows or MacOS, first check the drive where the EIH repository is downloaded to is accessible to Docker:

 ![](https://msdnshared.blob.core.windows.net/media/2016/06/d4w-shared-drives.png)
 
 make sure the target drive is checked before invoking the shell. 
 
 ```bash
 docker run -it -v [absolute path to EIH folder]:/EIH maciverlabnu/ergodic-information-harvesting
 ```
 replace the `[absolute path to EIH folder]` part with the absolute path to your local EIH repository folder, *e.g.* `C:/Ergodic-Information-Harvesting` (remember to replace `\` with `/` when in Windows) or `~/Ergodic-Information-Harvesting`.
 
### 4. Start Reproducing Simulation

We used [Cython](https://cython.org/) to accelerate the simulation which requires compiling some of the code before running the simulation. Compile the accelerated code by calling the following command (this only needs to be done once):
```bash
cd /EIH
chmod +x ./BuildCython.sh
. ./BuildCython.sh
```

You are all set for the environment setup. As the code will take several days to run (see benchmarks below), you may need to **prevent your system from going to sleep**. This is necessary with MacOS. With MacOS: Open a terminal, and type `caffeinate` and hit return. Your system will be prevented from sleeping until you hit Control-C. You can start reproducing all the simulation results by running the main simulation code:
```bash
cd /EIH/SimulationCode/
python3 RunAllSims.py
```
By default, `RunAllSims.py` will check the number of available CPU threads and automally run parallel simulation jobs with the maximum number of threads possible. Nonetheless, the number of threads can be manually specified by passing the desired parallel thread count argument, for example
```bash
python3 RunAllSims.py 20
```
will run 20 threads in parallel.

Once all the simulation jobs are done, exit the Singularity shell environment by calling the `exit` command. 

### 5. Reproduce Figure Results
The figure generation code is written in MATLAB so MATLAB R2017a or a more recent version is required. To start, open the `makeFigurePanels.m` code in MATLAB under the `Production-Figure-Code` folder. 
- Launch `Ergodic-Information-Harvesting/Production-Figure-Code/makeFigurePanels.m` using MATLAB. Note that the code has been tested with MATLAB `R2017a` and `R2018a`.
- Specify input parameters
  - Specify which figures to reproduce. Default: `targetFig = 'all'`. Example for one figure: set `targetFig = 'fig2'` to select figure 2 as the target.
  - Specify whether to generate figures on the published data, or your locally simulated data. Default: `USE_PUBLISHED_DATASET = 1`, which uses the published dataset included in the repository. Alternatively, set `USE_PUBLISHED_DATASET = 0`; this will force the code to use your simulation results located at `Ergodic-Information-Harvesting/SimulationCode/SimData/`. 
- Run the MATLAB code

PDF(s) will be saved under `Ergodic-Information-Harvesting/Production-Figure-Code/FigureOutput`.

### Benchmark Running Time
- **Benchmark on PC**: `~98.2 hours` on a 2015 MacOS desktop system (iMac 2015, Intel i7 Quad Core with 4.4GHz turboboost, running with Docker preferences set to use 4 cores (see above))
- **Benchmark on HPCC**: `~10.79 hours` on Northwestern University QUEST HPCC system's 8th generation computing node with single Intel Xeon Gold 6132 2.6 GHz Intel® QPI 2666 MHz with 3.70 GHz turboboost (Linux OS, 28 physical cores running with `nThreads = 56`)

### Additional Note for Linux and MacOS Users
#### Prevent System from Sleeping During Simulation
To prevent MacOS from sleeping in these instances, use `caffeinate` at a Terminal window running simulation jobs.
