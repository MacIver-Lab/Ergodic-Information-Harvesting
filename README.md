# Code and data to reproduce results from "Tuning movement for sensing in an uncertain world" by Chen Chen, Todd D. Murphey, and Malcolm A. MacIver, Northwestern University, Evanston IL, USA


[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2511)
[![Docker status](https://img.shields.io/docker/cloud/build/maciverlabnu/ergodic-information-harvesting.svg)](https://hub.docker.com/r/maciverlabnu/ergodic-information-harvesting)
[![Docker pulls](https://img.shields.io/docker/pulls/maciverlabnu/ergodic-information-harvesting.svg)](https://hub.docker.com/r/maciverlabnu/ergodic-information-harvesting)
[![nbviewer](https://camo.githubusercontent.com/bfeb5472ee3df9b7c63ea3b260dc0c679be90b97/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f72656e6465722d6e627669657765722d6f72616e67652e7376673f636f6c6f72423d66333736323626636f6c6f72413d346434643464)](https://nbviewer.jupyter.org/github/MacIver-Lab/Ergodic-Information-Harvesting/blob/master/Tutorial/Ergodic_Information_Harvesting_Tutorial.ipynb)
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/MacIver-Lab/Ergodic-Information-Harvesting/master?filepath=Tutorial%2FErgodic_Information_Harvesting_Tutorial.ipynb)
 
# Ergodic Information Harvesting (EIH) Video & Tutorial
- Video explainer of key concepts underlying ergodic information harvesting, and an example application of the EIH algorithm to control an underwater electrolocation robot: [Link to YouTube](https://youtu.be/eF6J-YmPdIA)
- Interactive Jupyter notebook tutorial, click to view online: [![nbviewer](https://camo.githubusercontent.com/bfeb5472ee3df9b7c63ea3b260dc0c679be90b97/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f72656e6465722d6e627669657765722d6f72616e67652e7376673f636f6c6f72423d66333736323626636f6c6f72413d346434643464)](https://nbviewer.jupyter.org/github/MacIver-Lab/Ergodic-Information-Harvesting/blob/master/Tutorial/Ergodic_Information_Harvesting_Tutorial.ipynb)
  
  or use  [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/MacIver-Lab/Ergodic-Information-Harvesting/master?filepath=Tutorial%2FErgodic_Information_Harvesting_Tutorial.ipynb) to run interactively through online Jupyter Notebook

# Steps to reproduce the results shown in the EIH paper
All of the simulation code is written with [Python 3](https://www.python.org/). All of the figure plotting files are written in MATLAB (R2017a+). The code can be run on:
- A local computer, which is very easy to set up but the performance is ultimately limited by the number of locally accessible CPU cores
- Cloud computing virtual servers through any Infrastructure as a Service (IaaS) provider, *e.g.* [Amazon Elastic Compute Cloud](https://aws.amazon.com/ec2/), [Google Cloud Compute Engine](https://cloud.google.com/compute/), or academic [HPCC (High Performance Computing Cluster)](https://en.wikipedia.org/wiki/HPCC) systems. Cloud computing is easy to setup and provides a way to scale up the total number of running threads (*e.g.* Google Cloud Compute Engine allows up to 96 CPU threads per instance). Our code's runtime environment and dependencies are fully containerized through [Docker](https://www.docker.com) to minimize the effort needed for environment setup, and for easy scaling to run the code faster (you can run the code on a single CPU over a few days, or on many CPUS on the cloud in a few hours). Setting the code up to run on a cloud service for the first time is somewhat involved if you are not used to doing this. We have a set of screencasts to walk through doing this from scratch on Amazon Web Services that we are happy to share.

Note that this repository has included all of the published data, including all the simulations, to reproduce all figures in the paper. To reproduce our figures from the published data, rather than re-run all the simulations from scratch, simply jump to step 5 below.

## Detailed Steps
To avoid possible [dependency hell](https://en.wikipedia.org/wiki/Dependency_hell) and minimize the effort of setting up the runtime environment we used for our results, we prebuilt a [container image](https://en.wikipedia.org/wiki/Container_(virtualization)) to be used for executing all the simulation code in Python using [Docker](https://docs.docker.com/get-started/overview/). Here is an article explaining the utility of containers for reproducibility of research: [Singularity: Scientific containers for mobility of compute](https://doi.org/10.1371/journal.pone.0177459). Note that this is only for reproducing simulations: for generation of the figures from the simulations, a local installation of MATLAB (not provided in the container) is still required.

### 1. Clone This Repository
To start, you will need to clone the most recent version of the EIH repository. Two tools are required:
- `git` - is used to pull all the non-data files from this repository. Go to [git's official release page](https://git-scm.com/downloads) to download and install git. Then use the following command to clone this repository:
 ```bash
 git clone --depth=1 https://github.com/MacIver-Lab/Ergodic-Information-Harvesting
 ```
- `git-lfs` - is used to pull all the published data (required to reproduce the results). Go to [git-lfs's official release page](https://git-lfs.github.com/) to download and install. Then run the following command **inside the root directory of the cloned EIH repo `./Ergodic-Information-Harvesting/`**
```bash
cd Ergodic-Information-Harvesting
git lfs install
git lfs pull
```
If you succeeded, you should see files being downloaded by `git-lfs`. Once it is setup, should you decide to delete the files and start again, you should only need to do the `git clone` step.

### 2. Install Docker and Pull the EIH Container Image

- Install Docker by following the official documentations: https://docs.docker.com/engine/install/

- (Optionally) For Linux or Linux-based HPCC/cloud computing environments, please additionally follow the [post-installation setup steps for Linux](https://docs.docker.com/engine/install/linux-postinstall/) to allow running docker without `sudo`. If you don't want or unable to follow this step, you will need to make sure to run docker commands with `sudo docker` rather than `docker`.

### 3. Invoke Shell in the EIH Container Image
The container image is a fully self-contained Linux OS image with Python 3 dependencies setup for generating the EIH simulations developed for the study. We will invoke the command line tool inside of the EIH container image to interact with the resources inside the container and start the simulations.

```bash
docker run -it -v [absolute path to EIH folder]:/EIH maciverlabnu/ergodic-information-harvesting
```

replace the `[absolute path to EIH folder]` part with the absolute path to your local EIH repository folder, *e.g.* `C:/Ergodic-Information-Harvesting` (remember to replace `\` with `/` when in Windows) or `~/Ergodic-Information-Harvesting`.

If you are already inside the Ergodic-Information-Harvesting folder, you can simply do
```bash
docker run -it -v $PWD:/EIH maciverlabnu/ergodic-information-harvesting
```
 
### 4. Start Reproducing Simulation

We used [Cython](https://cython.org/) to accelerate the simulation which requires compiling some of the code before running the simulation. Compile the accelerated code by calling the following command (this only needs to be done once):
```bash
# Command to run inside the container
cd /EIH
chmod +x ./BuildCython.sh
./BuildCython.sh
```

You are now all set for the environment setup. You can start reproducing all the simulation results by running the main simulation code:
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

The figure generation code is written in MATLAB and MATLAB R2017a or a more recent version is required. To start, open the `makeFigurePanels.m` code in MATLAB under the `Production-Figure-Code` folder. To reproduce figure 2, for example, use the following procedure:
- Launch `Ergodic-Information-Harvesting/Production-Figure-Code/makeFigurePanels.m` using MATLAB. Note that the code has been tested with MATLAB `R2017a` and `R2018a`.
- Specify input parameters
  - Set `targetFig = 'fig2'` to select figure 2 as the target
  - Set `USE_PUBLISHED_DATASET = 1` to use the published dataset included in the repository. Alternatively, if the local simulation jobs are completed, use of `USE_PUBLISHED_DATASET = 0` will force the code to use reproduced data located at `Ergodic-Information-Harvesting/SimulationCode/SimData/`
- Run the MATLAB code

You should see a new MATLAB figure containing Figure 2 panels. PDF(s) will be saved under `Ergodic-Information-Harvesting/Production-Figure-Code/FigureOutput/fig2/`.

To reproduce all the figures, follow the same steps, but set `targetFig = 'all'`

### Benchmark Running Time
- **Benchmark on AWS**: `~22.17 hours` on AWS EC2 c5a.4xlarge instance (Ubuntu 18.04 LTS 64-bit, AMD EPYC 7R32 with boost frequency up-to 3.3GHz, 8 Core/16 Threads available under HVM)

### Additional Note for Linux and MacOS Users
#### Prevent System from Sleeping During Simulation
To prevent MacOS from sleeping in these instances, use `caffeinate` at a Terminal window running simulation jobs.
