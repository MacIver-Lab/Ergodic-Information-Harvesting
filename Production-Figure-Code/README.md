# Ergodic-Infotaxis
Ergodic-Infotaxis algorithm

## Simulation
The simulation is done in Python 3.6 using [Jupyter Notebook](http://jupyter.org/) simply IPython. Depend on the need, the simulation can be runned on:
- Local computer, which is very easy to setup but limited by the number of accessible CPU cores.
- Cloud computing virtual servers throug any popular Infrastructure as a Service (IaaS) provider, *e.g.* [Amazon Elastic Compute Cloud](https://aws.amazon.com/ec2/) or [Google Cloud Compute Engine](https://cloud.google.com/compute/). This is still easy to setup and provides a way to scale up the total number of running threads (*e.g.* Google Cloud Compute Engine allows up to 96 CPU threads per instance).

### Setting up runtime environment
The code below assumes the cloud instance runs on Linux (Ubuntu 14.04 or higher).
- Update system dependencies
```bash
sudo apt update
sudo apt upgrade
```

- Install, and then update Anaconda3
```bash
wget "https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh"
bash ./Anaconda3-5.0.1-Linux-x86_64.sh
```
```bash
source ~/.bashrc
conda update --all
```

- Configure Jupyter Notebook
  1. generate configure file
  ```bash
  jupyter notebook --generate-config
  ```
  2. append the following line to the configuration file. Note that `c.NotebookApp.port = 5123` specifies the desired port for Jupyter Notebook is `5123` and hence it needs to be allowed in the firewall for inbound traffic
  ```python
  c = get_config()
  c.NotebookApp.ip = '*'
  c.NotebookApp.open_browser = False
  c.NotebookApp.port = 5123
  ```
  
- Start Jupyter Notebook
```bash
jupyter notebook
```

### References
1. [Docs - Installing Anaconda on Linux](https://docs.anaconda.com/anaconda/install/linux)
2. [Blog - Running Jupyter Notebook on Google Cloud Platform in 15 min](https://towardsdatascience.com/running-jupyter-notebook-in-google-cloud-platform-in-15-min-61e16da34d52)
