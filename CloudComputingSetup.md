# Use Cloud to Accelerate EIH Simulation
Cloud computing virtual servers can be used to accelerate EIH simulation by providing more available threads through any popular Infrastructure as a Service (IaaS) provider, *e.g.* [Amazon Elastic Compute Cloud](https://aws.amazon.com/ec2/) or [Google Cloud Compute Engine](https://cloud.google.com/compute/) (AWS now allows up to 128 CPU threads and Google Cloud Compute Engine allows up to 96 CPU threads per instance).

### Setting up runtime environment
The code below assumes the cloud instance runs on Linux (Ubuntu 14.04 or higher). Please follow Google or AWS's tutorial on how to start a virtual server instance and connect to it. Once you have connected to the instance, follow the steps below to setup EIH environment for simulation.

- Install Anaconda3, make sure you answer `yes` when asked whether to Anaconda to the path or not.
  ```bash
  wget "https://repo.continuum.io/archive/Anaconda3-5.1.0-Linux-x86_64.sh"
  bash ./Anaconda3-5.1.0-Linux-x86_64.sh
  ```
- Refresh shell to take effect and update conda
  ```bash
  source ~/.bashrc
  conda update -n base conda
  conda update --all -y
  ```
- Configure Jupyter Notebook
  1. generate configure file
  ```bash
  jupyter notebook --generate-config
  ```
  2. append the following line to the configuration file. Note that `c.NotebookApp.port = 5123` specifies the desired port for Jupyter Notebook is `5123` and hence it needs to be allowed in the firewall for inbound traffic, and `5123` is just an arbitrary port to work with.
  ```python
  c = get_config()
  c.NotebookApp.ip = '*'
  c.NotebookApp.open_browser = False
  c.NotebookApp.port = 5123
  ```
- Install `git` (for cloning this repo) and `tmux` (to session management)
  ```bash
  sudo apt install git tmux
  ```
  or use `yum` if you are using Amazon Linux
  ```bash
  sudo yum install git tmux
  ```
- Clone this repository (make you have `git`)
  ```bash
  git clone https://github.com/MacIver-Lab/Ergodic-Information-Harvesting ./EIH
  cd EIH
  ```

- Start a new `tmux` session and launch Jupyter Notebook
  ```bash
  tmux
  jupyter notebook&
  ```
- Copy the URL shown in the command line and then replace the IP address with the public IP address of your instance, *e.g.* `http://[Instance's public IPv4 address]:5123/?token=0519ae8b7405a2e1c558e3324d12532504c18f03a0512252`

### References
1. [Docs - Installing Anaconda on Linux](https://docs.anaconda.com/anaconda/install/linux)
2. [Blog - Running Jupyter Notebook on Google Cloud Platform in 15 min](https://towardsdatascience.com/running-jupyter-notebook-in-google-cloud-platform-in-15-min-61e16da34d52)
