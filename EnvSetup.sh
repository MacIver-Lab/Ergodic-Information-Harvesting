#!/bin/bash

echo "Ergodic-Infotaxis simulation environment setup script"

# Step 0 - Checkout Ergodic-Infotaxis repository and run this script
git config --global credential.helper 'cache --timeout=360000'
git config --global user.email "chenchen.bme@gmail.com"
git config --global user.name "Chen Chen"
git clone https://chenchen2015@github.com/chenchen2015/Ergodic-Infotaxis ./Ergodic-Infotaxis

# Step 1 - Update System
echo "Step 1 - Update System"
sudo apt update
sudo apt upgrade -y

# Step 2 - Install Anaconda3
echo "Step 2 - Install Anaconda3"
wget "https://repo.continuum.io/archive/Anaconda3-5.1.0-Linux-x86_64.sh"
bash ./Anaconda3-5.0.1-Linux-x86_64.sh
source ~/.bashrc
conda update --all -y
conda clean --packages -y

# Step 4 - Configure Jupyter Notebook
echo "Step 4 - Configure Jupyter Notebook"
jupyter notebook --generate-config
cp ./Ergodic-Infotaxis/jupyter_notebook_config.py .jupyter/jupyter_notebook_config.py

cd ./Ergodic-Infotaxis/Code/
jupyter notebook
