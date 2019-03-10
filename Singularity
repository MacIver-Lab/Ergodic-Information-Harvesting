BootStrap: docker
From: continuumio/miniconda3

%help
    This is a Singularity container for EIH simulations.
    It builds a minimal Python 3 environment running on Linux with all the libraries required for EIH simulations.

%labels
    Maintainer Chen Chen (chenchen.bme@gmail.com)
    Version v0.2
   
%environment
     export conda=/opt/conda/bin/conda
     export pip=/opt/conda/bin/pip
     export python3=/opt/conda/bin/python
     export python=python3

%post 
     # update system and install gcc
     apt-get update
     apt-get install gcc g++ -y
     # update conda
     /opt/conda/bin/conda update --all -y --quiet
     # Update pip
     /opt/conda/bin/pip install -U pip -q
     # Install dependencies
     /opt/conda/bin/conda install -c conda-forge -y -q tqdm cython numba scipy=1.0.1 numpy=1.16.1 
     # Clean up
     /opt/conda/bin/conda clean --all -y --quiet
     apt-get autoremove -y
     apt-get clean
     # create bind points for HPCC environment
     mkdir -p /EIH

%test  
     echo "Testing python..."
     /opt/conda/bin/python -V
