# EIH Simulation Docker image recipe
# Chen Chen
# chenchen.bme@gmail.com

# based on miniconda3
FROM continuumio/miniconda3

# get gcc and g++ for compiling Cython
RUN apt-get update --fix-missing && \
    apt-get install -y gcc g++ && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN export conda=/opt/conda/bin/conda && \
    export pip=/opt/conda/bin/pip && \
    export python3=/opt/conda/bin/python && \
    export python=python3

# install python dependencies
RUN conda install -c conda-forge -y -q tqdm cython numba scipy=1.0.1 numpy=1.16.1 && \
    conda clean --all -y -q

# create binding point for EIH repository
RUN mkdir /EIH
VOLUME /EIH

# entry point is bash
CMD [ "/bin/bash" ]
