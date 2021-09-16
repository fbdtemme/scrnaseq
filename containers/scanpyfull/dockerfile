FROM continuumio/miniconda3

RUN conda install --yes h5py
RUN conda install -c bioconda --yes loompy
RUN python -m pip install pybiomart
RUN conda install -c conda-forge --yes scanpy
