# hash:sha256:010cd9b0103b7b84416a6bba1f1c3952da5c31865f2b68ab502495489abb68e2
FROM registry.codeocean.com/codeocean/miniconda3:4.9.2-python3.8-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && add-apt-repository -y ppa:ngsolve/ngsolve \
    && add-apt-repository -y ppa:opencpn/opencpn \
    && add-apt-repository -y universe \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        libgl1-mesa-glx=21.2.6-0ubuntu0.1~20.04.1 \
        libopenmpi-dev=4.0.3-0ubuntu1 \
        netgen=6.2.1804+dfsg1-3build2 \
        openmpi-bin=4.0.3-0ubuntu1 \
        xvfb=2:1.20.13-1ubuntu1~20.04.2 \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y --channel=conda-forge --channel=ngsolve \
        matplotlib==3.5.0 \
        ngsolve==6.2.2008 \
        pyvista==0.33.2 \
        vtk==8.2.0 \
    && conda clean -ya

RUN pip install -U --no-cache-dir \
    cellorientationlinearelasticity==1.0
