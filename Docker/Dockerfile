ARG SCAMP_version=v2.6.7
ARG SWARP_version=2.38.0
ARG Sextractor_version=2.25.0
ARG PSFEX_version=3.21.1

# Building image from
FROM ubuntu:18.04

LABEL maintainer="David Corre <david.corre.fr@gmail.com>"

# Create a variable corresponding to the name of the new user that will be created inside the Docker.
ENV USR newuser

# Set up working directory
WORKDIR /home

# Install updates to base image
RUN \
  apt-get update -y \
  && apt-get install -y 

# Install requirements
RUN \
  apt-get install unzip wget curl git build-essential gfortran autoconf libtool libfftw3-dev libplplot-dev vim -y \ 
  && apt-get install libatlas3-base libatlas-base-dev libcurl4-openssl-dev -y \
  && apt-get autoremove -y \
  && apt-get clean -y 

# Create new user
RUN useradd -ms /bin/bash ${USR}

# Download cfitsio
RUN \ 
   curl -OL -m 3600 http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio_latest.tar.gz \
   && tar zxf cfitsio_latest.tar.gz \
   && rm -f cfitsio_latest.tar.gz

# Build cfitsio
RUN \
   cd cfitsio-* \
   && ./configure \
   && make \
   && make fpack \
   && make funpack \
   && make install \
   && cp fpack funpack /usr/local/bin/ \
   && cp fitsio*.h longnam.h /usr/local/include/ \
   && cp libcfitsio.* /usr/lib/

# Create astromatic folder
RUN \
   mkdir astromatic

# Download and build cdsclient used by SCAMP
RUN \
   curl -OL -m 3600 http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient.tar.gz \
   && tar xfz cdsclient.tar.gz \
   && rm -f cdsclient.tar.gz \
   && cd cdsclient-* \
   && ./configure \
   && make \
   && make install

# Downloads and build astromatic softwares
ARG SCAMP_version
ARG SWARP_version
ARG Sextractor_version
ARG PSFEX_version
RUN \
   curl -m 3600 -L -o astromatic/psfex_$PSFEX_version.zip https://github.com/astromatic/psfex/archive/$PSFEX_version.zip \
   && unzip -q astromatic/psfex_$PSFEX_version.zip -d astromatic/ \
   && rm -f astromatic/psfex_$PSFEX_version.zip \
   && cd astromatic/psfex-* \
   && sh autogen.sh \
   && ./configure \
   && make \
   && make install

RUN \
   curl -m 7200 -L -o astromatic/scamp_$SCAMP_version.zip https://github.com/astromatic/scamp/archive/$SCAMP_version.zip \
   && unzip -q astromatic/scamp_$SCAMP_version.zip -d astromatic/ \
   && rm -f astromatic/scamp_$SCAMP_version.zip \
   && cd astromatic/scamp-* \
   && sh autogen.sh \
   && ./configure \
   && make \
   && make install

RUN \
   curl -m 3600 -L -o astromatic/sextractor_$Sextractor_version.zip https://github.com/astromatic/sextractor/archive/$Sextractor_version.zip \
   && unzip -q astromatic/sextractor_$Sextractor_version.zip -d astromatic/ \
   && rm -f astromatic/sextractor_$Sextractor_version.zip \
   && cd astromatic/sextractor-* \
   && sh autogen.sh \
   && ./configure \
   && make \
   && make install

RUN \
   curl -m 7200 -L -o astromatic/swarp_$SWARP_version.zip https://github.com/astromatic/swarp/archive/$SWARP_version.zip \
   && unzip -q astromatic/swarp_$SWARP_version.zip -d astromatic/ \
   && rm -f astromatic/swarp_$SWARP_version.zip \ 
   && cd astromatic/swarp-* \
   && ./configure \
   && make \
   && make install

# Clone and build hotpants
RUN \
   git clone https://github.com/acbecker/hotpants.git \
   && cd hotpants \
   && make \
   && cp hotpants /usr/local/bin/

# switch ownership
# (all commands are root until a USER command changes that)
USER ${USR}

# Set the working directory
WORKDIR /home/${USR}

# change ownership from root to USR:
RUN chown -R ${USR}:${USR}  /home/${USR}

# Create directory to link on volume with host machine
RUN \
   mkdir /home/${USR}/gmadet/

# Install python libraries through miniconda
ENV PATH="/home/${USR}/miniconda3/bin:$PATH"
RUN \
   curl -m 7200 -OL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \ 
   && bash Miniconda3-latest-Linux-x86_64.sh -b \
   && rm -f Miniconda3-latest-Linux-x86_64.sh 

# Install python libraries
RUN \
   conda install python=3 numpy scipy matplotlib astropy pandas shapely requests h5py sphinx sphinx_rtd_theme scikit-image 

RUN \ 
   pip install lacosmic hjson voevent-parse xmltodict astroML regions photutils keras keras-vis opencv-python \
   && pip install --pre astroquery

# Change working directory to gmadet
WORKDIR /home/${USR}/gmadet

# define entrypoint which is the default executable at launch
ENTRYPOINT ["bash"]
