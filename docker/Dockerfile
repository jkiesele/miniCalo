FROM ubuntu:18.04

RUN apt-get update
RUN apt-get install -y build-essential
RUN export DEBIAN_FRONTEND=noninteractive && apt-get install -y qt5-default

############  root part
RUN apt-get install -y libx11-dev libxpm-dev libxft-dev libxext-dev
RUN apt-get install -y cmake git
#RUN apt-get install -y python3 
RUN apt-get install -y dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev \
libxft-dev libxext-dev python libssl-dev

RUN export DEBIAN_FRONTEND=noninteractive && apt-get install -y software-properties-common
#RUN add-apt-repository ppa:rock-core/qt4 &&  apt update

RUN export DEBIAN_FRONTEND=noninteractive && apt-get install -y gfortran libpcre3-dev \
xlibmesa-glu-dev libglew1.5-dev libftgl-dev \
libmysqlclient-dev libfftw3-dev libcfitsio-dev \
graphviz-dev libavahi-compat-libdnssd-dev \
libldap2-dev python-dev libxml2-dev libkrb5-dev \
libgsl0-dev libqt4-dev

#RUN mkdir /temproot && \
#    cd /temproot && \
#    git clone http://github.com/root-project/root.git && \
#    cd root && \
#    git checkout -b v6-22-00-patches
#RUN cd /temproot && \
#    mkdir build && \
#    cd build && \
#    cmake ../root 
#RUN cd /temproot/build && \
#    make -j5 && \
#    make install && \
#    rm -rf /temproot


#ensure import ROOT
#ENV PYTHONPATH="/usr/local/lib:${PYTHONPATH}"
RUN ldconfig

RUN apt-get install -y wget apt libxmu-dev

RUN mkdir /tempgeant &&\
    cd /tempgeant &&\
    wget https://geant4-data.web.cern.ch/releases/geant4.10.07.tar.gz &&\
    tar xfz *.gz 


RUN ls /tempgeant/ &&\
    cd /tempgeant/geant4.10.07 &&\
    mkdir ../build &&\
    cd ../build &&\
    cmake -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_OPENGL_X11=ON  -DGEANT4_USE_QT=ON  ../geant4.10.07/

RUN cd /tempgeant/build &&\
    make -j5 &&\
    make install &&\
    rm -rf /tempgeant


RUN mkdir -p /usr/local/G4Rhadrons_with_Regge_model
COPY G4Rhadrons_with_Regge_model /usr/local/G4Rhadrons_with_Regge_model


SHELL ["/bin/bash", "-c"]

RUN apt-get install -y libxi-dev

RUN cd /usr/local/share/Geant4-10.7.0/geant4make &&\
    . geant4make.sh &&\
    cd /usr/local/G4Rhadrons_with_Regge_model &&\
    make 

RUN cp /root/geant4_workdir/tmp/Linux-g++/exampleRhad/libexampleRhad.so /usr/local/G4Rhadrons_with_Regge_model/

ENV LD_LIBRARY_PATH="/usr/local/G4Rhadrons_with_Regge_model/:${LD_LIBRARY_PATH}"
ENV G4RHADRONS="/usr/local/G4Rhadrons_with_Regge_model"


RUN mkdir -p /usr/local/ &&\
    cd /usr/local/ &&\
    wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8303.tgz &&\
    tar xzf pythia8303.tgz &&\
    mv pythia8303 pythia8

RUN cd /usr/local/pythia8 &&\
    ./configure  &&\
    make -j5
RUN rm -f /usr/local/pythia8/pythia8303.tgz     

ENV PYTHIA8_HOME="/usr/local/pythia8"
ENV PYTHIA8__HOME="/usr/local/pythia8"
ENV LD_LIBRARY_PATH="/usr/local/pythia8/lib:${LD_LIBRARY_PATH}"


RUN apt-get install -y libfastjet-dev

ENV G4DATAINSTALL=/usr/local/share/Geant4-10.7.0/data
ENV G4LEVELGAMMADATA=/usr/local/share/Geant4-10.7.0/data/PhotonEvaporation5.7
ENV G4NEUTRONXSDATA=/usr/local/share/Geant4-10.7.0/data/G4NEUTRONXS1.4
ENV G4LEDATA=/usr/local/share/Geant4-10.7.0/data/G4EMLOW7.13
ENV G4SAIDXSDATA=/usr/local/share/Geant4-10.7.0/data/G4SAIDDATA2.0
ENV G4RADIOACTIVEDATA=/usr/local/share/Geant4-10.7.0/data/RadioactiveDecay5.6
ENV G4NEUTRONHPDATA=/usr/local/share/Geant4-10.7.0/data/G4NDL4.6
ENV G4ENSDFSTATEDATA=/usr/local/share/Geant4-10.7.0/data/G4ENSDFSTATE2.3
ENV G4PARTICLEXSDATA=/usr/local/share/Geant4-10.7.0/data/G4PARTICLEXS3.1

RUN apt-get install -y nano x11-apps
