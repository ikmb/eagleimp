FROM ubuntu:jammy
LABEL authors="Eike Wacker & Lars Wienbrandt" \
      description="Docker image containing EagleImp"

RUN apt-get update && apt-get install -y wget cmake zlib1g-dev libboost-dev libboost-filesystem-dev libboost-program-options-dev libtbb2-dev build-essential wget libbz2-dev liblzma-dev libhts-dev && rm -rf /var/lib/apt/lists/*

ADD . /opt/eagleimp

RUN cd /opt && cd eagleimp && cmake -DCMAKE_BUILD_TYPE=Release -S src -B build && cd build && make -j 

ENV PATH="/opt/eagleimp/build:/opt/eagleimp:${PATH}"
ARG PATH="/opt/eagleimp/build:/opt/eagleimp:${PATH}"