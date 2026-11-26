ARG UBUNTU_VERSION=24.04

FROM ubuntu:${UBUNTU_VERSION} AS build

ARG MARCH=x86-64-v2
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential cmake git \
        zlib1g-dev libboost-filesystem-dev libboost-program-options-dev libtbb-dev libhts-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src

COPY . .

RUN cmake -S src -B build -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG -march=${MARCH} -mtune=generic" \
    && cmake --build build -j"$(nproc)"

FROM ubuntu:${UBUNTU_VERSION}
LABEL org.opencontainers.image.source="https://github.com/ikmb/eagleimp" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        libboost-program-options1.83.0 libgomp1 libhts3t64 libtbb12 bcftools gawk \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/build/eagleimp /src/launch_eagleimp /usr/local/bin/
COPY --from=build /src/*.awk /opt/eagleimp/

RUN sed -i -E '1,/^###$/ s/^([A-Z0-9]+)=(.*)$/\1=${\1-\2}/' /usr/local/bin/launch_eagleimp
ENV REPOPATH=/opt/eagleimp
WORKDIR /data
CMD ["eagleimp", "-v"]
