# MoDELib container image: the library's command-line tools, the input templates and the tutorials.
#
#   docker build -t modelib .
#   docker run -it --rm modelib
#
# Contents of the image, under /opt/MoDELib:
#   build/tools/DDomp/DDomp                                      (also on PATH as DDomp)
#   build/tools/MicrostructureGenerator/microstructureGenerator  (also on PATH)
#   Library/  python/  lib/  tutorials/
#   simulations/  (notebook, driver, geometries and the reference outputs of the two verification cases)
#   testsPy/clusterDynamics/  (checks of the cluster-dynamics rate equations)
# DDqt (Qt 6 and VTK) and pyMoDELib (pybind11) are not in the image.

############################
# Stage 1: build
FROM ubuntu:24.04 AS build

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake \
        libeigen3-dev libfftw3-dev libboost-dev libsuitesparse-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/MoDELib
COPY . .

# Build options (see the header of CMakeLists.txt):
#   MARCH             instruction set of the binaries. "native" would tie them to the CPU of the build machine.
#   CD_MSIZE CD_ISIZE numbers of mobile and immobile cluster-dynamics species, fixed at compile time.
#                     The defaults (4 and 8) are those of the spatialCDtest tutorial and of Zr4_Fitted.txt.
ARG MARCH=x86-64-v2
ARG CD_MSIZE=4
ARG CD_ISIZE=8

# pyMoDELib (pybind11) and DDqt (Qt 6 and VTK) are switched off.
RUN cmake -S . -B build \
        -DUSE_PYBIND11=OFF -DBUILD_DDQT=OFF \
        -DMODELIB_MARCH=${MARCH} \
        -DMODELIB_CD_MSIZE=${CD_MSIZE} -DMODELIB_CD_ISIZE=${CD_ISIZE} \
 && cmake --build build -j"$(nproc)"

############################
# Stage 2: runtime
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        libfftw3-double3 libgomp1 libcholmod5 libumfpack6 \
        python3 python3-numpy \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/MoDELib/build/tools/DDomp/DDomp \
                  /opt/MoDELib/build/tools/DDomp/DDomp
COPY --from=build /opt/MoDELib/build/tools/MicrostructureGenerator/microstructureGenerator \
                  /opt/MoDELib/build/tools/MicrostructureGenerator/microstructureGenerator
COPY Library   /opt/MoDELib/Library
COPY python    /opt/MoDELib/python
COPY lib       /opt/MoDELib/lib
COPY tutorials /opt/MoDELib/tutorials
COPY simulations /opt/MoDELib/simulations
COPY testsPy/clusterDynamics /opt/MoDELib/testsPy/clusterDynamics

RUN ln -s /opt/MoDELib/build/tools/DDomp/DDomp /usr/local/bin/DDomp \
 && ln -s /opt/MoDELib/build/tools/MicrostructureGenerator/microstructureGenerator /usr/local/bin/microstructureGenerator

# Smoke tests. Each one shortens a tutorial to two steps; the build fails if a pattern no longer matches
# or if the second configuration is not written.
#   1. dipoleNoise:   dislocation dynamics with glide-plane noise, periodic domain.
#   2. spatialCDtest: cluster dynamics alone on the finite-element mesh (mobile and immobile species).
#                     It needs the default species counts; leave it out of SMOKE_TESTS for other counts.
ARG SMOKE_TESTS="dipoleNoise spatialCDtest"
RUN set -e; \
    for tutorial in ${SMOKE_TESTS}; do \
        cp -r /opt/MoDELib/tutorials/${tutorial} /opt/MoDELib/tutorials/smokeTest; \
        cd /opt/MoDELib/tutorials/smokeTest; \
        python3 generateInputFiles.py; \
        sed -i 's/^Nsteps=[^;]*;/Nsteps=2;/; s/^outputFrequency=[^;]*;/outputFrequency=1;/' inputFiles/DefectiveCrystal.txt; \
        grep -q '^Nsteps=2;' inputFiles/DefectiveCrystal.txt; \
        grep -q '^outputFrequency=1;' inputFiles/DefectiveCrystal.txt; \
        microstructureGenerator .; \
        DDomp .; \
        test -f evl/evl_1.txt; \
        cd /; rm -rf /opt/MoDELib/tutorials/smokeTest; \
    done

LABEL org.opencontainers.image.title="MoDELib" \
      org.opencontainers.image.description="Mechanics of Defects Evolution Library: command-line tools, input templates and tutorials" \
      org.opencontainers.image.source="https://github.com/MechanicsOfDefects/MoDELib" \
      org.opencontainers.image.licenses="GPL-2.0-only"

WORKDIR /opt/MoDELib/tutorials
CMD ["bash"]
