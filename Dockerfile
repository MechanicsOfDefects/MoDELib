# MoDELib container image: the library's command-line tools, the input templates and the tutorials.
#
#   docker build -t modelib .
#   docker run -it --rm modelib
#
# Contents of the image, under /opt/MoDELib:
#   build/tools/DDomp/DDomp                                      (also on PATH as DDomp)
#   build/tools/MicrostructureGenerator/microstructureGenerator  (also on PATH)
#   Library/  python/  tutorials/
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

# Instruction set of the binaries. "native" would tie them to the CPU of the build machine.
ARG MARCH=x86-64-v2

# Three changes to the build files, made inside the image only:
#   1. Eigen: CMakeLists.txt assigns the MacPorts path.
#   2. Architecture: -march=native is replaced by -march=${MARCH}.
#   3. DDqt: removed from the tools, since it requires Qt 6 and VTK.
# The checks that follow stop the build if a pattern no longer matches.
RUN sed -i 's#set(EIGEN3_INCLUDE_DIRS /opt/local/include/eigen3)#set(EIGEN3_INCLUDE_DIRS /usr/include/eigen3)#' CMakeLists.txt \
 && sed -i "s#-march=native#-march=${MARCH}#" CMakeLists.txt \
 && sed -i '/DDqt DDqt/d' tools/CMakeLists.txt \
 && grep -q 'set(EIGEN3_INCLUDE_DIRS /usr/include/eigen3)' CMakeLists.txt \
 && ! grep -q 'march=native' CMakeLists.txt \
 && ! grep -q 'DDqt' tools/CMakeLists.txt

RUN cmake -S . -B build -DUSE_PYBIND11=OFF \
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
COPY tutorials /opt/MoDELib/tutorials

RUN ln -s /opt/MoDELib/build/tools/DDomp/DDomp /usr/local/bin/DDomp \
 && ln -s /opt/MoDELib/build/tools/MicrostructureGenerator/microstructureGenerator /usr/local/bin/microstructureGenerator

# Smoke test: two steps of the dipoleNoise tutorial. The build fails if no configuration is written.
RUN cp -r /opt/MoDELib/tutorials/dipoleNoise /opt/MoDELib/tutorials/smokeTest \
 && cd /opt/MoDELib/tutorials/smokeTest \
 && python3 generateInputFiles.py \
 && sed -i 's/^Nsteps=.*/Nsteps=2;/; s/^outputFrequency=.*/outputFrequency=1;/' inputFiles/DD.txt \
 && microstructureGenerator . \
 && DDomp . \
 && test -f evl/evl_1.txt \
 && cd / && rm -rf /opt/MoDELib/tutorials/smokeTest

LABEL org.opencontainers.image.title="MoDELib" \
      org.opencontainers.image.description="Mechanics of Defects Evolution Library: command-line tools, input templates and tutorials" \
      org.opencontainers.image.source="https://github.com/MechanicsOfDefects/MoDELib" \
      org.opencontainers.image.licenses="GPL-2.0-only"

WORKDIR /opt/MoDELib/tutorials
CMD ["bash"]
