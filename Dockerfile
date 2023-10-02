FROM gcc:latest

RUN apt update
RUN apt -y install cmake/stable

RUN ln -s /usr/local/bin/gfortran /usr/bin/gfortran
RUN apt -y install mpich/stable

RUN git clone https://github.com/pletzer/fidibench
RUN cd fidibench && \
    mkdir build && \
    cd build && \
    cmake .. && \
    cmake --build . && \
    cmake --install . --prefix /usr/local/fidibench

#CMD ["ls", "upwind"]
