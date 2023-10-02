FROM gcc:latest

WORKDIR /app

COPY . .

RUN apt update
RUN apt -y install cmake/stable

RUN ln -s /usr/local/bin/gfortran /usr/bin/gfortran
RUN apt -y install mpich/stable

RUN git clone https://github.com/pletzer/fidibench
RUN cd fidibench && \
    mkdir build && \
    cd build && \
    cmake -D CMAKE_INSTALL_PREFIX=/usr/local .. && \
    cmake --build . && \
    cmake --install .

ENV PATH="/usr/local/bin:./:${PATH}"

