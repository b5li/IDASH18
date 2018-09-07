FROM ubuntu
MAINTAINER Baiyu Li email: baiyu@ucsd.edu

# Add the tar file of the web site 
ADD gmp /usr/src/gmp
ADD ntl /usr/src/ntl
ADD IDASH18 /usr/src/IDASH18

RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y m4
RUN apt-get install -y time

WORKDIR /usr/src
WORKDIR gmp
RUN ./configure
RUN make
RUN make install

WORKDIR /usr/src
WORKDIR ntl
WORKDIR src
RUN ./configure
RUN make
RUN make install

WORKDIR /usr/src
WORKDIR IDASH18
WORKDIR lib
RUN make clean
RUN make all

WORKDIR /usr/src
WORKDIR IDASH18
WORKDIR gwas
RUN make clean
RUN make all
RUN make testlogreg
RUN cp testlogreg /usr/bin/testlogreg


