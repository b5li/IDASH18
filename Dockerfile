FROM gcc
MAINTAINER Baiyu Li email: baiyu@ucsd.edu

# Add the tar file of the web site 
ADD gmp /usr/src/gmp
ADD ntl /usr/src/ntl
ADD IDASH18 /usr/src/IDASH18

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
RUN make all

WORKDIR /usr/src
WORKDIR IDASH18
WORKDIR gwas
RUN make all
RUN make test
RUN cp foo /usr/bin/foo

# CMD ["/usr/src/IDASH18/gwas/foo","idashdata/covariates.txt","idashdata/snpMat.txt"]

