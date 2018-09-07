# To create a Docker image:

1. Download gmp and ntl libraries, unpack to directories "gmp" and "ntl" relative to the root of this repo


```
wget https://ftp.gnu.org/gnu/gmp/gmp-6.1.2.tar.bz2
tar xvvjf gmp-6.1.2.tar.bz2
mv gmp-6.1.2 gmp
wget http://www.shoup.net/ntl/ntl-11.3.0.tar.gz
tar xvvzf ntl-11.3.0.tar.gz
mv ntl-11.3.0 ntl
```

2. Run the docker build command to create an image called "idash18linreg", which can be uploaded to Docker Hub:

```
docker build --rm -t idash18linreg .
docker tag idash18linreg idash18ucsd/helinreg
docker push idash18ucsd/helinreg
```

# To run the docker image:

1. Install the docker software

..* on Ubuntu/Debian, do `apt-get install docker`
..* on Mac, the docker software can be downloaded here <https://store.docker.com/editions/community/docker-ce-desktop-mac>

2. Now we can pull the uploaded docker image (on linux we must use `sudo docker ...` but on Mac we can just do `docker ...` with regular user account):

..* To pull and run the docker image for linear regression:
```
sudo docker pull idash18ucsd/helinreg

sudo docker run -it --rm idash18ucsd/helinreg
```

..* To pull and run the docker image for logistic regression:
```
sudo docker pull idash18ucsd/helogreg

sudo docker run -it --rm idash18ucsd/helogreg
```

3. Then we shall see a command line prompt of the docker environment:

```
root@...:/usr/src/IDASH18/gwas# 
```

We can run the test program inside the docker environment like this:

```
root@...:/usr/src/IDASH18/gwas# ./testlinreg idashdata/covariates.txt idashdata/snpMat.txt
```

We can also run the test program from outside of the docker environment, e.g., in a terminal:

```
sudo docker run -it --rm \
 -v "$PWD"/IDASH18/gwas/idashdata/covariates.txt:/data/covariates.txt \
 -v "$PWD"/IDASH18/gwas/idashdata/snpMat.txt:/data/snpMat.txt \
 idash18 testlinreg /data/covariates.txt /data/snpMat.txt
```

The "-v" option mounts a file to a path inside the docker environment. 
The test programs "testlinreg" and "testlogreg" are copied to `/usr/bin/`, 
so we can call it just as "testlinreg" and "testlogreg".

