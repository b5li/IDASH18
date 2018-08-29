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

2. Run the docker build command to create an image called "idash18", which can be uploaded to Docker Hub:

```
docker build --rm -t idash18 .

docker tag idash18 baiyuli/idash18
```

# To run the docker image:

```
docker run -it --rm \
 -v "$PWD"/IDASH18/gwas/idashdata/covariates.txt:/data/covariates.txt \
 -v "$PWD"/IDASH18/gwas/idashdata/snpMat.txt:/data/snpMat.txt \
 idash18 foo /data/covariates.txt /data/snpMat.txt
```
