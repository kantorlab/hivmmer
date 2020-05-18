FROM ubuntu:16.04
LABEL maintainer "Mark Howison <mhowison@ripl.org>"
LABEL repository kantorlab
LABEL image hivmmer
LABEL tag v0.2.0

RUN apt-get update -y
RUN apt-get install -y bzip2 util-linux wget

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
 && bash Miniconda3-latest-Linux-x86_64.sh -b \
 && rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH /root/miniconda3/bin:$PATH

RUN conda install -y -c kantorlab hivmmer
RUN conda clean -ay

RUN apt-get remove -y wget
RUN apt-get autoremove -y
RUN apt-get clean -y
