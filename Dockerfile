FROM phusion/baseimage:0.9.12
MAINTAINER Justin Payne, justin.payne@fda.hhs.gov


WORKDIR /tmp/
RUN apt-get update -y
RUN apt-get install wget default-jre gcc make g++ python python-dev samtools bowtie2 -y

RUN wget https://bootstrap.pypa.io/get-pip.py -q
RUN python get-pip.py


RUN wget http://downloads.sourceforge.net/project/varscan/VarScan.v2.3.7.jar -q
RUN cp VarScan.v2.3.7.jar /usr/bin/VarScan.jar
ENV CLASSPATH /usr/bin/VarScan.jar
ENV NUMCORES 1


RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.3.5-2/sratoolkit.2.3.5-2-ubuntu64.tar.gz -q
RUN tar -zxf /tmp/sratoolkit.2.3.5-2-ubuntu64.tar.gz
RUN cp /tmp/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump.2.3.5.2 /usr/bin/fastq-dump
RUN rm -r /tmp/*
RUN pip install snp-pipeline biopython

RUN apt-get clean