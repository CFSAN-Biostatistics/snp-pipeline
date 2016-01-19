FROM phusion/baseimage:0.9.12
MAINTAINER Justin Payne, justin.payne@fda.hhs.gov


WORKDIR /tmp/
RUN apt-get update -y \
	&& apt-get install \
		bowtie2 -y \
		default-jre \
		g++ \
		gcc \
		gsl-bin \
		libgsl0-dev \
		make \
		python \
		python-dev \
		samtools \
		wget \
	&& apt-get clean 


#install varscan, art and sra
RUN wget https://bootstrap.pypa.io/get-pip.py -q \
	&& python get-pip.py


RUN wget http://downloads.sourceforge.net/project/varscan/VarScan.v2.3.7.jar -q \
	&& cp VarScan.v2.3.7.jar /usr/bin/VarScan.jar \
	&& wget http://www.niehs.nih.gov/research/resources/assets/docs/artsrcchocolatecherrycake031915linuxtgz.tgz -q \
	&& tar -zxf /tmp/artsrcchocolatecherrycake031915linuxtgz.tgz \
	&& cd /tmp/art_src_ChocolateCherryCake_Linux \
	&& ./configure \
	&& make \
	&& make install \
	&& cd /tmp/ \
	&& wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.3.5-2/sratoolkit.2.3.5-2-ubuntu64.tar.gz -q \
	&& tar -zxf /tmp/sratoolkit.2.3.5-2-ubuntu64.tar.gz \
	&& cp /tmp/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump.2.3.5.2 /usr/bin/fastq-dump \
	&& rm -r /tmp/*

#install snp-pipeline and snp-mutator	
RUN pip install numpy snp-pipeline biopython snp-mutator 

ENV CLASSPATH /usr/bin/VarScan.jar
ENV NUMCORES 1
	
#Test snp_pipeline
WORKDIR /test/
RUN copy_snppipeline_data.py lambdaVirusInputs testLambdaVirus \
	&& cd testLambdaVirus \
	&& run_snp_pipeline.sh -s samples reference/lambda_virus.fasta \
	&& copy_snppipeline_data.py lambdaVirusExpectedResults expectedResults \
	&& diff -q snplist.txt expectedResults/snplist.txt \
	&& diff -q snpma.fasta expectedResults/snpma.fasta \
	&& diff -q referenceSNP.fasta expectedResults/referenceSNP.fasta
	
ENTRYPOINT ["run_snp_pipeline.sh"]
CMD ["-h"]