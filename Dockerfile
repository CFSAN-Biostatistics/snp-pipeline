FROM phusion/baseimage:0.9.22
MAINTAINER Justin Payne, justin.payne@fda.hhs.gov


WORKDIR /tmp/
RUN add-apt-repository ppa:webupd8team/java \
	&& echo debconf shared/accepted-oracle-license-v1-1 select true | debconf-set-selections \
	&& echo debconf shared/accepted-oracle-license-v1-1 seen true | debconf-set-selections \
	&& apt-get update -y \
	&& apt-get install \
		build-essential \
		bzip2 \
		bowtie2 -y \
		oracle-java8-installer \
		g++ \
		gcc \
		git \
		gsl-bin \
		libgsl0-dev \
		make \
		ncurses-dev \
		python \
		python-dev \
		tabix \
		wget \
		zlib1g-dev \
	&& apt-get clean 

WORKDIR /tmp/
#install samtools
RUN	wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 -O     - | tar xj ; (cd htslib-1.3.2  ; make; make install; cd /tmp) \
 && wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 -O     - | tar xj ; (cd samtools-1.2  ; make; make install; cd /tmp) \
 && wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 -O - | tar xj ; (cd bcftools-1.3.1; make; make install; cd /tmp) 

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
	&& cp /tmp/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump.2.3.5.2 /usr/bin/fastq-dump 
	
#install picard
WORKDIR /tmp/
RUN	wget https://github.com/broadinstitute/picard/releases/download/2.17.0/picard.jar -q \
	&&	cp picard.jar /usr/bin/picard.jar

#install snp-pipeline and snp-mutator	
RUN pip install numpy snp-pipeline biopython snp-mutator 

ENV PATH "$PATH:/tmp/samtools-1.2/bin:/tmp/bcftools-1.3.1/bin"
ENV CLASSPATH "/usr/bin/VarScan.jar:/usr/bin/picard.jar"
ENV NUMCORES 1
	
#Test snp_pipeline
WORKDIR /test/
RUN cfsan_snp_pipeline data lambdaVirusInputs testLambdaVirus \
	&& cd testLambdaVirus \
	&& cfsan_snp_pipeline run -s samples reference/lambda_virus.fasta \
	&& copy_snppipeline_data.py lambdaVirusExpectedResults expectedResults \
	&& diff -q snplist.txt expectedResults/snplist.txt \
	&& diff -q snpma.fasta expectedResults/snpma.fasta \
	&& diff -q referenceSNP.fasta expectedResults/referenceSNP.fasta
	
ENTRYPOINT ["run_snp_pipeline.sh"]
CMD ["-h"]