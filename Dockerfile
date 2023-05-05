FROM amazoncorretto:8
MAINTAINER Justin Payne, justin.payne@fda.hhs.gov


WORKDIR /tmp/
RUN yum groupinstall -y 'Development Tools' \
	&& yum install -y \
		bzip2-devel \
		gcc-c++ \
		git \
		hostname \
		make \
		ncurses-devel \
		python3 \
		python3-devel \
		tar \
		wget \
		which \
		xz-devel \
		zlib-devel \
	&& yum clean all

WORKDIR /tmp/

#Dependency versions, can be updated in the build with build_args
#https://docs.docker.com/engine/reference/builder/#using-arg-variables
ARG BCFTOOLS_VER
ENV BCFTOOLS_VER=${BCFTOOLS_VER:-1.8}
ARG BOWTIE2_VER
ENV BOWTIE2_VER=${BOWTIE2_VER:-2.5.1}
ARG HTSLIB_VER
ENV HTSLIB_VER=${HTSLIB_VER:-1.3.2}
ARG GATK_VER
ENV GATK_VER=${GATK_VER:-3.8-1-0-gf15c1c3ef}
ARG PICARD_VER
ENV PICARD_VER=${PICARD_VER:-2.27.5}
ARG SAMTOOLS_VER
ENV SAMTOOLS_VER=${SAMTOOLS_VER:-1.8}
ARG SRATOOLKIT_VER
ENV SRATOOLKIT_VER=${SRATOOLKIT_VER:-2.8.1}
ARG VARSCAN_VER
ENV VARSCAN_VER=${VARSCAN_VER:-2.3.9}

#install bowtie2
RUN wget https://github.com/BenLangmead/bowtie2/archive/v$BOWTIE2_VER.tar.gz -qO						 			- | tar xz && (cd bowtie2-$BOWTIE2_VER  && make && make install && cd /tmp)

#install samtools
RUN	wget https://github.com/samtools/htslib/releases/download/$HTSLIB_VER/htslib-$HTSLIB_VER.tar.bz2 -qO       		- | tar xj && (cd htslib-$HTSLIB_VER  	 && make && make install && cd /tmp)
RUN wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VER/samtools-$SAMTOOLS_VER.tar.bz2 -qO  	- | tar xj && (cd samtools-$SAMTOOLS_VER && make && make install && cd /tmp)
RUN wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VER/bcftools-$BCFTOOLS_VER.tar.bz2 -qO 	- | tar xj && (cd bcftools-$BCFTOOLS_VER && make && make install && cd /tmp)

#install varscan, art and sra
RUN wget https://bootstrap.pypa.io/get-pip.py -q \
	&& python3 get-pip.py

#install VARSCAN, ART, SRA Toolkit, GATK, Picard
RUN wget http://downloads.sourceforge.net/project/varscan/VarScan.v$VARSCAN_VER.jar -q \
	&& cp VarScan.v$VARSCAN_VER.jar /usr/bin/VarScan.jar 
# RUN wget https://www.niehs.nih.gov/research/resources/assets/docs/artsrcchocolatecherrycake031915linuxtgz.tgz -q \
# 	&& tar -zxf /tmp/artsrcchocolatecherrycake031915linuxtgz.tgz \
# 	&& cd /tmp/art_src_ChocolateCherryCake_Linux \
# 	&& ./configure \
# 	&& make \
# 	&& make install \
# 	&& cd /tmp/ 
RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/$SRATOOLKIT_VER/sratoolkit.$SRATOOLKIT_VER-ubuntu64.tar.gz -q \
	&& tar -zxf /tmp/sratoolkit.$SRATOOLKIT_VER-ubuntu64.tar.gz \
	&& cp /tmp/sratoolkit.$SRATOOLKIT_VER-ubuntu64/bin/fastq-dump.$SRATOOLKIT_VER /usr/bin/fastq-dump 

# https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
RUN wget --content-disposition https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-$GATK_VER.tar.bz2 -q \
	&& tar -jxf /tmp/GenomeAnalysisTK-$GATK_VER.tar.bz2 \
	&& cp /tmp/GenomeAnalysisTK-$GATK_VER/GenomeAnalysisTK.jar /usr/bin/GenomeAnalysisTK.jar 

	
RUN wget https://github.com/broadinstitute/picard/releases/download/$PICARD_VER/picard.jar -q \
	&& cp picard.jar /usr/bin/picard.jar

#install snp-pipeline and snp-mutator	
RUN pip install numpy biopython snp-mutator
WORKDIR /src/
COPY ./ /src/ 
RUN pip install .

ENV PATH "$PATH:/tmp/samtools-$SAMTOOLS_VER/bin:/tmp/bcftools-$BCFTOOLS_VER/bin:/tmp/bowtie2-$BOWTIE2_VER/bin"
ENV CLASSPATH "/usr/bin/VarScan.jar:/usr/bin/picard.jar:/usr/bin/GenomeAnalysisTK.jar"
ENV NUMCORES 4
	
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