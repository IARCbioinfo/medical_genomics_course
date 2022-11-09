################## BASE IMAGE #####################
FROM continuumio/miniconda3:4.10.3

################## METADATA #######################

LABEL base_image="continuumio/miniconda3"
LABEL version="4.10.3"
LABEL software="medical_genomics"
LABEL software.version="2022"
LABEL about.summary="Container image containing all requirements for practical 1"
LABEL about.home="http://github.com/IARCbioinfo/medical_genomics_course"
LABEL about.documentation="http://github.com/IARCbioinfo/medical_genomics_course/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/medical_genomics_course/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER **nalcala** <**alcalan@iarc.fr**>

################## INSTALLATION ######################
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda install -c conda-forge mamba
COPY Practical1/environment.yml /
RUN mamba env create -n practical1 
RUN mamba install -n practical1 nextflow=22.10.0
RUN mamba install -n practical1 salmon=1.9.0
RUN mamba install -n practical1 fastqc=0.11.9
RUN mamba install -n practical1 multiqc=1.13
RUN mamba install -n practical1 samtools=1.16.1
RUN mamba install -n practical1 tbb=2021.6.0
RUN mamba clean -a
ENV PATH /opt/conda/envs/practical1/bin:$PATH
