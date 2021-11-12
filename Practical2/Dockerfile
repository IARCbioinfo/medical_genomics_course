################## BASE IMAGE #####################
FROM continuumio/miniconda3:4.10.3

################## METADATA #######################

LABEL base_image="continuumio/miniconda3"
LABEL version="4.10.3"
LABEL software="medical_genomics"
LABEL software.version="2021"
LABEL about.summary="Container image containing all requirements for practical 2"
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
RUN mamba env create -n practical2
RUN mamba install -n practical2 nextflow=21.04.0
RUN mamba install -n practical2 -c r rstudio=1.1.456
RUN mamba install -n practical2 -c bioconda bioconductor-mofa2=1.4.0
RUN mamba clean -a
RUN pip install mofapy2
ENV PATH /opt/conda/envs/practical2/bin:$PATH
