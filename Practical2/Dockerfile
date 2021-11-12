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
COPY Practical2/environment.yml /
COPY Practical2/MOFA2_test.R /
RUN conda env create -n practical2
RUN R -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN R -e 'BiocManager::install("MOFA2")'
RUN R < MOFA2_test.R --no-save
#RUN mamba install -n practical2 -c r rstudio=1.1.456
#RUN mamba install -n practical2 -c bioconda bioconductor-mofa2=1.4.0
RUN mamba clean -a
ENV PATH /opt/conda/envs/practical2/bin:$PATH
