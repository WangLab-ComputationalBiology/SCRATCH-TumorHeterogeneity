# Use a specific version of Ubuntu as the base image
FROM --platform=linux/x86_64 rocker/verse:latest

# Set the working directory inside the container
WORKDIR /opt

# Timezone settings
ENV TZ=US/Central
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone
# pass your PAT at build time so remotes::install_github can auth
ARG GITHUB_PAT
ENV GITHUB_PAT=${GITHUB_PAT}

# Install system dependencies
RUN apt-get update && apt-get install -y \
    software-properties-common \
    dirmngr \
    gnupg \
    apt-transport-https \
    ca-certificates \
    wget \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    default-jre \
    libgfortran5 \
    liblapack-dev \
    libopenblas-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff5-dev \
    zlib1g-dev \
    libxt-dev


# Updating quarto to Quarto v1.4.553
RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.553/quarto-1.4.553-linux-amd64.deb -O quarto-1.4.553-linux-amd64.deb
RUN dpkg -i quarto-1.4.553-linux-amd64.deb

# # Install remotes package
# RUN R -e "install.packages('remotes')"

# Install Python3
# RUN apt-get install -y \
#     python3 \
#     python3-pip
# RUN apt-get update && apt-get install -y python3 python3-pip python3-venv
RUN apt-get update && apt-get install -y python3 python3-pip python3-venv python3-dev build-essential

# Install core R packages
RUN Rscript -e "install.packages(c('R.utils','rmarkdown','devtools','tidyverse','readr', 'dplyr', 'ggplot2', 'cowplot', 'remotes', 'BiocManager','reticulate', 'HGNChelper'), repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages(c('leiden', 'RColorBrewer', 'reshape2', 'scales', 'NMF', 'colorspace', 'tibble', 'data.table', 'stringr', 'Matrix', 'bigmemory', 'patchwork', 'pheatmap'), repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('viridis',dependencies = TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('bigmemory',dependencies = TRUE)"
RUN Rscript -e "install.packages('doMC', dependencies = TRUE, repos='http://R-Forge.R-project.org')"
RUN Rscript -e "install.packages('optparse', dependencies = TRUE)"
RUN Rscript -e "install.packages('pheatmap', dependencies = TRUE)"


RUN Rscript -e "BiocManager::install(c('S4Vectors','DelayedMatrixStats','BiocGenerics','Biobase', 'SummarizedExperiment', 'AnnotationDbi', 'org.Hs.eg.db'), ask=FALSE, update=TRUE)"
RUN Rscript -e "BiocManager::install(c( \
    'HDF5Array','rhdf5','rhdf5lib', \
    'SingleCellExperiment', \
    'GOSemSim','MatrixGenerics','treeio','DOSE','ggtree','enrichplot', \
    'clusterProfiler','DirichletMultinomial','rtracklayer','GenomicFeatures', \
    'BSgenome','ensembldb','TFBSTools', \
    'BSgenome.Hsapiens.UCSC.hg38','EnsDb.Hsapiens.v86', 'Homo.sapiens'), \
  ask=FALSE, update=FALSE )"

# Setting repository URL
ARG R_REPO="http://cran.us.r-project.org"

# # Caching R-lib on the building process
# RUN Rscript -e "install.packages(${R_DEPS}, Ncpus = 8, repos = '${R_REPO}', clean = TRUE)"
# RUN Rscript -e "install.packages(${WEB_DEPS}, Ncpus = 8, repos = '${R_REPO}', clean = TRUE)"

# # Install BiocManager
# RUN Rscript -e "BiocManager::install(${R_BIOC_DEPS})"
# RUN Rscript -e 'remotes::install_github("ctlab/fgsea")'
RUN sed -i 's/-Werror=format-security//g' /usr/local/lib/R/etc/Makeconf
RUN Rscript -e "remotes::install_github('jlaffy/scalop', dependencies=TRUE)"
RUN Rscript -e "install.packages('NMF')"


RUN Rscript -e 'BiocManager::install("readr", dependencies = TRUE)'
# RUN Rscript -e 'BiocManager::install("dplyr", dependencies = TRUE)'
# RUN Rscript -e 'BiocManager::install("ggplot2", dependencies = TRUE)'
# RUN Rscript -e 'BiocManager::install("Seurat", dependencies = TRUE)'
# RUN Rscript -e 'BiocManager::install("DT", dependencies = TRUE)'
# RUN Rscript -e 'BiocManager::install("SingleCellExperiment", dependencies = TRUE)'
# RUN Rscript -e 'BiocManager::install("scDblFinder", dependencies = TRUE, force = TRUE)'
# RUN Rscript -e 'BiocManager::install("lpsymphony", dependencies = TRUE, force = TRUE)'
# RUN Rscript -e 'BiocManager::install("IHW", dependencies = TRUE, force = TRUE)'
# RUN Rscript -e 'BiocManager::install("scp", dependencies = TRUE, force = TRUE)'
# RUN Rscript -e 'BiocManager::install(c("DOSE", "enrichplot", "clusterProfiler"), force = TRUE)'

# Install Seurat Wrappers
RUN wget https://github.com/satijalab/seurat/archive/refs/heads/seurat5.zip -O /opt/seurat-v5.zip
RUN wget https://github.com/satijalab/seurat-data/archive/refs/heads/seurat5.zip -O /opt/seurat-data.zip
RUN wget https://github.com/satijalab/seurat-wrappers/archive/refs/heads/seurat5.zip -O /opt/seurat-wrappers.zip

RUN Rscript -e "devtools::install_local('/opt/seurat-v5.zip')"
RUN Rscript -e "devtools::install_local('/opt/seurat-data.zip')"
RUN Rscript -e "devtools::install_local('/opt/seurat-wrappers.zip')"
RUN Rscript -e "devtools::install_github('BlakeRMills/MetBrewer')"
RUN Rscript -e "devtools::install_github('sjmgarnier/viridis')"



# # Install SCP package from GitHub
# # RUN Rscript -e "remotes::install_github('bnprks/BPCells/r')"
# RUN Rscript -e "devtools::install_github('PaulingLiu/ROGUE', dependencies = TRUE, force = TRUE)"
# # RUN Rscript -e "devtools::install_github('zhanghao-njmu/SCP', dependencies = TRUE, force = TRUE)"
# RUN Rscript -e "remotes::install_github('zhanghao-njmu/SCP', upgrade = 'always', dependencies = TRUE, force = TRUE)"
# RUN Rscript -e "remotes::install_github('cellgeni/sceasy', upgrade = 'always', dependencies = TRUE, force = TRUE)"


# # Download the Miniconda installer
# RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
#     chmod +x /tmp/miniconda.sh && \
#     /tmp/miniconda.sh -b -p /opt/miniconda && \
#     rm /tmp/miniconda.sh

# # Update PATH environment variable
# ENV PATH=/opt/miniconda/bin:$PATH


# # Install R packages
# RUN Rscript -e 'install.packages("remotes")' && \
#     Rscript -e 'remotes::install_github("zhanghao-njmu/SCP", upgrade = "always", force = TRUE, quiet = TRUE)' \
#     Rscript -e 'SCP::PrepareEnv( \
#             miniconda_repo = "https://mirrors.bfsu.edu.cn/anaconda/miniconda", \
#             pip_options = "-i https://pypi.tuna.tsinghua.edu.cn/simple")'

# Set the conda binary path and prepare the SCP environment
# RUN Rscript -e 'options(reticulate.conda_binary = "/opt/miniconda/bin/conda"); SCP::PrepareEnv(force = TRUE)'


# RUN Rscript -e 'renv::activate()'
# RUN wget https://github.com/zhanghao-njmu/SCP/archive/refs/heads/main.zip -O /opt/SCP.zip
# RUN unzip -o /opt/SCP.zip -d /opt/SCP
# RUN Rscript -e "devtools::install_local('/opt/SCP/SCP-main')"


# RUN Rscript -e "devtools::install_local('/opt/SCP.zip')"
#  RUN Rscript -e 'devtools::install_github("zhanghao-njmu/SCP")'
# RUN Rscript -e 'remotes::install_github("zhanghao-njmu/SCP", dependencies = TRUE, force = TRUE)'


# Install packages on Github
# RUN Rscript -e "devtools::install_github(${DEV_DEPS})"


# # Create and activate virtual environment
# RUN python3 -m venv /opt/venv
# ENV PATH="/opt/venv/bin:$PATH"

# # # Upgrade pip and install Python packages in venv
# # RUN pip install --upgrade pip && \
# #     pip install numpy pandas scikit-learn matplotlib seaborn jupyter jupyter-cache papermill

# # Create and use /opt/venv for everything Python
# RUN python3 -m venv /opt/venv
# ENV PATH="/opt/venv/bin:$PATH"
# ENV RETICULATE_PYTHON="/opt/venv/bin/python" \
#     PYTHONNOUSERSITE=1


# # Create and activate a virtual environment before installing Python packages
# RUN python3 -m venv /opt/venv \
#     && /opt/venv/bin/pip install --no-cache-dir numpy pandas scikit-learn matplotlib seaborn jupyter jupyter-cache papermill anndata scanpy scipy session_info scSpectra metatime python-igraph leidenalg && \
#     python - <<'PY'
#     && ln -s /opt/venv/bin/python /usr/local/bin/python \
#     && ln -s /opt/venv/bin/pip /usr/local/bin/pip

# import importlib, sys
# for m in ("igraph","leidenalg"):
#     importlib.import_module(m)
# print("OK: igraph/leidenalg present in", sys.executable)
# PY

#     # Create and activate virtual environment
# RUN python -m venv /opt/venv
# ENV PATH="/opt/venv/bin:$PATH"

# Build tools so pip can compile wheels if needed (arm64, etc.)
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential python3-dev \
  && rm -rf /var/lib/apt/lists/*

# Create venv and install Python deps (including python-igraph + leidenalg)
RUN python3 -m venv /opt/venv \
  && /opt/venv/bin/pip install --upgrade pip setuptools wheel \
  && /opt/venv/bin/pip install --no-cache-dir \
       numpy pandas scipy scikit-learn matplotlib seaborn \
       jupyter jupyter-cache papermill \
       anndata scanpy session_info scSpectra metatime celltypist \
       python-igraph==0.11.* leidenalg==0.10.* \
  && /opt/venv/bin/python - <<'PY'
import importlib, sys
for m in ("igraph","leidenalg"):
    importlib.import_module(m)
print("OK: igraph/leidenalg present in", sys.executable)
PY

# Make the venv the default for PATH and for R/reticulate
ENV PATH="/opt/venv/bin:${PATH}"
ENV RETICULATE_PYTHON=/opt/venv/bin/python


# Setting celltypist variable
# ENV CELLTYPIST_FOLDER=/opt/celltypist

# # Installing celltypist models
# COPY setup.py /opt/
# RUN python3 /opt/setup.py

# # Download and install SCP manually
# RUN wget https://github.com/zhanghao-njmu/SCP/archive/refs/heads/main.zip -O /opt/SCP.zip
# RUN unzip /opt/SCP.zip -d /opt/SCP

# # Force SCP to use system Python (modify reticulate before installing SCP)
# RUN Rscript -e "install.packages('reticulate', repos='http://cran.us.r-project.org')" && \
#     Rscript -e "library(reticulate); use_python('/usr/bin/python3', required=TRUE); options(reticulate.conda_binary=NULL, SCP_env_name=NULL)" && \
#     Rscript -e "devtools::install_local('/opt/SCP/SCP-main')"
    
# Install Python packages for data science
# RUN python3 -m pip install --no-cache-dir numpy pandas scikit-learn matplotlib seaborn jupyter
# RUN python3 -m pip install --no-cache-dir jupyter-cache
# RUN python3 -m pip install --no-cache-dir papermill

# # Install Python packages for data science
# RUN python3 -m pip install --no-cache-dir numpy pandas scikit-learn matplotlib seaborn jupyter
# RUN python3 -m pip install --no-cache-dir jupyter-cache
# RUN python3 -m pip install --no-cache-dir papermill


# Additional packages
# RUN apt-get install -y libhdf5-dev
# RUN Rscript -e "install.packages('hdf5r')"
# Additional packages
RUN apt-get update && apt-get install -y libhdf5-dev && rm -rf /var/lib/apt/lists/*
RUN Rscript -e "install.packages('hdf5r')"


# Java + Fortran 
RUN apt-get update && apt-get install -y default-jre libgfortran5

# JAGS
RUN apt-get install -y jags

# Install infercnv (fallback if BiocManager fails)
# RUN Rscript -e "remotes::install_github('broadinstitute/inferCNV', dependencies = TRUE, upgrade = 'never')"

# Install SCEVAN
# RUN Rscript -e "devtools::install_github('AntonioDeFalco/SCEVAN', dependencies = TRUE, upgrade = 'never')"


# Install copykat
# RUN Rscript -e "devtools::install_github('navinlabcode/copykat')"

# Install annotables
# RUN Rscript -e "devtools::install_github('stephenturner/annotables')"
# RUN Rscript -e "install.packages('rJava', repos = 'http://cran.us.r-project.org')"
# RUN Rscript -e "devtools::install_github('miccec/yaGST', dependencies = TRUE, upgrade = 'never')"

# RUN Rscript -e "devtools::install_github('AntonioDeFalco/SCEVAN', dependencies = TRUE, upgrade = 'never')"
# # install SCP
# RUN Rscript -e "remotes::install_github( \
#   'zhanghao-njmu/SCP', \
#   dependencies=TRUE, \
#   upgrade='always', \
#   auth_token = Sys.getenv('GITHUB_PAT'))"


RUN apt-get update && \
    apt-get install -y --no-install-recommends \
       libgsl-dev \
    && rm -rf /var/lib/apt/lists/*  
  # Install DirichletMultinomial + TFBSTools (and all of their R & Bioc deps)
# RUN Rscript -e "BiocManager::install('DirichletMultinomial', ask = FALSE, update      = FALSE, dependencies = TRUE)"
# RUN Rscript -e "BiocManager::install('TFBSTools', ask = FALSE, update      = FALSE, dependencies = TRUE)"
# RUN Rscript -e "install.packages(c('DirichletMultinomial','TFBSTools'), dependencies = TRUE, repos = BiocManager::repositories())"
# # install Azimuth
# RUN Rscript -e "remotes::install_github( \
#   'satijalab/azimuth', ref = 'master', \
#   dependencies=TRUE, \
#   upgrade='always', \
#   auth_token = Sys.getenv('GITHUB_PAT'))"


#RUN Rscript -e "remotes::install_version('Matrix', version = '1.6-1')"
#RUN Rscript -e "install.packages('SeuratObject')"
#RUN Rscript -e "install.packages('scCustomize')"

# Cleaning apt-get cache
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

# 5) validate loads at build time
# RUN Rscript -e "library(SingleCellExperiment); library(SCP)"

# Command to run on container start
CMD ["bash"]

