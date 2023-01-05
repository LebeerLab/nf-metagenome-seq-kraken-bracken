FROM condaforge/mambaforge

ENV TDYAMPVER="v0.2.2"

LABEL org.opencontainers.image.authors="tim.van.rillaer@hotmail.com"

# Install all conda cloud available tools
COPY conda.yml .
RUN mamba env update -n root -f conda.yml

RUN apt-get update && apt-get install --yes --no-install-recommends

# Install tidyamplicons
RUN R -e 'remotes::install_github("SWittouck/tidyamplicons", ref = Sys.getenv("TDYAMPVER"))'
