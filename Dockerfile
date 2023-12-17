# Dockerfile to automate setup of the CoreDetector multiple genome
# alignment tool into a working out-of-box Docker container.
#
# Assuming that you are running this Dockerfile from the GitHub
# working directory and that your genome data files are stored
# under the example/ subdirectory, build and run this container
# using Docker/Podman with:
#   docker build -t coredetector .
#   docker run -it -v $(pwd)/example:/example coredetector
#
# The genome aligner can then be run in the expected way in the 
# container shell:
#   ./pipeline_Minimap.sh -g example/quick_genomes.txt -o example/output -d 20
#
# Code author: Russell A. Edson, Biometry Hub
# Date last modified: 17/12/2023
FROM docker.io/ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get -y upgrade && apt-get -y install \
  build-essential \
  git \
  wget \
  openjdk-11-jdk \
  zlib1g-dev

# Download + install Minimap2-2.26 release
RUN wget "https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26.tar.bz2" \
  && tar -xjf minimap2-2.26.tar.bz2 \
  && cd minimap2-2.26 \
  && make \
  && cp minimap2 misc/paftools.js /usr/bin/ \
  && cd ..

# Download + install k8-1.0 release
RUN wget "https://github.com/attractivechaos/k8/releases/download/v1.0/k8-1.0.tar.bz2" \
  && tar -xjf k8-1.0.tar.bz2 \
  && cp k8-1.0/k8-x86_64-Linux /usr/bin/k8

# Download + install GSAlign (latest version from git repository)
RUN git clone "https://github.com/hsinnan75/GSAlign" \
  && cd GSAlign \
  && make \
  && cp bin/bwt_index bin/GSAlign /usr/bin/ \
  && cd ..

# Pull CoreDetector and install
RUN mkdir -p /CoreDetector
COPY . /CoreDetector/
RUN cp CoreDetector/MFbio.jar CoreDetector/pipeline_Minimap.sh \
  CoreDetector/pipeline_GSAlign.sh / \
  && chmod u+x pipeline_Minimap.sh pipeline_GSAlign.sh

# Helpful message for container startup
RUN echo 'echo -e "Run the CoreDetector Multiple Genome Alignment tool using:\n"\
  "  ./pipeline_Minimap.sh -g example/quick_genomes.txt -o example/output"\
  "-d 20 -n 4\nfor example. Refer to the usage guide available at"\
  "https://github.com/mfruzan/CoreDetector for the options and input formats."\
  "\nNote: make sure to write the output files to a directory under the"\
  "example/ bind-mount, so that you can access them from outside of the"\
  "Docker container afterward!\n"' > /etc/profile.d/welcome.sh

CMD [ "/bin/bash", "-l" ]
