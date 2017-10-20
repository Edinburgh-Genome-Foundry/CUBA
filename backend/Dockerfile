FROM python:3.6
ENV PYTHONUNBUFFERED 1
ENV USING_DOCKER 1


RUN apt-get update
RUN apt-get install -y libgraphviz-dev graphviz

# The next line enables to build NumberJack (used by GoldenHinges)
RUN apt-get install -y python-dev swig libxml2-dev zlib1g-dev libgmp-dev

# The next line installs NCBI BLAST (for GeneBlocks, DNAWeaver, etc.)
RUN apt-get install -y ncbi-blast+

# The next lines install wkhtmltopdf (for Caravagene)
ENV QT_QPA_PLATFORM offscreen
# ENV QT_QPA_FONTDIR=/usr/share/fonts
RUN apt-get install -y fontconfig libfontconfig1 libfreetype6 libx11-6 libxext6 libxrender1 xfonts-base xfonts-75dpi
RUN apt-get install -y xvfb wkhtmltopdf
RUN apt-get install -y unzip
ENV QT_QPA_FONTDIR=/root/.fonts/
RUN wget -O "Raleway.zip" "https://fonts.google.com/download?family=Raleway" \
    && unzip -u "Raleway.zip" -d /root/.fonts \
    && rm -f "Raleway.zip"
RUN wget -O "OpenSans.zip" "https://fonts.google.com/download?family=Open+Sans" \
    && unzip -u "OpenSans.zip" -d /root/.fonts \
    && rm -f "OpenSans.zip"
RUN wget -O "CabinSketch.zip" "https://fonts.google.com/download?family=Cabin+Sketch" \
    && unzip -u "CabinSketch.zip" -d /root/.fonts \
    && rm -f "CabinSketch.zip"
RUN fc-cache -f


RUN mkdir /root/.ssh && chmod 700 /root/.ssh
RUN ssh-keyscan -H -p 10022 git.genomefoundry.org >> /root/.ssh/known_hosts
RUN ssh-keyscan -H -p 10022 178.62.46.213 >> /root/.ssh/known_hosts

# Add a deploy key
COPY ./.ssh-deploy-keys/dna_cauldron_deploy /root/.ssh/id_rsa
COPY ./.ssh-deploy-keys/dna_cauldron_deploy.pub /root/.ssh/id_rsa.pub
RUN chmod 600 /root/.ssh/id_rsa


RUN mkdir /code
ADD requirements.txt /code/
WORKDIR /code
RUN pip install -r requirements.txt

# Remove SSH keys
RUN rm -rf /root/.ssh/

ADD . /code/

ARG CACHEBUST=51
RUN pip install --upgrade git+git://github.com/Edinburgh-Genome-Foundry/BandWitch.git
RUN pip install --upgrade git+git://github.com/Edinburgh-Genome-Foundry/Plateo.git
