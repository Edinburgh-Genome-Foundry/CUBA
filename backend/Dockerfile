FROM python:3.6
ENV PYTHONUNBUFFERED 1
ENV USING_DOCKER 1


RUN apt-get update
RUN apt-get install -y libgraphviz-dev graphviz

# The next line enables to build NumberJack (used by GoldenHinges)
RUN apt-get install -y python-dev swig libxml2-dev zlib1g-dev libgmp-dev

# The next line installs NCBI BLAST (for GeneBlocks, DNAWeaver, etc.)
RUN apt-get install -y ncbi-blast+

# Numberjack (requirement of GoldenHinges) is veeery slow to install,
# It has its own Docker file, line so changes in other lines don't reinstall it.
RUN pip install Numberjack
RUN mkdir /code
ADD requirements.txt /code/
WORKDIR /code
RUN pip install -r requirements.txt

# The next lines install wkhtmltopdf (for Caravagene)
ENV QT_QPA_PLATFORM offscreen
# ENV QT_QPA_FONTDIR=/usr/share/fonts
RUN apt-get install -y fontconfig libfontconfig1 libfreetype6 libx11-6 libxext6 libxrender1
RUN apt-get install -y xfonts-base xfonts-75dpi fonts-font-awesome fonts-lato
RUN apt-get install -y xvfb wkhtmltopdf
RUN apt-get install -y unzip
ENV QT_QPA_FONTDIR=/root/.fonts/
RUN wget -O "Raleway.zip" "https://fonts.google.com/download?family=Raleway" \
    && unzip -u "Raleway.zip" -d /root/.fonts \
    && rm -f "Raleway.zip"
RUN wget -O "Inconsolata.zip" "https://fonts.google.com/download?family=Inconsolata" \
    && unzip -u "Inconsolata.zip" -d /root/.fonts \
    && rm -f "Inconsolata.zip"
RUN wget -O "OpenSans.zip" "https://fonts.google.com/download?family=Open+Sans" \
    && unzip -u "OpenSans.zip" -d /root/.fonts \
    && rm -f "OpenSans.zip"
RUN wget -O "CabinSketch.zip" "https://fonts.google.com/download?family=Cabin+Sketch" \
    && unzip -u "CabinSketch.zip" -d /root/.fonts \
    && rm -f "CabinSketch.zip"
RUN wget -O "Lato.zip" "https://fonts.google.com/download?family=Lato" \
    && unzip -u "Lato.zip" -d /root/.fonts \
    && rm -f "Lato.zip"
RUN fc-cache -f

ADD . /code/

ARG CACHEBUST=2
RUN pip install --upgrade git+git://github.com/Edinburgh-Genome-Foundry/geneblocks
RUN pip install biopython==1.70
