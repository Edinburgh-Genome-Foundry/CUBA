FROM python:3.6
ENV PYTHONUNBUFFERED 1
ENV USING_DOCKER 1

RUN apt-get update
###############################################################################
# Graphviz
RUN apt-get install -y graphviz graphviz-dev libgraphviz-dev
###############################################################################
# The next line enables to build NumberJack (used by GoldenHinges)
RUN apt-get install -y python-dev swig libxml2-dev zlib1g-dev libgmp-dev
# Note: for updates using a newer OS, use the below lines:
# Ubuntu 20.04 has swig version 4 but we need version 3 to install Numberjack properly:
# RUN apt-get install -y python-dev libxml2-dev zlib1g-dev libgmp-dev
# RUN apt-get remove -y swig
# RUN apt-get install -y swig3.0
# RUN ln /usr/bin/swig3.0 /usr/bin/swig
###############################################################################
# For python-Levenshtein:
RUN apt-get install -y gcc python3-dev
RUN apt-get install -y libxslt1-dev g++
###############################################################################
# The next line installs NCBI BLAST (for GeneBlocks, DNAWeaver, etc.)
RUN apt-get install -y ncbi-blast+
###############################################################################
# Numberjack (requirement of GoldenHinges) is veeery slow to install,
# It has its own Dockerfile line so changes in other lines don't reinstall it.
RUN pip install Numberjack
# If Numberjack doesn't install from pip, use:
# Numberjack is built from source because pip doesn't install it properly:
# USER jovyan
# RUN wget https://github.com/Edinburgh-Genome-Foundry/Numberjack/archive/v1.2.0.tar.gz
# RUN tar -zxvf v1.2.0.tar.gz
# WORKDIR $HOME/Numberjack-1.2.0
# RUN python setup.py build -solver Mistral
# RUN python setup.py install
###############################################################################
# The next lines install wkhtmltopdf (for Caravagene)
ENV QT_QPA_PLATFORM offscreen
# ENV QT_QPA_FONTDIR=/usr/share/fonts
RUN apt-get install -y fontconfig libfontconfig1 libfreetype6 libx11-6 libxext6 libxrender1
RUN apt-get install -y xfonts-base xfonts-75dpi fonts-font-awesome fonts-lato
RUN apt-get install -y xvfb
RUN wget https://github.com/wkhtmltopdf/wkhtmltopdf/releases/download/0.12.4/wkhtmltox-0.12.4_linux-generic-amd64.tar.xz
RUN tar xvf wkhtmltox-0.12.4_linux-generic-amd64.tar.xz
RUN mv wkhtmltox/bin/wkhtmlto* /usr/bin/
RUN ln -nfs /usr/bin/wkhtmltopdf /usr/local/bin/wkhtmltopdf
###############################################################################
# Install fonts
RUN apt-get install -y unzip
ENV QT_QPA_FONTDIR=/root/.fonts/
RUN wget -O "Raleway.zip" "https://fonts.google.com/download?family=Raleway" \
    && unzip -u "Raleway.zip" -d /root/.fonts \
    && rm -f "Raleway.zip" \
    && rm /root/.fonts/README.txt
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

###############################################################################
# If default approach cannot find graphviz, then path is specified:
#RUN pip install pygraphviz==1.5 --install-option="--include-path=/usr/include/graphviz" --install-option="--library-path=/usr/lib/x86_64-linux-gnu/graphviz"
###############################################################################
# EGF packages:
RUN mkdir /code
ADD requirements.txt /code/
WORKDIR /code
RUN pip install -r requirements.txt

ADD . /code/
