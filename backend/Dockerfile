FROM python:3.6
ENV PYTHONUNBUFFERED 1
ENV USING_DOCKER 1

RUN apt-get update
RUN apt-get install -y libgraphviz-dev graphviz

# The next line enables to build NumberJack (used by GoldenHinges)
RUN apt-get install -y python-dev swig libxml2-dev zlib1g-dev libgmp-dev

# The next line installs NCBI BLAST (for GeneBlocks, DNAWeaver, etc.)
RUN apt-get install -y ncbi-blast+

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

ARG CACHEBUST=29
RUN pip install --upgrade git+git://github.com/Edinburgh-Genome-Foundry/Geneblocks.git
# RUN pip install --upgrade git+ssh://git@git.genomefoundry.org:10022/valentin/geneblocks.git
# RUN pip install --upgrade git+ssh://git@git.genomefoundry.org:10022/valentin/goldenhinges.git


# Remove SSH keys
RUN rm -rf /root/.ssh/

ADD . /code/
