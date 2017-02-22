FROM python:3.6
ENV PYTHONUNBUFFERED 1
ENV USING_DOCKER 1

RUN mkdir ~/.ssh

RUN apt-get update && apt-get install -y redis-server

RUN mkdir /code
WORKDIR /code
ADD requirements.txt /code/

RUN pip install -r requirements.txt
ADD . /code/
