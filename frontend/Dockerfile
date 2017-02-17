FROM node:wheezy

RUN npm install http-server -g
RUN mkdir /code
ADD . /code
WORKDIR /code
RUN npm install
