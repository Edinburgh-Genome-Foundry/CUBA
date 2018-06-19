FROM node:9.11.2-jessie

RUN mkdir /code
WORKDIR /code
ADD package.json .
RUN npm install
ADD . .
