FROM node:wheezy

RUN mkdir /code
WORKDIR /code
ADD package.json .
RUN yarn install
ADD . .
