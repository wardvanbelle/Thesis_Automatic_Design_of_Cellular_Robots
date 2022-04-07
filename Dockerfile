FROM julia:latest

RUN mkdir -p /home/app

COPY ./app /home/app
