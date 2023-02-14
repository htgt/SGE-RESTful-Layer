FROM python:3.8-slim-bullseye

EXPOSE 8081

ENV PYTHONUNBUFFERED: 1

WORKDIR /usr/src/app
COPY . .
RUN apt-get update && apt-get install build-essential -y
RUN apt-get install -y git sudo

RUN make
RUN sudo make install

COPY . .

CMD [ "make", "setup-venv"]


# syntax=docker/dockerfile:1
