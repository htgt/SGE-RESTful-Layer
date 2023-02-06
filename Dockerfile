FROM python:3.8-slim-bullseye

EXPOSE 8081

ENV PYTHONUNBUFFERED: 1

WORKDIR /usr/src/app
COPY . .
RUN apt-get update && apt-get install build-essential -y
RUN apt-get install -y git sudo

RUN make
RUN sudo make install
RUN make test

COPY . .

CMD [ "python3", "-m" , "gunicorn", "--bind=0.0.0.0:8081", "src.app:app"]


# syntax=docker/dockerfile:1
