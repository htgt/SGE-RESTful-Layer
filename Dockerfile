FROM python:3.8-slim-bullseye

EXPOSE 8081

ENV PYTHONUNBUFFERED: 1

WORKDIR /
COPY . .
RUN apt update && apt install build-essential -y --no-install-recommends


RUN make
RUN make install
RUN make test

COPY . .

CMD [ "python3", "-m" , "gunicorn", "--bind=0.0.0.0:8081", "src.app:app"]


# syntax=docker/dockerfile:1
