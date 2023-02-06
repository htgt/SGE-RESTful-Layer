FROM python:3.8-slim-bullseye as base

EXPOSE 8081
WORKDIR /app

ENV PYTHONUNBUFFERED: 1

COPY . .
RUN make
RUN make install
RUN make test

COPY . .

CMD [ "python3", "-m" , "gunicorn", "--bind=0.0.0.0:8081", "src.app:app"]


# syntax=docker/dockerfile:1
