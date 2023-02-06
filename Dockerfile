FROM python:3.8-slim-bullseye

EXPOSE 8081
WORKDIR /app

ENV PYTHONUNBUFFERED: 1

RUN apt install build-essential -y --no-install-recommends
COPY . .
RUN make
RUN make install
RUN make test

COPY . .

CMD [ "python3", "-m" , "gunicorn", "--bind=0.0.0.0:8081", "src.app:app"]


# syntax=docker/dockerfile:1
