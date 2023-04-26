FROM python:3.8.0 as base

EXPOSE 8081
WORKDIR /app

COPY requirements.txt requirements.txt
COPY Makefile Makefile
COPY src src
COPY tests tests
COPY schemas schemas

RUN make install
RUN make setup-venv

FROM base as gunicorn
ARG GUNICORN_CONF_FILE=gunicorn.conf.py
COPY gunicorn.conf.py gunicorn.conf.py

FROM base as unittest
ENV BENCHLING_TENANT=${BENCHLING_TENANT:-unittest}
ENV DOCKER_ENV=${DOCKER_ENV:-unittest}
CMD [ "sh", "-c", "make test"]

FROM base as local
COPY .env .env
CMD [ "sh", "-c", "make run-gunicorn", "-c", "gunicorn.conf.py"]

FROM base as remote
CMD [ "sh", "-c", "make run-gunicorn", "-c", "gunicorn.conf.py"]