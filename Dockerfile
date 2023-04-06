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

COPY gunicorn.prod.conf.py gunicorn.prod.conf.py
COPY gunicorn.conf.py gunicorn.conf.py

FROM base as prod
CMD [ "sh", "-c", "make run-gunicorn", "-c", "gunicorn.prod.conf.py"]

FROM base as test
CMD [ "sh", "-c", "make run-gunicorn", "-c", "gunicorn.conf.py"]

FROM base as unittest
ENV BENCHLING_TENANT=unittest
ENV DOCKER_ENV=unittest
CMD [ "sh", "-c", "make test"]


# CMD [ "sh", "-c", "make run-gunicorn", "-c", "gunicorn.conf.py"]