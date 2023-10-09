FROM python:3.8.0 as base

EXPOSE 8081
WORKDIR /app

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .

CMD [ "python3", "-m" , "gunicorn", "src.app:app"]


FROM base as unittest
ENV BENCHLING_TENANT=${BENCHLING_TENANT:-unittest}
ENV DOCKER_ENV=${DOCKER_ENV:-unittest}
CMD [ "sh", "-c", "make test"]

