FROM python:3.8.0

EXPOSE 8081
WORKDIR /app

COPY requirements.txt requirements.txt
COPY Makefile Makefile
COPY src src
COPY tests tests

RUN make install
RUN make setup-venv


CMD [ "sh", "-c", "make run-gunicorn" ]