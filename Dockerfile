FROM python:3.8.0

EXPOSE 8081
WORKDIR /app

COPY requirements.txt requirements.txt
COPY Makefile Makefile
COPY src src
COPY tests tests

RUN make install
RUN make setup-venv

COPY benchling_schema_ids.json benchling_schema_ids.json
COPY gunicorn.prod.conf.py gunicorn.prod.conf.py
COPY gunicorn.conf.py gunicorn.conf.py

CMD [ "sh", "-c", "make run-gunicorn", "-c", "gunicorn.conf.py"]