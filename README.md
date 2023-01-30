# SGE RESTful Layer

## Install

```
python3 -m venv venv
. venv/bin/activate
pip install -r requirements.txt
```

## Run Application
```
. venv/bin/activate
flask --app src/app run --host=0.0.0.0 --port=8080
```

## Run with Gunicorn locally
```
gunicorn --bind 0.0.0.0:5000 src.app:app
```

## Controlling the Gunicorn service
```
sudo systemctl stop sge
sudo systemctl status sge
sudo systemctl start sge

Service config can be found here:
/etc/systemd/system/sge.service
```

## Run in Docker

Build image

```docker build -t sge-restful-layer . ```


Run container

```docker run -d -p 8081:8081 sge-restful-layer ```

## Run Unit Tests
First make sure you are running venv! 

If not:
```
. venv/bin/activate
```

Then 

```
python -m unittest discover
```

## Run linter
```
pycodestyle --statistics src
```
or to run and show the error and description how to solve the error
```
pycodestyle --show-source --show-pep8 src
```

