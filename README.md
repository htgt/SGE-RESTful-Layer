# SGE RESTful Layer

## Install
With Makefile 
```sh
make
make install
make setup-venv
```
```make``` sets up the git hooks that run unittests and pycodestyle on /src and /tests on ```git push```.
```make install``` installs dependancies below.
```make setup-venv``` creates a venv at ./venv and installs requirements.txt(s)


OR 

```
python3 -m venv venv
. venv/bin/activate
pip install -r requirements.txt
```

## Environmental variables
Ensure that the service you are using has the following enviromental variables or that you have a .env file in the top level of the restful layer directory.
```
BENCHLING_SECRET_KEY="secret_key" or "unittest"
GUNICORN_ENV="prod" or "test" or "unittest"
BENCHLING_TENANT="prod" or "test" or "unittest"
DOCKER_ENV="prod" or "test" or "unittest"
```

## Run Application
There are many ways to run the RESTful layer, it can be run locally with Gunicorn or Flask or in a Docker with Gunicorn.
The easiest way is using the Makefile
```sh
make run
```
```make run-gunicorn``` Launches the app with Gunicorn.
```make run-flask``` Launches the app with flask.
```make run-flask-debug``` Launches the app with flask with debug enabled.
```make run-docker-local``` Creates a Docker container and launches the app with Gunicorn and copies the local .env file.
```make run-docker-remote``` Creates a Docker container and launches the app with Gunicorn.

OR

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

The easiest way to run in a docker container is to use the makefile.
```make run-docker-local```
This builds the container if not created, copies the .env and then runs with gunicorn.

Alternatively you can run without a .env or just run the unittests.
```make run-docker-remote```
```make run-docker-test```

Build image

```docker build -t sge-restful-layer . ```

Or with makefile
```make build-docker-local```
or 
```make build-docker-remote```


Run container

```docker run -p 8081:8081 sge-restful-layer ```

Or with makefile (also builds)
```make run-docker-local```
```make run-docker-remote```

clean docker with
```make clean-docker```

## Run Unit Tests

With Makefile 
```make test``` Runs the various unittests.
```make run-docker-test``` Runs in a Ubuntu python 3.8 based docker container.


OR

First make sure you are running venv! 

If not:
```
. venv/bin/activate
```

Then 

```
python -m unittest discover -v
```

## Run linter
```
pycodestyle --statistics src
```
or to run and show the error and description how to solve the error
```
pycodestyle --show-source --show-pep8 src
```

Or with the makefile

```
make check-lint
```

This follows the pep8 convention with modifications documented in setup.cfg

## Run Autolinter

Run autopep8 inplace with recursion with the makefile

```
make auto-lint-tests
make auto-lint-src
```

This follows the pep8 convention with modifications documented in setup.cfg