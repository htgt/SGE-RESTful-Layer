.EXPORT_ALL_VARIABLES:
.ONESHELL:
SHELL := /bin/bash

VENV = venv
PYTHON = $(VENV)/bin/python
PIP = $(VENV)/bin/pip

APP = $(PREFIX)/src/app
ENVIRONMENTAL_VARIABLE_FILE := .env

-include ${ENVIRONMENTAL_VARIABLE_FILE}

export GUNICORN_ENV ?= prod
export DOCKER_ENV ?= ${GUNICORN_ENV}
$(info $$DOCKER_ENV = ${DOCKER_ENV})
MAKE_VERSION := $(shell make --version | grep '^GNU Make' | sed 's/^.* //g')
$(info "make version = ${MAKE_VERSION}, minimum version 3.82 required for multiline.")

# Docker
name = "sge-restful-layer"
tag=$$DOCKER_ENV

init: 
	@git config core.hooksPath .githooks
	@chmod +x .githooks/*

install: install-basics
	@echo "Installing..."
	@sudo apt-get update

	@if [ "$(shell which curl)" = "" ]; then
		$(MAKE) install-curl;
	fi
	@if [ "$(shell which python3.8-dev)" = "" ]; then
		$(MAKE) install-python3.8-dev;
	fi
	@if [ "$(shell which python3.8-venv)" = "" ]; then
		$(MAKE) install-python3.8-venv;
	fi
	@if [ "$(shell which libglib2.0-dev)" = "" ]; then
		$(MAKE) install-libglib2.0-dev;
	fi
	@if [ "$(shell which autoconf)" = "" ]; then
		$(MAKE) install-autoconf;
	fi
install-sudo:
	@echo "Installing sudo..."
	@apt-get update
	@apt-get -y install sudo

install-curl:
	@echo "Installing curl..."
	@apt-get -y install curl

install-python3.8-venv:
	@echo "Installing python3.8-venv..."
	@sudo apt-get -y install python3.8-venv

install-python3.8-dev: 
	@if [ "$(shell which python3)" = "" ]; then
		# Python3 not installed.
		@echo "Installing python3.8-dev..."
		@sudo apt-get -y install python3.8-dev
	else
		@PYTHONPATH = which python
		@ver=$$(python3 -V 2>&1 | sed 's/.* \([0-9]\).\([0-9]\).*/\1\2/')
		@if [ "$$ver" -ge 38 ]; then
			PYTHONPATH38 = which python3
		else
			@echo "Installing python3.8-dev..."
			@sudo apt-get -y install python3.8-dev
			PYTHONPATH38 = which python3.8
		fi
		@sudo update-alternatives --install ${PYTHONPATH} python ${PYTHONPATH38} 2 
		@sudo update-alternatives --config python 
	fi

install-libglib2.0-dev: 
	@echo "Installing libglib2.0-dev..."
	@sudo apt-get -y install libglib2.0-dev 

install-autoconf: 
	@echo "Installing autoconf..."
	@sudo apt-get -y install autoconf libtool

install-docker:
	@echo "Installing docker..."
	@curl -fsSL https://get.docker.com -o get-docker.sh
	@sh get-docker.sh
	@sudo groupadd docker
	@sudo usermod -aG docker $$USER
	@newgrp docker


install-basics: install-sudo install-curl install-autoconf
	@sudo apt-get -y install build-essential

venv/bin/activate:
	@python -m venv venv

setup-venv: venv/requirements_run

venv/requirements_run: venv/bin/activate requirements.txt
	@./venv/bin/pip install -U pip wheel setuptools 
	@./venv/bin/pip install -r requirements.txt
	@touch venv/requirements_run

clean-venv/requirements_run:
	@rm -f venv/requirements_run

activate-venv: setup-venv
	@. venv/bin/activate

test: setup-venv
	@. venv/bin/activate
	@python -m unittest

run-flask: setup-venv
	@. venv/bin/activate
	@flask --app src/app run --host=0.0.0.0 --port=8081

run-flask-debug: setup-venv
	@. venv/bin/activate
	@flask --app src/app --debug run --host=0.0.0.0 --port=8081

run-gunicorn: setup-venv
	@. venv/bin/activate 
	@python -m gunicorn src.app:app

docker-touch:
	@ver=$$(docker version --format '{{.Server.Version}}' 2>&1 | sed -E 's/([0-9]+).*/\1/')
	if [ "$$ver" -lt 23 ]; then
		echo Warning Docker engine version $$ver \< 23, changing build to buildx.
		docker buildx install
		DOCKER_BUILDKIT=1
	fi
	@docker build --pull -t "${name}:${tag}" --target base .;
	@touch docker-touch

build-docker: docker-touch

build-docker-gunicorn: build-docker
	@echo Gunicorn tenant = $$GUNICORN_ENV
	@if [[ $$GUNICORN_ENV == prod ]]; then
		@docker build --pull -t "${name}:${tag}" --target gunicorn --build-arg GUNICORN_CONF_FILE=gunicorn.prod.conf.py .;
	else
		@docker build --pull -t "${name}:${tag}" --target gunicorn --build-arg GUNICORN_CONF_FILE=gunicorn.conf.py .;
	fi

build-docker-local: build-docker-gunicorn
	@echo Benchling tenant = $$BENCHLING_TENANT
	@if [[ $$BENCHLING_TENANT == unittest ]]; then
		$(MAKE) build-docker-test
	else
		@docker build --pull -t "${name}:${tag}" --target local .;
	fi

build-docker-remote: build-docker-gunicorn
	@echo Benchling tenant = $$BENCHLING_TENANT
	@if [[ $$BENCHLING_TENANT == unittest ]]; then
		$(MAKE) build-docker-test
	else
		@docker build --pull -t "${name}:${tag}" --target remote .;
	fi

build-docker-test: build-docker
	@docker build --pull -t "${name}:${tag}" --target unittest .;

run-docker: build-docker
	@docker run --name "${name}" -p 8081:8081 -t "${name}:${tag}"

run-docker-local: build-docker-local run-docker

run-docker-remote: build-docker-remote run-docker

run-docker-test: build-docker-test run-docker

run-docker-interactive: build-docker
	@docker run -i --name "${name}" -t "${name}:${tag}" bash

connect-docker-interactive: run-docker
	@docker exec -it ${name} bash

clean-docker:
	@docker builder prune -af
	@docker container prune -f
	@docker image prune -af
	@rm -f docker-touch


check-lint: activate-venv
	@echo "Running pycodestyle for src/"
	@pycodestyle --statistics -qq src || true
	@echo "Running pycodestyle for tests/"
	@pycodestyle --statistics -qq tests || true


auto-lint-tests: activate-venv
	@python -m autopep8 -r -i tests/

auto-lint-src: activate-venv
	@python -m autopep8 -r -i src/

clean:
clean:
	@rm -rf __pycache__
	@rm -rf venv

