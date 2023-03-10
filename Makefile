.ONESHELL:
SHELL := /bin/bash

VENV = venv
PYTHON = $(VENV)/bin/python
PIP = $(VENV)/bin/pip

# Docker
name = "sge-restful-layer"
tag = "test"

ifeq ($(PREFIX),)
	PREFIX := /usr/local
endif

APP = $(PREFIX)/src/app
BENCHLING_CONFIG_FILE = $(PREFIX)/src/benchling/config.cfg


init: 
	@git config core.hooksPath .githooks
	@chmod +x .githooks/*

install: 
	@echo "Installing..."
	@if [ "$(shell which sudo)" = "" ]; then
		$(MAKE) install-sudo;
	fi
	@sudo apt-get update

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

install-@sudo:
	@echo "Installing @sudo..."
	@apt-get -y install @sudo

install-python3.8-venv:
	@echo "Installing python3.8-venv..."
	@sudo apt-get -y install python3.8-venv

install-python3.8-dev: 
	@if [ "$(shell which python3)" = "" ]; then
		# Python3 not installed.
		@echo "Installing python3.8-dev..."
		@sudo apt-get -y install python3.8-dev
	else
		PYTHONPATH = which python
		ver=$(python3 -V 2>&1 | sed 's/.* \([0-9]\).\([0-9]\).*/\1\2/')
		@if [ "${ver}" -ge "38" ]; then
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

$(BENCHLING_CONFIG_FILE): 
	if  [[ ! -f "${BENCHLING_CONFIG_FILE}" ]]; then
		echo "THISISASECRETKEYAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" > $(BENCHLING_CONFIG_FILE)
	fi

test: setup-venv
	@. venv/bin/activate
	@python -m unittest

run-flask: ./src/benchling/config.cfg setup-venv
	@. venv/bin/activate
	@flask --app src/app run --host=0.0.0.0 --port=8081

run-gunicorn: ./src/benchling/config.cfg setup-venv
	@. venv/bin/activate 
	@python -m gunicorn --bind 0.0.0.0:8081 src.app:app

docker_touch:
	@docker build --pull -t "${name}:${tag}" .
	touch docker_touch
build_docker: docker_touch

run-docker: build_docker
	@docker run -p 8081:8081 -t "${name}:${tag}"

clean: 
	@rm -rf __pycache__
	@rm -rf venv