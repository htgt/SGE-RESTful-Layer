.ONESHELL:
SHELL := /bin/bash

VENV = venv
PYTHON = $(VENV)/bin/python
PIP = $(VENV)/bin/pip

ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

APP = $(PREFIX)/src/app


init: 
	git config core.hooksPath .githooks
	chmod +x .githooks/*

install: 
	echo "Installing..."
	if [ "$(shell which sudo)" = "" ]; then
        $(MAKE) install-sudo;
    fi
	sudo apt-get update

	if [ "$(shell which python3.8-dev)" = "" ]; then
        $(MAKE) install-python3.8-dev;
    fi
	if [ "$(shell which python3.8-venv)" = "" ]; then
        $(MAKE) install-python3.8-venv;
    fi
	if [ "$(shell which libglib2.0-dev)" = "" ]; then
        $(MAKE) install-libglib2.0-dev;
    fi
	if [ "$(shell which autoconf)" = "" ]; then
        $(MAKE) install-autoconf;
    fi

install-sudo:
	echo "Installing sudo..."
	apt-get -y install sudo

install-python3.8-venv:
	echo "Installing python3.8-venv..."
	sudo apt-get -y install python3.8-venv

install-python3.8-dev: 
	if [ "$(shell which python3)" != "" ] && ["$(shell python3 -v)" >= 3.8 ]; then
		PYTHONPATH = which python
		PYTHONPATH38 = which python3
		sudo update-alternatives --install ${PYTHONPATH} python ${PYTHONPATH38} 2
	else
		echo "Installing python3.8-dev..."
		sudo apt-get -y install python3.8-dev
		PYTHONPATH = which python
		PYTHONPATH38 = which python3.8
	fi
	sudo update-alternatives --install ${PYTHONPATH} python ${PYTHONPATH38} 2 
	sudo update-alternatives --config python 


install-libglib2.0-dev: 
	echo "Installing libglib2.0-dev..."
	sudo apt-get -y install libglib2.0-dev 

install-autoconf: 
	echo "Installing autoconf..."
	sudo apt-get -y install autoconf libtool

setup-venv:
	echo "running venv"
	python -m venv venv
	./venv/bin/pip install -U pip wheel setuptools 
	./venv/bin/pip install -r requirements.txt
	./venv/bin/pip install -r sge-primer-scoring/requirements.txt

test:
	. venv/bin/activate
	&& python -m unittest

run:
	. venv/bin/activate
	&& flask --app src/app run --host=0.0.0.0 --port=8080
	&& gunicorn --bind 0.0.0.0:5000 src.app:app


clean: 
	rm -rf __pycache__
	rm -rf venv