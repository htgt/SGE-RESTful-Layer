# SGE RESTful Layer

## Install

```
python3 -m venv venv
. venv/bin/activate
pip install -r requirements.txt
```

## Run
```
. venv/bin/activate
flask --app app run --host=0.0.0.0 --port=8080
```

## Run linter
```
pycodestyle --statistics src
```
or to run and show the error and description how to solve the error
```
pycodestyle --show-source --show-pep8 src
```