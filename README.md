# riborepa_web
## Install
```
git clone git@github.com:autosome-org/riborepa_web.git riborepa
cd riborepa
python3 -m venv .venv
source .venv/bin/activate
pip install wheel
pip install -r requirements.txt
deactivate
```

## Run
```
source .venv/bin/activate
export FLASK_APP=routes.py
export FLASK_ENV=development
flask run --host 0.0.0.0 --port 9042
```
