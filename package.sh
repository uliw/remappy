#!/bin/bash
# Note that the pypi config and auth info is in n~.pypirc
cd /home/uliw/user/python-scripts/remappy
rm -rf docs/*
pdoc -o docs esbmtk
rm dist/*
python setup.py clean --all
python3 setup.py sdist bdist_wheel
echo "Uploading with"
echo "twine upload --repository pypi ./dist/*"
twine upload --repository remappy ./dist/*
