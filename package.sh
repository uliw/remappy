#!/bin/bash
# Note that the pypi config and auth info is in n~.pypirc
cd /home/uliw/user/python-scripts/remappy
rm -rf docs/*
pdoc -o docs ./remappy
rm dist/*
python -m build
echo "Uploading with"
echo "twine upload --repository pypi ./dist/*"
twine upload --repository remappy ./dist/*
