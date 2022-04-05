#!/bin/bash
#
cd /home/uliw/user/python-scripts/esbmtk
rm -rf docs/*
pdoc -o docs esbmtk
rm dist/*
python setup.py clean --all
python3 setup.py sdist bdist_wheel
echo "Uploading with"
echo "twine upload --repository pypi ./dist/*"
twine upload --repository pypi ./dist/*
