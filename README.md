# epsscan
This is a python wrapper around several anomaly detection algorithms written in c++. To compile this you will need:
boost.python
python 2.7
scons

#Instructions for OSX
Assuming you have brew installed:

> brew install python
> brew install boost --with-python
> brew install scons

Then you can cd to the directory and run:

> scons

This will create a build directory containing build/release and build/debug builds. You can use this library by adding the
build/release/pywrapper directory to your python path. So:

> export PYTHONPATH=build/release/pywrapper:${PYTHONPATH}

You should then be able to use this library as a standard python module by doing:
import eps_scan

# Calling eps scan 


