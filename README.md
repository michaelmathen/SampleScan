# SampleScan
This is a python wrapper around several anomaly detection algorithms written in c++. To compile this you will need:
boost.python
python 2.7
scons

#Instructions for OSX 
Assuming you have brew installed:
```
> brew install python
> brew install boost --with-python
> brew install scons
```
Then you can cd to the directory and run:
```
> scons
```
This will create a build directory containing build/release and build/debug builds. You can use this library by adding the
build/release/pywrapper directory to your python path. So:
```
> export PYTHONPATH=build/release/pywrapper:${PYTHONPATH}
```
You should then be able to use this library as a standard python module by doing:
import eps_scan
Linux instructions should be basically identical. (I have run this on arch linux and ubuntu)
# Calling eps scan 
Example code:
```
import eps_scan
import random
#Generate a list of random points with random anomalies
pts = [eps_scan.Point(random.random(), random.random(), bool(random.randint(0, 1))) for i in xrange(1000)]
#Run the algorith with a net size of 100 and a sample size of 1000 and print the found region.
#The last parameter prevents us from considering regions that are too large or small since the approximation breaks down there.
print eps_scan.netDisks(pts, 100, 1000, .01) # runs on disks
print eps_scan.netRectsG(pts, 100, 1000, .01) # runs on axis aligned rectangles. 
# For the rectangle algorithm run:
#export OMP_NUM_THREADS=n
# where n is the number of threads you want to use with the algorithm.
```
Try checking the example directory or the testing directory for more examples. The api is defined in the pywrapper directory 
so you can also check that out. 

