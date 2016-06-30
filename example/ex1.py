import eps_scan
import random
import time
#Generate a list of random points with random anomalies
pts = [eps_scan.Point(random.random(), random.random(), bool(random.randint(0, 1))) for i in xrange(1000)]
#Run the algorith with a net size of 100 and a sample size of 1000 and print the found region.
#The last parameter prevents us from considering regions that are too large or small since the approximation breaks down there.
st =  time.time()
print eps_scan.netDisks(pts, 100, 4000, .01)
en = time.time()
print(en - st)


st =  time.time()
print eps_scan.netDisks2(pts, 100, 4000, .01)
en = time.time()
print(en - st)
