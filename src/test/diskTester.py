import eps_scan
import random
from itertools import combinations_with_replacement, combinations
import matplotlib.pyplot as pyplot
import matplotlib.patches as patches

Npts = [eps_scan.Point(random.random(), random.random(), bool(random.randint(0, 1))) for i in xrange(10)]
sA = [eps_scan.Point(random.random(), random.random(), True) for i in xrange(100)]
sB = [eps_scan.Point(random.random(), random.random(), False) for i in xrange(100)]


def findMaxDisk(Npts, sA, sB):
    """
    Computes the max disk in a brain dead slow way to 
    ensure correctness.
    """
    maxReg = None
    maxVal = 0
    for (i, j, k) in combinations(Npts, 3):
        newReg = eps_scan.makeDisk3(i, j, k)
        # Slow step
        sACount = sum(int(newReg.inside(pt)) for pt in sA)
        sBCount = sum(int(newReg.inside(pt)) for pt in sB)
        newReg.setNumAnomalies(sACount, len(sA))
        newReg.setNumPoints(sBCount, len(sB))
        if newReg.statistic() >= maxVal:
            maxReg = newReg
            maxVal = newReg.statistic()
    return maxReg
        
reg1 = eps_scan.netDisksSample(Npts, sA, sB, .00001)
reg2 = findMaxDisk(Npts, sA, sB)

sACount = sum(int(reg1.inside(pt)) for pt in sA)
sBCount = sum(int(reg1.inside(pt)) for pt in sB)
print reg1.numAnomalies(), sACount
print reg1.numPoints(), sBCount
print reg1.totalAnomalies(), len(sA)
print reg1.totalPoints(), len(sB)

print "alldisks", reg1.statistic()
print "maxDisks", reg2.statistic()

print(reg1)
print(reg2)

fig = pyplot.figure()
ax = fig.add_subplot(111)

xa = [pt.x for pt in sA]
xb = [pt.x for pt in sB]
ya = [pt.y for pt in sA]
yb = [pt.y for pt in sB]
ax.scatter(xa, ya, marker='.', color='red')
ax.scatter(xb, yb, marker='.', color='blue')
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])        
ax.add_patch(
    patches.Circle(
        (reg1.a, reg1.b),   # (x,y)
        radius = reg1.radius,
        fill=False,
        color="green",
        linewidth=3)          # height
)
ax.add_patch(
    patches.Circle(
        (reg2.a, reg2.b),   # (x,y)
        radius = reg2.radius,
        fill=False,
        color="black",
        linewidth=3)          # height
)
ax.axis("tight")
pyplot.show()
