#include <random>
#include <cmath>
#include <fstream>
#include <string>

#include "algorithms/Point.hpp"
#include "algorithms/Utility.hpp"
#include "algorithms/DiskScan.hpp"
#include "algorithms/Region.hpp"
using namespace anomaly;
using namespace std;
int main() {
  AllDisks allDisks;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> udist(-1.0, 1.0);
  vector<Point> testNet;
  vector<Point> testSampleA;
  vector<Point> testSampleB;

  for (int i = 0; i < 100; i++) {
    Point pt(udist(generator), udist(generator));
    testNet.push_back(pt);
  }

  for (int i = 0; i < 4000; i++) {
    Point pt(udist(generator), udist(generator));
    testSampleA.push_back(pt);
  }

  for (int i = 0; i < 4000; i++) {
    Point pt(udist(generator), udist(generator));
    testSampleB.push_back(pt);
  }
  auto reg = allDisks.run(testNet.begin(), testNet.end(),
			  testSampleA.begin(), testSampleA.end(),
			  testSampleB.begin(), testSampleB.end(), .001);
  cout << reg << endl;
}
