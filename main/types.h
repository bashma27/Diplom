#pragma once
#include <algorithm>
#include <iostream>
#include <fstream>
#include <set>
#include <iomanip>
#include <vector>
#include <unordered_set>
#include <math.h>

using namespace std;

struct Coord2 {
    double x, y;
    Coord2(double _x, double _y) : x(_x), y(_y) {}
    Coord2() : x(0.), y(0.) {}
};

struct Coord3 {
    double x, y, z;
    Coord3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    Coord3() : x(0.), y(0.), z(0.) {}
};