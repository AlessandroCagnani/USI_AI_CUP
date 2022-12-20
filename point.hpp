#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <ctime>
#include <string>

using namespace std;

struct point{
    double coord_x;
    double coord_y;
};

int euclidean_distance (const point x, const point y) {
    double ml_cord_x = (x.coord_x - y.coord_x)*(x.coord_x - y.coord_x);
    double ml_cord_y = (x.coord_y - y.coord_y)*(x.coord_y - y.coord_y);
    double cost = sqrt(ml_cord_x + ml_cord_y);
    // cout << "cost: " << cost << endl;
    if (cost == 0)
        return 1;
    else
        return round(cost);
    // return round(cost);
}
