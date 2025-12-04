#include <cmath>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern "C" double calculateAngularDistance(double RA1, double DEC1, double RA2, double DEC2){
    // Convert degrees to radians
    double radRA1 = RA1 * M_PI / 180.0;
    double radDEC1 = DEC1 * M_PI / 180.0;
    double radRA2 = RA2 * M_PI / 180.0;
    double radDEC2 = DEC2 * M_PI / 180.0;
    // Calculate angular distance using the haversine formula
    double deltaRA = radRA2 - radRA1;
    double deltaDEC = radDEC2 - radDEC1;
    double havTheta = sin(deltaDEC / 2.0) * sin(deltaDEC / 2.0) + cos(radDEC1) * cos(radDEC2) * sin(deltaRA / 2.0) * sin(deltaRA / 2.0);
    double radAngDist = 2.0 * asin(sqrt(havTheta));
    // Convert back to degrees
    double angDist = radAngDist * 180.0 / M_PI;
    return angDist;
}