#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>
#include <QtCore>
#include "ogrsf_frmts.h"

/*
 * THINGS TO FIX
 *
 * 1. in calcAzimuth, add conditions for horizontal and vertical azimuths
 *
 */

const double PI = 3.14159265;

class Geometry
{
public:
    Geometry();

    static double addDegrees(double base, double addValue);
    static double angleBetweenLines(double x1, double y1, double x2, double y2, double x3, double y3);
    static double calcAzimuth(double startX, double startY, double endX, double endY);
    static int calcCoords(double startX, double startY, double azimuth, double distance, double &newX, double &newY);
    static bool pointInPolygon(OGRLinearRing *pRing, double x, double y);
};

#endif // GEOMETRY_H
