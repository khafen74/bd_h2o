#include "geometry.h"

Geometry::Geometry()
{

}

double Geometry::addDegrees(double base, double addValue)
{
    double value, remainder;

    value = base + addValue;

    if (value > 360.0)
    {
        remainder = value - 360.0;
        value = remainder;
    }
    else if (value < 0.0)
    {
        remainder = fabs(value);
        value = 360.0 - remainder;
    }

    return value;
}

double Geometry::angleBetweenLines(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double angle;

    angle = atan2(y1-y3, x1-x3) - atan2(y2-y3, x2-x3);

    while (angle < (-1.0*PI))
    {
        angle += 2.0*PI;
    }
    while (angle > PI)
    {
        angle -= 2.0*PI;
    }

    //qDebug()<<angle;
    return angle;
}

double Geometry::calcAzimuth(double startX, double startY, double endX, double endY)
{
    double theta, azimuth;

    if (startX > endX && startY > endY)
    {
        theta = atan2((startY-endY),(startX-endX)) * 180.0 / PI;
        azimuth = 180.0 + theta;
    }
    else if (startX < endX && startY > endY)
    {
        theta = atan2((startY-endY),(endX-startX)) * 180.0 / PI;
        azimuth = 360.0 - theta;
    }
    else if (startX < endX && startY < endY)
    {
        theta = atan2((endY-startY),(endX-startX)) * 180.0 / PI;
        azimuth = theta;
    }
    else if (startX > endX && startY < endY)
    {
        theta = atan2((startX-endX),(endY-startY)) * 180.0 / PI;
        azimuth = 90.0 + theta;
    }
    else
    {
        qDebug()<<"AZIMUTH ERROR: may be a straight line";
    }

    return azimuth;
}

int Geometry::calcCoords(double startX, double startY, double azimuth, double distance, double &newX, double &newY)
{
    double deltaX, deltaY, theta;

    if (azimuth > 0.0 && azimuth < 90.0)
    {
        theta = azimuth;
        deltaY = sin(theta*(PI/180.0))*distance;
        deltaX = cos(theta*(PI/180.0))*distance;
    }
    else if (azimuth > 90.0 && azimuth < 180.0)
    {
        theta = azimuth - 90.0;
        deltaX = sin(theta*(PI/180.0))*distance*(-1.0);
        deltaY = cos(theta*(PI/180.0))*distance;
    }
    else if (azimuth > 180.0 && azimuth < 270.0)
    {
        theta = azimuth - 180.0;
        deltaY = sin(theta*(PI/180.0))*distance*(-1.0);
        deltaX = cos(theta*(PI/180.0))*distance*(-1.0);
    }
    else if (azimuth > 270.0 && azimuth < 360.0)
    {
        theta = 360.0 - azimuth;
        deltaY = sin(theta*(PI/180.0))*distance*(-1.0);
        deltaX = cos(theta*(PI/180.0))*distance;
    }
    else
    {

    }

    newX = startX + deltaX;
    newY = startY + deltaY;

    return 0;
}

double Geometry::distance_point(double x1, double y1, double x2, double y2)
{
    double xdiff = fabs(x1-x2);
    double ydiff = fabs(y1-y2);
    double dist = sqrt(pow(xdiff, 2.0) + pow(ydiff, 2.0));
    return dist;
}

bool Geometry::pointInPolygon(OGRLinearRing *pRing, double x, double y)
{
    double angle = 0.0;

    for (int i=0; i<pRing->getNumPoints()-1; i++)
    {
        angle += angleBetweenLines(pRing->getX(i), pRing->getY(i), pRing->getX(i+1), pRing->getY(i+1), x, y);
    }
    //qDebug()<<angle;

    if (angle > PI)
    {
        return true;
    }
    else
    {
        return false;
    }
}

