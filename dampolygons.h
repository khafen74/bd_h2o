#ifndef DAMPOLYGONS_H
#define DAMPOLYGONS_H

#include "dampoints.h"

class DamPolygons
{
public:
    DamPolygons(DamPoints pts, double maxDamHt);

    void setMaxDamHeight(double maxDamHt);

private:
    double m_maxDamHt;
};

#endif // DAMPOLYGONS_H
