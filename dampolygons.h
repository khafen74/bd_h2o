#ifndef DAMPOLYGONS_H
#define DAMPOLYGONS_H

#include "dampoints.h"

const double ANGLE_OFFSET[5] = {-90.0, -45.0, 0.0, 45.0, 90.0};

class DamPolygons
{
public:
    DamPolygons(DamPoints pts);

    void init(DamPoints pondPts);

    void loadDriver();
    void createFields(OGRLayer *pLayer);
    void createPondPolygons(OGRLayer *pPts, OGRLayer *pPolys);
    double getUpstreamDistance(double dist);
    void setDemPath(const char *demPath);
    void setFieldValues(OGRFeature *pFeat, int bratID, int pondID, double maxWSE);
    void setMaxDamHeight(double maxDamHt);
    void setMaxDistance(double distance);
    void setOutDirPath(const char *outDir);

private:
    OGRSFDriver *m_pDriverShp;

    double m_maxDist;
    const char *m_outDir, *m_demPath, *m_layerName;
};

#endif // DAMPOLYGONS_H
