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
    void loadDriver_GDAL();
    void calculateWaterDepth(OGRLayer *pPts, OGRLayer *pPolys);
    void createDepthRasters();
    void createFields(OGRLayer *pLayer);
    void createPondPolygons(OGRLayer *pPts, OGRLayer *pPolys);
    double getUpstreamDistance(double dist);
    void setDemPath(const char *demPath);
    void setFieldValues(OGRFeature *pFeat, int bratID, int pondID);
    void setMaxDamHeight(double maxDamHt);
    void setMaxDistance(double distance);
    void setOutDirPath(const char *outDir);
    void setRasterPaths();

private:
    OGRSFDriver *m_pDriverShp;
    GDALDriver *m_pDriverTiff;

    double m_maxDist, m_cellWidth, m_cellHeight;
    const char *m_outDir, *m_demPath, *m_layerName;
    QString m_qsMid, m_qsLo, m_qsHi;
};

#endif // DAMPOLYGONS_H
