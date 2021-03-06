#ifndef DAMPOLYGONS_H
#define DAMPOLYGONS_H

#include "dampoints.h"

const double ANGLE_OFFSET[5] = {-90.0, -45.0, 0.0, 45.0, 90.0};

class DamPolygons
{
public:
    DamPolygons(DamPoints pts, int type, const char *fdirPath);

    void init(DamPoints pondPts);

    void loadDriver();
    void loadDriver_GDAL();
    void calculateWaterDepth(OGRLayer *pPts, OGRLayer *pPolys);
    void calculateWaterDepth(OGRLayer *pPts, const char *pondIdPath, const char *htAbovePath);
    void createDepthRasters();
    void createHandInput_ponds(OGRLayer *pPts);
    void createFields(OGRLayer *pLayer);
    void createFields_BRAT(OGRLayer *pLayer);
    void createPondPolygons(OGRLayer *pPts, OGRLayer *pPolys);
    QString getHiDepthPath();
    QString getLoDepthPath();
    QString getMidDepthPath();
    QString getHiPondPath();
    QString getLoPondPath();
    QString getMidPondPath();
    double getUpstreamDistance(double dist);
    void setDemPath(const char *demPath);
    void setFieldValues(OGRFeature *pFeat, int bratID, int pondID);
    void setMaxDamHeight(double maxDamHt);
    void setMaxDistance(double distance);
    void setOutDirPath(const char *outDir);
    void setRasterPaths();
    void summarizePondDepths(OGRLayer *pPts);
    void summarizePondDepths_raster(OGRLayer *pPts);
    void summarizeReachDepths(OGRLayer *pPts);

private:
    OGRSFDriver *m_pDriverShp;
    GDALDriver *m_pDriverTiff;

    int m_nType, m_nSummaryIter;
    double m_maxDist, m_cellWidth, m_cellHeight;
    const char *m_outDir, *m_demPath, *m_layerName, *m_fdirPath;
    QString m_qsMid, m_qsLo, m_qsHi,
            m_qsMidHt, m_qsLoHt, m_qsHiHt,
            m_qsMidHtOut, m_qsLoHtOut, m_qsHiHtOut,
            m_qsHandOut,
            m_qsPondID,
            m_qsMidPond, m_qsLoPond, m_qsHiPond,
            m_qsBratDir, m_qsBratName,
            m_qsDamID,
            m_qsHtAbove;
    QList<QVector<double> > m_qlDams;
};

#endif // DAMPOLYGONS_H
