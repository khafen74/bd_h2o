#ifndef DAMPOINTS_H
#define DAMPOINTS_H

#include "ogrsf_frmts.h"
#include "ogr_core.h"
#include "ogr_api.h"
#include "geometry.h"
#include "statistics.h"
#include "raster.h"

class DamPoints
{
public:
    DamPoints(const char *pointsPath);
    DamPoints(const char *demPath, const char *bratPath, const char *outDirPath, double modCap);

    void init(const char *bratPath);

    void createDamPoints_BRAT(OGRLayer *pBratLyr, OGRLayer *pDamsLyr);
    void createDamPoints_Copy();
    void createFields(OGRLayer *pLayer);
    const char* getDemPath();
    const char* getLayerName();
    const char* getOutDirPath();
    void loadDriver();
    void setBratCapacity(double capacity);
    void setDemPath(const char *demPath);
    void setDamHeights(OGRFeature *pFeat, double low, double mid, double high, double max);
    void setFieldValues(OGRFeature *pFeat, int bratID, double groundElev, double slope, double azimuth, double ptX, double ptY);
    void setOutDir(const char *outDirPath);

    static void setPondAttributes(OGRFeature *pFeat, double lowarea, double midarea, double hiarea, double lowvol, double midvol, double hivol);

private:
    OGRSFDriver *m_pDriverShp;

    double m_modCap, m_meanDamHeight, m_confHi, m_confLo;
    const char *m_outDir, *m_layerName, *m_demPath;
};

#endif // DAMPOINTS_H
