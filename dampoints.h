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
    DamPoints(const char *demPath, const char *bratPath, const char *outDirPath, double modCap);
    DamPoints(const char *demPath, const char *bratPath, const char *outDirPath, double modCap, const char *exPath);

    void init(const char *bratPath);
    void init(const char *bratPath, const char *exPath);

    void compareArea(const char *damsIn, const char *csvOut);
    void createDamPoints_BRAT(OGRLayer *pBratLyr, OGRLayer *pDamsLyr);
    void createDamPoints_Copy(OGRLayer *pBratLyr, OGRLayer *pDamsLyr, OGRLayer *pExLyr);
    void createFields(OGRLayer *pLayer);
    QString getBratDir();
    QString getBratName();
    const char* getDemPath();
    const char* getLayerName();
    const char* getOutDirPath();
    void loadDriver();
    void setBratCapacity(double capacity);
    void setBratPath(QString path);
    void setBratName(QString name);
    void setDemPath(const char *demPath);
    void setDamHeights(OGRFeature *pFeat, double low, double mid, double high, double max);
    void setFieldValues(OGRFeature *pFeat, int bratID, double groundElev, double slope, double azimuth, double ptX, double ptY);
    void setOutDir(const char *outDirPath);

    static void setPondAttributes(OGRFeature *pFeat, double lowarea, double midarea, double hiarea, double lowvol, double midvol, double hivol);

private:
    OGRSFDriver *m_pDriverShp;

    double m_modCap, m_meanDamHeight, m_confHi, m_confLo;
    const char *m_outDir, *m_layerName, *m_demPath;
    QString m_qsBratDir, m_qsBratName;
};

#endif // DAMPOINTS_H
