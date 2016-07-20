#ifndef STORAGEMODEL_H
#define STORAGEMODEL_H

#include "dampolygons.h"
#include "reachlines.h"

class StorageModel
{
public:
    StorageModel(const char *bratPath, const char *outPath, const char *demPath, const char *fdirPath, const char *facPath, double capacity, int type);

    void init(const char *bratPath, const char *outPath, const char *demPath, const char *fdirPath, const char *facPath, double capacity);

    void calcFinalWSE(DamPolygons pondExtents);
    void calcSurfaceWSE();
    void calcWSEChange();
    void cleanOutDir();
    void createHandInputs();
    void createModflowInputs(DamPolygons pondExtents);
    void run();
    void runFromPoints(const char *damsIn, const char *csvOut, int nRunType=1);
    void runFromPointsWithHeights(const char *damsIn, const char *csvOut);
    void setOutputPaths(DamPolygons pondExtents);

private:
    const char *m_bratPath, *m_outPath, *m_demPath, *m_fdirPath, *m_facPath;
    QString m_absPath;
    QVector<QString> m_qvPondPaths, m_qvSurfaceDepthPaths, m_qvSurfaceWSEPaths, m_qvWSEPaths, m_qvHandIn, m_qvGWPondID, m_qvGWChange, m_qvHead;
    double bratCap;
    int m_nType;
};

#endif // STORAGEMODEL_H
