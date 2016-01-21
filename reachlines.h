#ifndef REACHLINES_H
#define REACHLINES_H

#include "dampoints.h"

class ReachLines
{
public:
    ReachLines(DamPoints pondPts);

    void init(DamPoints pondPts);

    void createFields(OGRLayer *pLayer);
    void loadDriver();
    void setOutDirPath(const char *outDir);
    void summarizeReachDepths(OGRLayer *pPts, OGRLayer *pReach, OGRLayer *pBrat);

private:
    OGRSFDriver *m_pDriverShp;
    QString m_qsOutDir, m_qsBratDir, m_qsLayerName, m_qsBratName;
};

#endif // REACHLINES_H
