#include "storagemodel.h"

StorageModel::StorageModel(const char *bratPath, const char *outPath, const char *demPath, const char *fdirPath, const char *facPath, double capacity)
{
    init(bratPath, outPath, demPath, fdirPath, facPath, capacity);
}

void StorageModel::init(const char *bratPath, const char *outPath, const char *demPath, const char *fdirPath, const char *facPath, double capacity)
{
    m_bratPath = bratPath;
    m_outPath = outPath;
    m_demPath = demPath;
    m_fdirPath = fdirPath;
    m_facPath = facPath;
    bratCap = capacity;
}

void StorageModel::calcFinalWSE(DamPolygons pondExtents)
{
    qDebug()<<"setting output paths";
    setOutputPaths(pondExtents);
    qDebug()<<"starting surface WSE";
    calcSurfaceWSE();
    qDebug()<<"starting final WSE";
}

void StorageModel::calcSurfaceWSE()
{
    Raster raster;
    for (int i=0; i<m_qvSurfaceDepthPaths.length(); i++)
    {
        raster.add(m_demPath, m_qvSurfaceDepthPaths[i].toStdString().c_str(), m_qvSurfaceWSEPaths[i].toStdString().c_str());
    }
}

void StorageModel::cleanOutDir()
{
    QString path = QString::fromUtf8(m_outPath);
    QDir dir(path);
    dir.setNameFilters(QStringList() << "*.*");

    foreach (QString dirFile, dir.entryList())
    {
        dir.remove(dirFile);
    }
}

void StorageModel::run()
{
    cleanOutDir();
    DamPoints pondPoints(m_demPath, m_bratPath, m_outPath, bratCap);
    DamPolygons pondPolys(pondPoints);
    calcFinalWSE(pondPolys);
    ReachLines reachStorage(pondPoints);
}

void StorageModel::runFromPoints(const char *damsIn, const char *csvOut)
{
    cleanOutDir();
    DamPoints pondPoints(m_demPath, m_bratPath, m_outPath, bratCap, damsIn);
    DamPolygons pondPolys(pondPoints);
    ReachLines reachStorage(pondPoints);
    pondPoints.compareArea(damsIn, csvOut);
}

void StorageModel::setOutputPaths(DamPolygons pondExtents)
{
    m_qvPondPaths.clear(), m_qvSurfaceDepthPaths.clear(), m_qvSurfaceWSEPaths.clear(), m_qvWSEPaths.clear();
    m_qvPondPaths.append(pondExtents.getLoPondPath()), m_qvPondPaths.append(pondExtents.getMidPondPath()), m_qvPondPaths.append(pondExtents.getHiPondPath());
    m_qvSurfaceDepthPaths.append(pondExtents.getLoDepthPath()), m_qvSurfaceDepthPaths.append(pondExtents.getMidDepthPath()), m_qvSurfaceDepthPaths.append(pondExtents.getHiDepthPath());
    QFileInfo fi(pondExtents.getHiDepthPath());
    QString absPath = fi.absolutePath();
    m_qvSurfaceWSEPaths.append(absPath+"/WSESurf_lo.tif"), m_qvSurfaceWSEPaths.append(absPath+"/WSESurf_mid.tif"), m_qvSurfaceWSEPaths.append(absPath+"/WSESurf_hi.tif");
    m_qvWSEPaths.append(absPath+"/WSE_lo.tif"), m_qvWSEPaths.append(absPath+"/WSE_mid.tif"), m_qvWSEPaths.append(absPath+"/WSE_hi.tif");
}
