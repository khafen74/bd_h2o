#include "storagemodel.h"

StorageModel::StorageModel(const char *bratPath, const char *outPath, const char *demPath, const char *fdirPath, const char *facPath, double capacity, int type)
{
    init(bratPath, outPath, demPath, fdirPath, facPath, capacity);
    m_nType = type;
    qDebug()<<"type"<<m_nType;
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
    setOutputPaths(pondExtents);
    createHandInputs();
    calcSurfaceWSE();

    Raster_BeaverPond rasterBP;
    for (int i=0; i<m_qvHandIn.length(); i++)
    {
        rasterBP.heightAboveNetwork(m_qvSurfaceWSEPaths[i].toStdString().c_str(), m_fdirPath, m_qvHandIn[i].toStdString().c_str(),
                                       m_qvWSEPaths[i].toStdString().c_str(), m_qvGWPondID[i].toStdString().c_str());
    }

    calcWSEChange();
}

void StorageModel::calcSurfaceWSE()
{
    Raster raster;
    for (int i=0; i<m_qvSurfaceDepthPaths.length(); i++)
    {
        raster.add(m_demPath, m_qvSurfaceDepthPaths[i].toStdString().c_str(), m_qvSurfaceWSEPaths[i].toStdString().c_str());
    }
}

void StorageModel::calcWSEChange()
{
    QString startHAND = m_absPath + "/WSE_start.tif";
    Raster raster;
    Raster_BeaverPond rasterBP;
    raster.heightAboveNetwork(m_demPath, m_fdirPath, m_facPath, startHAND.toStdString().c_str());

    for (int i=0; i<m_qvGWChange.length(); i++)
    {
        rasterBP.subtractHAND(startHAND.toStdString().c_str(), m_qvWSEPaths[i].toStdString().c_str(), m_qvGWChange[i].toStdString().c_str());
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

void StorageModel::createHandInputs()
{
    Raster_BeaverPond rasterBP;

    for (int i=0; i<m_qvPondPaths.length(); i++)
    {
        rasterBP.createHANDInput(m_qvPondPaths[i].toStdString().c_str(), m_facPath, m_qvHandIn[i].toStdString().c_str());
    }
}

void StorageModel::createModflowInputs(DamPolygons pondExtents)
{
    setOutputPaths(pondExtents);
    qDebug()<<"creating modflow inputs";
    calcSurfaceWSE();
    qDebug()<<"surface wse created";
    Raster_BeaverPond raster_bp;
    qDebug()<<"starting head calculations";
    raster_bp.head(m_demPath, m_facPath, m_qvHead[0].toStdString().c_str());
    qDebug()<<"starting loop";
    for (int i=0; i<m_qvSurfaceWSEPaths.length(); i++)
    {
        raster_bp.head(m_qvSurfaceWSEPaths[i].toStdString().c_str(), m_facPath, m_qvPondPaths[i].toStdString().c_str(), m_qvHead[i+1].toStdString().c_str());
    }
}

void StorageModel::run()
{
    cleanOutDir();
    DamPoints pondPoints(m_demPath, m_bratPath, m_facPath, m_outPath, bratCap);
    DamPolygons pondPolys(pondPoints, m_nType, m_fdirPath);
    //calcFinalWSE(pondPolys);
    ReachLines reachStorage(pondPoints);
    //setOutputPaths(pondPolys);
    createModflowInputs(pondPolys);
}

void StorageModel::runFromPoints(const char *damsIn, const char *csvOut, int nRunType)
{
    cleanOutDir();
    if (nRunType == 2)
    {
        m_nType = 2;
    }
    qDebug()<<"starting dam points";
    DamPoints pondPoints(m_demPath, m_bratPath, m_facPath, m_outPath, bratCap, damsIn, nRunType);
    qDebug()<<"starting pond polys";
    DamPolygons pondPolys(pondPoints, m_nType, m_fdirPath);
    ReachLines reachStorage(pondPoints);
    createModflowInputs(pondPolys);
    //pondPoints.compareArea(damsIn, csvOut);
}

void StorageModel::runFromPointsWithHeights(const char *damsIn, const char *csvOut)
{
    cleanOutDir();
    int nDamRunType = 4;
    if (nDamRunType == 4)
    {
        m_nType = 2;
    }
    DamPoints pondPoints(m_demPath, m_bratPath, m_facPath, m_outPath, bratCap, damsIn, nDamRunType);
    DamPolygons pondPolys(pondPoints, m_nType, m_fdirPath);
    ReachLines reachStorage(pondPoints);
    setOutputPaths(pondPolys);
    createModflowInputs(pondPolys);
}

void StorageModel::setOutputPaths(DamPolygons pondExtents)
{
    m_qvPondPaths.clear(), m_qvSurfaceDepthPaths.clear(), m_qvSurfaceWSEPaths.clear(), m_qvWSEPaths.clear()
            , m_qvHandIn.clear(), m_qvGWPondID.clear(), m_qvGWChange.clear(), m_qvHead.clear();
    m_qvPondPaths.append(pondExtents.getLoPondPath()), m_qvPondPaths.append(pondExtents.getMidPondPath()), m_qvPondPaths.append(pondExtents.getHiPondPath());
    m_qvSurfaceDepthPaths.append(pondExtents.getLoDepthPath()), m_qvSurfaceDepthPaths.append(pondExtents.getMidDepthPath()), m_qvSurfaceDepthPaths.append(pondExtents.getHiDepthPath());
    QFileInfo fi(pondExtents.getHiDepthPath());
    QString absPath = fi.absolutePath();
    m_absPath = absPath;
    m_qvSurfaceWSEPaths.append(absPath+"/WSESurf_lo.tif"), m_qvSurfaceWSEPaths.append(absPath+"/WSESurf_mid.tif"), m_qvSurfaceWSEPaths.append(absPath+"/WSESurf_hi.tif");
    m_qvWSEPaths.append(absPath+"/WSE_lo.tif"), m_qvWSEPaths.append(absPath+"/WSE_mid.tif"), m_qvWSEPaths.append(absPath+"/WSE_hi.tif");
    m_qvHandIn.append(absPath+"/HAND_lo.tif"), m_qvHandIn.append(absPath+"/HAND_mid.tif"), m_qvHandIn.append(absPath+"/HAND_hi.tif");
    m_qvGWPondID.append(absPath+"/GWPondID_lo.tif"),  m_qvGWPondID.append(absPath+"/GWPondID_mid.tif"),  m_qvGWPondID.append(absPath+"/GWPondID_hi.tif");
    m_qvGWChange.append(absPath+"/WSEChange_lo.tif"), m_qvGWChange.append(absPath+"/WSEChange_mid.tif"), m_qvGWChange.append(absPath+"/WSEChange_hi.tif");
    m_qvHead.append(absPath+"/head_start.tif"), m_qvHead.append(absPath+"/head_lo.tif"), m_qvHead.append(absPath+"/head_mid.tif")
            , m_qvHead.append(absPath+"/head_hi.tif");
}
