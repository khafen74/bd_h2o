#include "storagemodel.h"

StorageModel::StorageModel(const char *bratPath, const char *outPath, const char *demPath, double capacity)
{
    init(bratPath, outPath, demPath, capacity);
}

void StorageModel::init(const char *bratPath, const char *outPath, const char *demPath, double capacity)
{
    m_bratPath = bratPath;
    m_outPath = outPath;
    m_demPath = demPath;
    bratCap = capacity;
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
