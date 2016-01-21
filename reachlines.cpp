#include "reachlines.h"

ReachLines::ReachLines(DamPoints pondPts)
{
    init(pondPts);
}

void ReachLines::init(DamPoints pondPts)
{
    loadDriver();
    setOutDirPath(pondPts.getOutDirPath());
    m_qsBratDir = pondPts.getBratDir();
    m_qsBratName = pondPts.getBratName();
    m_qsLayerName = "ModeledReachStorage";

    OGRDataSource *pOutDs, *pBratDs;
    OGRLayer *pBrat_lyr, *pReach_lyr, *pPts_lyr;
    pOutDs = m_pDriverShp->CreateDataSource(m_qsOutDir.toStdString().c_str(), NULL);
    pBratDs = m_pDriverShp->CreateDataSource(m_qsBratDir.toStdString().c_str(), NULL);
    pPts_lyr = pOutDs->GetLayerByName(pondPts.getLayerName());
    pBrat_lyr = pBratDs->GetLayerByName(m_qsBratName.toStdString().c_str());
    pReach_lyr = pOutDs->CreateLayer(m_qsLayerName.toStdString().c_str(), pBrat_lyr->GetSpatialRef(), wkbMultiLineString, NULL);
    summarizeReachDepths(pPts_lyr, pReach_lyr, pBrat_lyr);

    OGRDataSource::DestroyDataSource(pOutDs);
    OGRDataSource::DestroyDataSource(pBratDs);
}

void ReachLines::createFields(OGRLayer *pLayer)
{
    OGRFieldDefn field("brat_ID", OFTInteger);
    pLayer->CreateField(&field);
    field.SetName("area_lo");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("area_mid");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("area_hi");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_lo");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_mid");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_hi");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
}

void ReachLines::loadDriver()
{
    OGRRegisterAll();
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    m_pDriverShp = registrar->GetDriverByName("ESRI Shapefile");
}

void ReachLines::setOutDirPath(const char *outDir)
{
    m_qsOutDir = QString::fromUtf8(outDir);
}

void ReachLines::summarizeReachDepths(OGRLayer *pPts, OGRLayer *pReach, OGRLayer *pBrat)
{
    createFields(pReach);
    int nReaches = pBrat->GetFeatureCount();

    QVector<int> reachIDs(nReaches);
    QVector<double> areasLo(nReaches, 0.0), volumesLo(nReaches, 0.0), areasMid(nReaches, 0.0), volumesMid(nReaches, 0.0), areasHi(nReaches, 0.0), volumesHi(nReaches, 0.0);

    for (int i=0; i<nReaches; i++)
    {
        OGRFeature *pBratFeat = pBrat->GetFeature(i);
        reachIDs[i] = pBratFeat->GetFID();
        OGRFeature::DestroyFeature(pBratFeat);
    }

    for (int i=0; i<pPts->GetFeatureCount(); i++)
    {
        OGRFeature *pDamFeat = pPts->GetFeature(i);
        int index = reachIDs.indexOf(pDamFeat->GetFieldAsInteger("brat_ID"));
        areasLo[index] += pDamFeat->GetFieldAsDouble("area_lo");
        areasMid[index] += pDamFeat->GetFieldAsDouble("area_mid");
        areasHi[index] += pDamFeat->GetFieldAsDouble("area_hi");
        volumesLo[index] += pDamFeat->GetFieldAsDouble("vol_lo");
        volumesMid[index] += pDamFeat->GetFieldAsDouble("vol_mid");
        volumesHi[index] += pDamFeat->GetFieldAsDouble("vol_hi");
        OGRFeature::DestroyFeature(pDamFeat);
    }

    for (int i=0; i<nReaches; i++)
    {
        OGRFeature *pBratFeat = pBrat->GetFeature(i);
        OGRFeature *pReachFeat = OGRFeature::CreateFeature(pReach->GetLayerDefn());
        pReachFeat->SetGeometry(pBratFeat->GetGeometryRef());
        int index = reachIDs.indexOf(pBratFeat->GetFID());
        pReachFeat->SetField("brat_ID", pBratFeat->GetFID()*1.0);
        pReachFeat->SetField("area_lo", areasLo[index]);
        pReachFeat->SetField("area_mid", areasMid[index]);
        pReachFeat->SetField("area_hi", areasHi[index]);
        pReachFeat->SetField("vol_lo", volumesLo[index]);
        pReachFeat->SetField("vol_mid", volumesMid[index]);
        pReachFeat->SetField("vol_hi", volumesHi[index]);
        pReach->CreateFeature(pReachFeat);
        OGRFeature::DestroyFeature(pReachFeat);
    }
}

