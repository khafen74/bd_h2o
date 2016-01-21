#include "dampolygons.h"

DamPolygons::DamPolygons(DamPoints pts)
{
    init(pts);
}

void DamPolygons::init(DamPoints pondPts)
{
    loadDriver();
    setOutDirPath(pondPts.getOutDirPath());
    setDemPath(pondPts.getDemPath());
    setMaxDistance(150.0);
    m_layerName = "PondPolygons";
    m_qsBratDir = pondPts.getBratDir();
    m_qsBratName = pondPts.getBratName();

    OGRDataSource *pOutDs;
    OGRLayer *pPoly_lyr, *pPts_lyr;

    pOutDs = m_pDriverShp->CreateDataSource(m_outDir, NULL);
    pPts_lyr = pOutDs->GetLayerByName(pondPts.getLayerName());
    pPoly_lyr = pOutDs->CreateLayer(m_layerName, pPts_lyr->GetSpatialRef(), wkbPolygon, NULL);

    createPondPolygons(pPts_lyr, pPoly_lyr);
    calculateWaterDepth(pPts_lyr, pPoly_lyr);
    summarizePondDepths(pPts_lyr);

    OGRDataSource::DestroyDataSource(pOutDs);
}

void DamPolygons::loadDriver()
{
    OGRRegisterAll();
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    m_pDriverShp = registrar->GetDriverByName("ESRI Shapefile");
}

void DamPolygons::loadDriver_GDAL()
{
    GDALAllRegister();
    m_pDriverTiff = GetGDALDriverManager()->GetDriverByName("GTiff");
}

void DamPolygons::calculateWaterDepth(OGRLayer *pPts, OGRLayer *pPolys)
{
    createDepthRasters();
    OGRPolygon *pPolygon;
    OGRLinearRing *pRing;
    OGREnvelope bound;
    Raster rasDem, rasLo, rasMid, rasHi, loPond, midPond, hiPond, loReach, midReach, hiReach;
    rasDem.setProperties(m_demPath);
    rasLo.setProperties(m_qsLo.toStdString().c_str());
    rasMid.setProperties(m_qsMid.toStdString().c_str());
    rasHi.setProperties(m_qsHi.toStdString().c_str());
    rasLo.setProperties(m_qsLoPond.toStdString().c_str());
    rasMid.setProperties(m_qsMidPond.toStdString().c_str());
    rasHi.setProperties(m_qsHiPond.toStdString().c_str());
    int nPonds = pPolys->GetFeatureCount();
    int pondID;
    double lo, mid, hi, gelev, demVal, depValNew;

    for (int i=0; i<nPonds; i++)
    {
        OGRFeature *pDamFeat, *pPolyFeat;
        pPolyFeat = pPolys->GetFeature(i);
        pPolygon = (OGRPolygon*) pPolyFeat->GetGeometryRef();
        pRing = pPolygon->getExteriorRing();
        pRing->getEnvelope(&bound);
        pondID = pPolyFeat->GetFieldAsInteger("pond_ID");
        pDamFeat = pPts->GetFeature(pondID);
        gelev = pDamFeat->GetFieldAsDouble("g_elev");
        lo = gelev + pDamFeat->GetFieldAsDouble("ht_lo");
        mid = gelev + pDamFeat->GetFieldAsDouble("ht_mid");
        hi = gelev + pDamFeat->GetFieldAsDouble("ht_hi");
        qDebug()<<pondID<<" of "<<nPonds;
        int left = rasDem.getCol(bound.MinX)-1;
        int right = rasDem.getCol(bound.MaxX)+1;;
        int top = rasDem.getRow(bound.MaxY)-1;
        int bottom = rasDem.getRow(bound.MinY)+1;

        for (int j=top; j<bottom; j++)
        {
            for (int k=left; k<right; k++)
            {
                double x = rasDem.xCoordinate(k);
                double y = rasDem.yCoordinate(j);

                if (Geometry::pointInPolygon(pRing, x, y))
                {
                    demVal = rasDem.valueAtPoint(x,y);
                    if (demVal < lo)
                    {
                        depValNew = lo - demVal;
                        if (depValNew > rasLo.valueAtPoint(x,y) && depValNew > 0.0)
                        {
                            rasLo.writeCellValue(m_qsLo.toStdString().c_str(), x, y, depValNew);
                            loPond.writeCellValue(m_qsLoPond.toStdString().c_str(), x, y, pondID*1.0);
                        }

                    }
                    if (demVal < mid)
                    {
                        depValNew = mid - demVal;
                        if (depValNew > rasMid.valueAtPoint(x,y) && depValNew > 0.0)
                        {
                            rasMid.writeCellValue(m_qsMid.toStdString().c_str(), x, y, depValNew);
                            midPond.writeCellValue(m_qsMidPond.toStdString().c_str(), x, y, pondID*1.0);
                        }
                    }
                    if (demVal < hi)
                    {
                        depValNew = hi - demVal;
                        if (depValNew > rasHi.valueAtPoint(x,y) && depValNew > 0.0)
                        {
                            rasHi.writeCellValue(m_qsHi.toStdString().c_str(), x, y, depValNew);
                            hiPond.writeCellValue(m_qsHiPond.toStdString().c_str(), x, y, pondID*1.0);
                        }
                    }
                }
            }
        }

        OGRFeature::DestroyFeature(pDamFeat);
        OGRFeature::DestroyFeature(pPolyFeat);
    }
}

void DamPolygons::createDepthRasters()
{
    setRasterPaths();
    loadDriver_GDAL();
    GDALDataset *pLoPondDS, *pMidPondDS, *pHiPondDS, *pLoDS, *pMidDS, *pHiDS, *pDemDS;
    double geot[6];

    pDemDS = (GDALDataset*) GDALOpen(m_demPath, GA_ReadOnly);
    pDemDS->GetGeoTransform(geot);

    m_cellWidth = geot[1];
    m_cellHeight = fabs(geot[5]);

    pLoDS = m_pDriverTiff->Create(m_qsLo.toStdString().c_str(), pDemDS->GetRasterXSize(), pDemDS->GetRasterYSize(), 1, GDT_Float32, NULL);
    pLoDS->SetGeoTransform(geot);
    pLoDS->GetRasterBand(1)->Fill(-9999.0);
    pLoDS->GetRasterBand(1)->SetNoDataValue(-9999.0);
    GDALClose(pLoDS);

    pMidDS = m_pDriverTiff->Create(m_qsMid.toStdString().c_str(), pDemDS->GetRasterXSize(), pDemDS->GetRasterYSize(), 1, GDT_Float32, NULL);
    pMidDS->SetGeoTransform(geot);
    pMidDS->GetRasterBand(1)->Fill(-9999.0);
    pMidDS->GetRasterBand(1)->SetNoDataValue(-9999.0);
    GDALClose(pMidDS);

    pHiDS = m_pDriverTiff->Create(m_qsHi.toStdString().c_str(), pDemDS->GetRasterXSize(), pDemDS->GetRasterYSize(), 1, GDT_Float32, NULL);
    pHiDS->SetGeoTransform(geot);
    pHiDS->GetRasterBand(1)->Fill(-9999.0);
    pHiDS->GetRasterBand(1)->SetNoDataValue(-9999.0);
    GDALClose(pHiDS);

    pLoPondDS = m_pDriverTiff->Create(m_qsLoPond.toStdString().c_str(), pDemDS->GetRasterXSize(), pDemDS->GetRasterYSize(), 1, GDT_Float32, NULL);
    pLoPondDS->SetGeoTransform(geot);
    pLoPondDS->GetRasterBand(1)->Fill(-9999.0);
    pLoPondDS->GetRasterBand(1)->SetNoDataValue(-9999.0);
    GDALClose(pLoPondDS);

    pMidPondDS = m_pDriverTiff->Create(m_qsMidPond.toStdString().c_str(), pDemDS->GetRasterXSize(), pDemDS->GetRasterYSize(), 1, GDT_Float32, NULL);
    pMidPondDS->SetGeoTransform(geot);
    pMidPondDS->GetRasterBand(1)->Fill(-9999.0);
    pMidPondDS->GetRasterBand(1)->SetNoDataValue(-9999.0);
    GDALClose(pMidPondDS);

    pHiPondDS = m_pDriverTiff->Create(m_qsHiPond.toStdString().c_str(), pDemDS->GetRasterXSize(), pDemDS->GetRasterYSize(), 1, GDT_Float32, NULL);
    pHiPondDS->SetGeoTransform(geot);
    pHiPondDS->GetRasterBand(1)->Fill(-9999.0);
    pHiPondDS->GetRasterBand(1)->SetNoDataValue(-9999.0);
    GDALClose(pHiPondDS);

    GDALClose(pDemDS);
}

void DamPolygons::createFields(OGRLayer *pLayer)
{
    OGRFieldDefn field("brat_ID", OFTInteger);
    pLayer->CreateField(&field);
    field.SetName("pond_ID");
    field.SetType(OFTInteger);
    pLayer->CreateField(&field);
}

void DamPolygons::createFields_BRAT(OGRLayer *pLayer)
{
    OGRFieldDefn field("area_lo", OFTInteger);
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

void DamPolygons::createPondPolygons(OGRLayer *pPts, OGRLayer *pPolys)
{
    double xCoords[5], yCoords[5];
    double azimuthCurrent, azimuthStart, distance, slope, maxHeight;
    int nDams = pPts->GetFeatureCount();
    createFields(pPolys);

    for (int i=0; i<nDams; i++)
    {
        OGRFeature *pDamFeature;
        OGRPolygon polygon;
        OGRLinearRing ring;
        OGRPoint point, *damPoint;
        OGRFeature *pPolyFeature = OGRFeature::CreateFeature(pPolys->GetLayerDefn());

        pDamFeature = pPts->GetFeature(i);
        setFieldValues(pPolyFeature, pDamFeature->GetFieldAsInteger("brat_ID"), i);
        azimuthStart = pDamFeature->GetFieldAsDouble("az_us");
        slope = pDamFeature->GetFieldAsDouble("slope");
        maxHeight = pDamFeature->GetFieldAsDouble("ht_max");
        distance = getUpstreamDistance(maxHeight/slope);
        damPoint = (OGRPoint*) pDamFeature->GetGeometryRef();

        for (int j=0; j<5; j++)
        {
            azimuthCurrent = Geometry::addDegrees(azimuthStart, ANGLE_OFFSET[j]);
            Geometry::calcCoords(damPoint->getX(), damPoint->getY(), azimuthCurrent, distance, xCoords[j], yCoords[j]);

            point.setX(xCoords[j]);
            point.setY(yCoords[j]);

            ring.addPoint(&point);
        }

        ring.addPoint(xCoords[0], yCoords[0]);
        polygon.addRing(&ring);
        pPolyFeature->SetGeometry(&polygon);
        pPolys->CreateFeature(pPolyFeature);
        OGRFeature::DestroyFeature(pPolyFeature);
        OGRFeature::DestroyFeature(pDamFeature);
    }
}

double DamPolygons::getUpstreamDistance(double dist)
{
    if (dist > m_maxDist)
    {
        return m_maxDist;
    }
    else
    {
        return dist;
    }
}

void DamPolygons::setDemPath(const char *demPath)
{
    m_demPath = demPath;
}

void DamPolygons::setFieldValues(OGRFeature *pFeat, int bratID, int pondID)
{
    pFeat->SetField("brat_ID", bratID);
    pFeat->SetField("pond_ID", pondID);
}

void DamPolygons::setMaxDistance(double distance)
{
    m_maxDist = distance;
}

void DamPolygons::setOutDirPath(const char *outDir)
{
    m_outDir = outDir;
}

void DamPolygons::setRasterPaths()
{
    QString dirPath = QString::fromUtf8(m_outDir);
    m_qsLo = dirPath + "/depLo.tif";
    m_qsMid = dirPath + "/depMid.tif";
    m_qsHi = dirPath + "/depHi.tif";
    m_qsLoPond = dirPath + "/depLoPond.tif";
    m_qsMidPond = dirPath + "/depMidPond.tif";
    m_qsHiPond = dirPath + "/depHiPond.tif";
}

void DamPolygons::summarizePondDepths(OGRLayer *pPts)
{
    int nDams = pPts->GetFeatureCount();
    QVector<double> areasLo(nDams), areasMid(nDams), areasHi(nDams), volumesLo(nDams, 0.0), volumesMid(nDams, 0.0), volumesHi(nDams, 0.0);
    QVector<int> pondIDs(nDams), countsLo(nDams, 0), countsMid(nDams, 0), countsHi(nDams, 0);

    for (int i=0; i<nDams; i++)
    {
        OGRFeature *pDamFeature;
        pDamFeature = pPts->GetFeature(i);
        pondIDs[i] = pDamFeature->GetFID();
        OGRFeature::DestroyFeature(pDamFeature);
    }

    loadDriver_GDAL();
    GDALDataset *pLo, *pMid, *pHi, *pLoPond, *pMidPond, *pHiPond;
    pLo = (GDALDataset*) GDALOpen(m_qsLo.toStdString().c_str(), GA_ReadOnly);
    pMid = (GDALDataset*) GDALOpen(m_qsMid.toStdString().c_str(), GA_ReadOnly);
    pHi = (GDALDataset*) GDALOpen(m_qsHi.toStdString().c_str(), GA_ReadOnly);
    pLoPond = (GDALDataset*) GDALOpen(m_qsLoPond.toStdString().c_str(), GA_ReadOnly);
    pMidPond = (GDALDataset*) GDALOpen(m_qsMidPond.toStdString().c_str(), GA_ReadOnly);
    pHiPond = (GDALDataset*) GDALOpen(m_qsHiPond.toStdString().c_str(), GA_ReadOnly);

    float *loVals = (float*) CPLMalloc(sizeof(float) * pLo->GetRasterXSize());
    float *midVals = (float*) CPLMalloc(sizeof(float) * pLo->GetRasterXSize());
    float *hiVals = (float*) CPLMalloc(sizeof(float) * pLo->GetRasterXSize());
    float *loPondVals = (float*) CPLMalloc(sizeof(float) * pLo->GetRasterXSize());
    float *midPondVals = (float*) CPLMalloc(sizeof(float) * pLo->GetRasterXSize());
    float *hiPondVals = (float*) CPLMalloc(sizeof(float) * pLo->GetRasterXSize());

    int index = -1;

    for (int i=0; i<pLo->GetRasterYSize(); i++)
    {
        pLo->GetRasterBand(1)->RasterIO(GF_Read, 0, i, pLo->GetRasterXSize(), 1, loVals, pLo->GetRasterXSize(), 1, GDT_Float32, 0, 0);
        pMid->GetRasterBand(1)->RasterIO(GF_Read, 0, i, pLo->GetRasterXSize(), 1, midVals, pLo->GetRasterXSize(), 1, GDT_Float32, 0, 0);
        pHi->GetRasterBand(1)->RasterIO(GF_Read, 0, i, pLo->GetRasterXSize(), 1, hiVals, pLo->GetRasterXSize(), 1, GDT_Float32, 0, 0);
        pLoPond->GetRasterBand(1)->RasterIO(GF_Read, 0, i, pLo->GetRasterXSize(), 1, loPondVals, pLo->GetRasterXSize(), 1, GDT_Float32, 0, 0);
        pMidPond->GetRasterBand(1)->RasterIO(GF_Read, 0, i, pLo->GetRasterXSize(), 1, midPondVals, pLo->GetRasterXSize(), 1, GDT_Float32, 0, 0);
        pHiPond->GetRasterBand(1)->RasterIO(GF_Read, 0, i, pLo->GetRasterXSize(), 1, hiPondVals, pLo->GetRasterXSize(), 1, GDT_Float32, 0, 0);
        for (int j=0; j<pLo->GetRasterXSize(); j++)
        {
            if (loVals[j] > 0.0)
            {
                index = pondIDs.indexOf(round(loPondVals[j]));
                countsLo[index]++;
                volumesLo[index] += loVals[j];
            }
            if (midVals[j] > 0.0)
            {
                index = pondIDs.indexOf(round(midPondVals[j]));
                countsMid[index]++;
                volumesMid[index] += midVals[j];
            }
            if (hiVals[j] > 0.0)
            {
                index = pondIDs.indexOf(round(hiPondVals[j]));
                countsHi[index]++;
                volumesHi[index] += hiVals[j];
            }
        }
    }

    for (int i=0; i<pondIDs.length(); i++)
    {
        areasLo[i] = countsLo[i]*m_cellWidth*m_cellHeight;
        areasMid[i] = countsMid[i]*m_cellWidth*m_cellHeight;
        areasHi[i] = countsHi[i]*m_cellWidth*m_cellHeight;
        volumesLo[i] = volumesLo[i]*m_cellWidth*m_cellHeight;
        volumesMid[i] = volumesMid[i]*m_cellWidth*m_cellHeight;
        volumesHi[i] = volumesHi[i]*m_cellWidth*m_cellHeight;
        OGRFeature *pDamFeature;
        pDamFeature = pPts->GetFeature(pondIDs[i]);
        //qDebug()<<pondIDs[i]<<areasLo[i]<<areasMid[i]<<areasHi[i]<<volumesLo[i]<<volumesMid[i]<<volumesHi[i];
        DamPoints::setPondAttributes(pDamFeature, areasLo[i], areasMid[i], areasHi[i], volumesLo[i], volumesMid[i], volumesHi[i]);
        pPts->SetFeature(pDamFeature);
        OGRFeature::DestroyFeature(pDamFeature);
    }

    CPLFree(loVals);
    CPLFree(midVals);
    CPLFree(hiVals);
    CPLFree(loPondVals);
    CPLFree(midPondVals);
    CPLFree(hiPondVals);

    GDALClose(pLo);
    GDALClose(pMid);
    GDALClose(pHi);
    GDALClose(pLoPond);
    GDALClose(pMidPond);
    GDALClose(pHiPond);
}

void DamPolygons::summarizeReachDepths(OGRLayer *pPts)
{
    loadDriver();
    OGRDataSource *pInDs;
    OGRLayer *pBratLyr;

    qDebug()<<"opening data source";
    pInDs = m_pDriverShp->CreateDataSource(m_qsBratDir.toStdString().c_str(), NULL);
    qDebug()<<"data source retrieved";
    pBratLyr = pInDs->GetLayerByName(m_qsBratName.toStdString().c_str());
    qDebug()<<"brat layer retrieved";
    createFields_BRAT(pBratLyr);
    int nReaches = pBratLyr->GetFeatureCount();

    QVector<int> reachIDs(nReaches);
    QVector<double> areasLo(nReaches, 0.0), volumesLo(nReaches, 0.0), areasMid(nReaches, 0.0), volumesMid(nReaches, 0.0), areasHi(nReaches, 0.0), volumesHi(nReaches, 0.0);

    for (int i=0; i<nReaches; i++)
    {
        OGRFeature *pReachFeat = pBratLyr->GetFeature(i);
        reachIDs[i] = pReachFeat->GetFID();
        OGRFeature::DestroyFeature(pReachFeat);
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
        OGRFeature *pReachFeat = pBratLyr->GetFeature(i);
        int index = reachIDs.indexOf(pReachFeat->GetFID());
        pReachFeat->SetField("area_lo", areasLo[index]);
        pReachFeat->SetField("area_mid", areasMid[index]);
        pReachFeat->SetField("area_hi", areasHi[index]);
        pReachFeat->SetField("vol_lo", volumesLo[index]);
        pReachFeat->SetField("vol_mid", volumesMid[index]);
        pReachFeat->SetField("vol_hi", volumesHi[index]);
        pBratLyr->SetFeature(pReachFeat);
        OGRFeature::DestroyFeature(pReachFeat);
    }
    OGRDataSource::DestroyDataSource(pInDs);
}

