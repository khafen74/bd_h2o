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

    OGRDataSource *pOutDs;
    OGRLayer *pPoly_lyr, *pPts_lyr;

    pOutDs = m_pDriverShp->CreateDataSource(m_outDir, NULL);
    pPts_lyr = pOutDs->GetLayerByName(pondPts.getLayerName());
    pPoly_lyr = pOutDs->CreateLayer(m_layerName, pPts_lyr->GetSpatialRef(), wkbPolygon, NULL);

    createPondPolygons(pPts_lyr, pPoly_lyr);
    calculateWaterDepth(pPts_lyr, pPoly_lyr);

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
    Raster rasDem, rasLo, rasMid, rasHi;
    rasDem.setProperties(m_demPath);
    rasLo.setProperties(m_qsLo.toStdString().c_str());
    rasMid.setProperties(m_qsMid.toStdString().c_str());
    rasHi.setProperties(m_qsHi.toStdString().c_str());
    int nPonds = pPolys->GetFeatureCount();
    int pondID;
    double lo, mid, hi, gelev, demVal, depValNew;
    double lowVol, midVol, hiVol, lowArea, midArea, hiArea;
    int lowCount, midCount, hiCount;
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
        lowVol = 0.0, midVol = 0.0, hiVol = 0.0;
        lowCount = 0, midCount = 0, hiCount = 0;

        for (int j=top; j<bottom; j++)
        {
            for (int k=left; k<right; k++)
            {
                double x = rasDem.xCoordinate(k);
                double y = rasDem.yCoordinate(j);
                //int r = rasDem.getRow(y)
                //qDebug()<<QString::number(x, 'f', 2)<<QString::number(y, 'f', 2);

                if (Geometry::pointInPolygon(pRing, x, y))
                {
                    //qDebug()<<"true";
                    demVal = rasDem.valueAtPoint(x,y);
                    if (demVal < lo)
                    {
                        depValNew = lo - demVal;
                        if (depValNew > rasLo.valueAtPoint(x,y) && depValNew > 0.0)
                        {
                            rasLo.writeCellValue(x, y, depValNew);
                            lowCount++;
                            lowVol += depValNew;
                        }

                    }
                    if (demVal < mid)
                    {
                        depValNew = mid - demVal;
                        if (depValNew > rasMid.valueAtPoint(x,y) && depValNew > 0.0)
                        {
                            rasMid.writeCellValue(x, y, depValNew);
                            midCount++;
                            midVol += depValNew;
                        }
                    }
                    if (demVal < hi)
                    {
                        depValNew = hi - demVal;
                        if (depValNew > rasHi.valueAtPoint(x,y) && depValNew > 0.0)
                        {
                            rasHi.writeCellValue(x, y, depValNew);
                            hiCount++;
                            hiVol += depValNew;
                        }
                    }
                }
            }
        }

        DamPoints::setPondAttributes(pDamFeat, lowCount*m_cellWidth*m_cellHeight, midCount*m_cellWidth*m_cellHeight, hiCount*m_cellWidth*m_cellHeight,
                                     lowVol*m_cellWidth*m_cellHeight, midVol*m_cellWidth*m_cellHeight, hiVol*m_cellWidth*m_cellHeight);

        pPts->SetFeature(pDamFeat);
        OGRFeature::DestroyFeature(pDamFeat);
        OGRFeature::DestroyFeature(pPolyFeat);
    }
}

void DamPolygons::createDepthRasters()
{
    setRasterPaths();
    loadDriver_GDAL();
    GDALDataset *pLoDS, *pMidDS, *pHiDS, *pDemDS;
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
}

