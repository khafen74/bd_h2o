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
    //createDepthRasters();

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

void DamPolygons::createDepthRasters()
{
    setRasterPaths();
    qDebug()<<"paths set";
    loadDriver_GDAL();
    qDebug()<<"drivers loaded";
    GDALDataset *pLoDS, *pMidDS, *pHiDS, *pDemDS;
    double geot[6];

    qDebug()<<m_demPath;
    pDemDS = (GDALDataset*) GDALOpen(m_demPath, GA_ReadOnly);
    pDemDS->GetGeoTransform(geot);
    qDebug()<<"creating lo"<<m_rasLoPath;

    pLoDS = m_pDriverTiff->Create(m_rasLoPath, pDemDS->GetRasterXSize(), pDemDS->GetRasterYSize(), 1, GDT_Float32, NULL);
    pLoDS->SetGeoTransform(geot);
    pLoDS->GetRasterBand(1)->Fill(0.0);
    pLoDS->GetRasterBand(1)->SetNoDataValue(-9999.0);
    GDALClose(pLoDS);
    qDebug()<<"creating mid"<<m_rasMidPath;

    pMidDS = m_pDriverTiff->Create(m_rasMidPath, pDemDS->GetRasterXSize(), pDemDS->GetRasterYSize(), 1, GDT_Float32, NULL);
    pMidDS->SetGeoTransform(geot);
    pMidDS->GetRasterBand(1)->Fill(-9999.0);
    pMidDS->GetRasterBand(1)->SetNoDataValue(-9999.0);
    GDALClose(pMidDS);
    qDebug()<<"creating hi"<<m_rasHiPath;

//    pHiDS = m_pDriverTiff->Create(m_rasHiPath, pDemDS->GetRasterXSize(), pDemDS->GetRasterYSize(), 1, GDT_Float32, NULL);
//    pHiDS->SetGeoTransform(geot);
//    pHiDS->GetRasterBand(1)->Fill(-9999.0);
//    pHiDS->GetRasterBand(1)->SetNoDataValue(-9999.0);
//    qDebug()<<"closing";


    GDALClose(pDemDS);
//    GDALClose(pMidDS);
//    GDALClose(pHiDS);
}

void DamPolygons::createFields(OGRLayer *pLayer)
{
    OGRFieldDefn field("brat_ID", OFTInteger);
    pLayer->CreateField(&field);
    field.SetName("pond_ID");
    field.SetType(OFTInteger);
    pLayer->CreateField(&field);
    field.SetName("WSE_max");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
}

void DamPolygons::createPondPolygons(OGRLayer *pPts, OGRLayer *pPolys)
{
    double xCoords[5], yCoords[5];
    double azimuthCurrent, azimuthStart, distance, slope, maxHeight;
    int nDams = pPts->GetFeatureCount();

    for (int i=0; i<nDams; i++)
    {
        OGRFeature *pDamFeature;
        OGRPolygon polygon;
        OGRLinearRing ring;
        OGRPoint point, *damPoint;
        OGRFeature *pPolyFeature = OGRFeature::CreateFeature(pPolys->GetLayerDefn());

        pDamFeature = pPts->GetFeature(i);
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

void DamPolygons::setFieldValues(OGRFeature *pFeat, int bratID, int pondID, double maxWSE)
{
    pFeat->SetField("brat_ID", bratID);
    pFeat->SetField("pond_ID", pondID);
    pFeat->SetField("WSE_max", maxWSE);
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
    QString temp;
    temp = QString::fromUtf8(m_outDir) + "/depLo.tif";
    m_rasLoPath = temp.toStdString().c_str();
    //qDebug()<<temp<<m_rasLoPath;
    temp = QString::fromUtf8(m_outDir) + "/depMid.tif";
    m_rasMidPath =  temp.toStdString().c_str();
    //qDebug()<<temp<<m_rasMidPath;
    temp = QString::fromUtf8(m_outDir) + "/depHi.tif";
    m_rasHiPath =  temp.toStdString().c_str();
    //qDebug()<<temp<<m_rasHiPath;
}

