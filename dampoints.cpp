#include "dampoints.h"

DamPoints::DamPoints(const char *pointsPath)
{
    init(pointsPath);
}

DamPoints::DamPoints(const char *demPath, const char *bratPath, const char *outDirPath, double modCap)
{
    setDemPath(demPath);
    setOutDir(outDirPath);
    setBratCapacity(modCap);
    init(bratPath);
}

void DamPoints::init(const char *bratPath)
{
    m_layerName = "ModeledDamPoints";
    loadDriver();

    QString qsBratDir, qsBratName;
    QFileInfo fi(QString::fromUtf8(bratPath));
    qsBratDir = fi.absolutePath();
    qsBratName = fi.baseName();

    OGRDataSource *pInDs, *pOutDs;
    OGRLayer *pBratIn, *pDamsOut;

    pInDs = m_pDriverShp->CreateDataSource(qsBratDir.toStdString().c_str(), NULL);
    pOutDs = m_pDriverShp->CreateDataSource(m_outDir, NULL);
    pBratIn = pInDs->GetLayerByName(qsBratName.toStdString().c_str());
    pDamsOut = pOutDs->CreateLayer(m_layerName, pBratIn->GetSpatialRef(), wkbPoint, NULL);

    createFields(pDamsOut);
    createDamPoints_BRAT(pBratIn, pDamsOut);

    OGRDataSource::DestroyDataSource(pInDs);
    OGRDataSource::DestroyDataSource(pOutDs);
}

void DamPoints::createDamPoints_BRAT(OGRLayer *pBratLyr, OGRLayer *pDamsLyr)
{
    const char *slopeField = "iGeo_Slope";
    const char *densField = "oCC_EX";
    double sampleDist = 50.0;
    double heigthSum = 0.0;
    int nDams = 0;

    OGRFeature *pBratFeat;
    OGRFeature *pDamFeat = OGRFeature::CreateFeature(pDamsLyr->GetLayerDefn());
    Raster raster_dem;
    int nFeatures = pBratLyr->GetFeatureCount();

    for (int i=0; i<nFeatures; i++)
    {
        int nDamCount = 0;
        double length, damDens, slope, spacing, damHeight, damElev, elev, azimuthStart, dist, endx, endy, end_elev;
        OGRPoint point1, point2;
        pBratFeat = pBratLyr->GetFeature(i);;
        OGRGeometry *pGeom = pBratFeat->GetGeometryRef();
        OGRLineString *pBratLine = (OGRLineString*) pGeom;
        int nPoints = pBratLine->getNumPoints();

        point1.setX(pBratLine->getX(0));
        point1.setY(pBratLine->getY(0));
        point2.setX(pBratLine->getX(nPoints-1));
        point2.setY(pBratLine->getY(nPoints-1));

        length = pBratLine->get_Length();
        damDens = pBratFeat->GetFieldAsDouble(densField);
        slope = pBratFeat->GetFieldAsDouble(slopeField);

        nDamCount = round(length * (damDens/1000.0) * m_modCap);

        if (nDamCount > 0)
        {
            spacing = length / (nDamCount*1.0);
        }
        else
        {
            spacing = 0.0;
        }
        azimuthStart = Geometry::calcAzimuth(point2.getX(), point2.getY(), point1.getX(), point1.getY());
        end_elev = raster_dem.sampleAlongLine_LowVal(m_demPath, pBratLine->getX(0), pBratLine->getY(0), azimuthStart, sampleDist, endx, endy);
        dist = (damHeight/slope) * 1.0;
        if (dist > spacing)
        {
            dist = spacing;
        }
        nDams += nDamCount;

        for (int j=0; j<nDamCount; j++)
        {
            OGRPoint damPoint;
            double pointDist = length - (spacing * (j * 1.0));
            damHeight = Random::random_lognormal(-0.0999, 0.42);
            heigthSum += damHeight;
            pBratLine->Value(pointDist, &damPoint);
            double x = damPoint.getX();
            double y = damPoint.getY();
            elev = raster_dem.sampleAlongLine_LowVal(m_demPath, damPoint.getX(), damPoint.getY(), azimuthStart, sampleDist, x, y);
            damPoint.setX(x);
            damPoint.setY(y);
            damElev = elev + damHeight;
            setFieldValues(pDamFeat, i, damElev, elev, slope, Geometry::calcAzimuth(damPoint.getX(), damPoint.getY(), endx, endy), x, y);
            pDamFeat->SetGeometry(&damPoint);
            pDamsLyr->CreateFeature(pDamFeat);
        }
    }
    OGRFeature::DestroyFeature(pDamFeat);
}

void DamPoints::createFields(OGRLayer *pLayer)
{
    OGRFieldDefn field("brat_ID", OFTInteger);
    pLayer->CreateField(&field);
    field.SetName("endx");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("endy");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("az_us");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("g_elev");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("d_elev");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("slope");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("area_low");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("area_mean");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("area_hi");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_low");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_mean");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_hi");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
}

const char *DamPoints::getDemPath()
{
    return m_demPath;
}

const char *DamPoints::getLayerName()
{
    return m_layerName;
}

const char *DamPoints::getOutDirPath()
{
    return m_outDir;
}


void DamPoints::loadDriver()
{
    OGRRegisterAll();
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    m_pDriverShp = registrar->GetDriverByName("ESRI Shapefile");
}

void DamPoints::setBratCapacity(double capacity)
{
    m_modCap = capacity;
}

void DamPoints::setDemPath(const char *demPath)
{
    m_demPath = demPath;
}

void DamPoints::setFieldValues(OGRFeature *pFeat, int bratID, double damElev, double groundElev, double slope, double azimuth, double ptX, double ptY)
{
    pFeat->SetField("brat_ID", bratID);
    pFeat->SetField("d_elev", damElev);
    pFeat->SetField("g_elev", groundElev);
    pFeat->SetField("slope", slope);
    pFeat->SetField("az_us", azimuth);
    pFeat->SetField("endx", ptX);
    pFeat->SetField("endy", ptY);
}

void DamPoints::setOutDir(const char *outDirPath)
{
    m_outDir = outDirPath;
}

