#include "dampoints.h"

DamPoints::DamPoints(const char *demPath, const char *bratPath, const char *facPath, const char *outDirPath, double modCap)
{
    setDemPath(demPath);
    setFacPath(facPath);
    setOutDir(outDirPath);
    setBratCapacity(modCap);
    init(bratPath);
}

DamPoints::DamPoints(const char *demPath, const char *bratPath, const char *facPath, const char *outDirPath, double modCap, const char *exPath, int type)
{
    setDemPath(demPath);
    setFacPath(facPath);
    setOutDir(outDirPath);
    setBratCapacity(modCap);
    init(bratPath, exPath, type);
}

void DamPoints::init(const char *bratPath)
{
    m_layerName = "ModeledDamPoints";
    loadDriver();

    QFileInfo fi(QString::fromUtf8(bratPath));
    setBratPath(fi.absolutePath());
    setBratName(fi.baseName());

    OGRDataSource *pInDs, *pOutDs;
    OGRLayer *pBratIn, *pDamsOut;

    pInDs = m_pDriverShp->CreateDataSource(m_qsBratDir.toStdString().c_str(), NULL);
    pOutDs = m_pDriverShp->CreateDataSource(m_outDir, NULL);
    pBratIn = pInDs->GetLayerByName(m_qsBratName.toStdString().c_str());
    pDamsOut = pOutDs->CreateLayer(m_layerName, pBratIn->GetSpatialRef(), wkbPoint, NULL);

    createFields(pDamsOut);
    createDamPoints_BRAT(pBratIn, pDamsOut);

    OGRDataSource::DestroyDataSource(pInDs);
    OGRDataSource::DestroyDataSource(pOutDs);
}

void DamPoints::init(const char *bratPath, const char *exPath, int type)
{
    m_layerName = "ModeledDamPoints";
    loadDriver();

    QString qsBratDir, qsBratName;
    QFileInfo fi(QString::fromUtf8(bratPath));
    setBratPath(fi.absolutePath());
    setBratName(fi.baseName());
    QFileInfo fi2(QString::fromUtf8(exPath));
    QString exFilePath = fi2.absolutePath();
    QString exName = fi2.baseName();

    OGRDataSource *pInDs, *pOutDs, *pExDS;
    OGRLayer *pBratIn, *pDamsOut, *pDamsIn;

    pInDs = m_pDriverShp->CreateDataSource(m_qsBratDir.toStdString().c_str(), NULL);
    pExDS = m_pDriverShp->CreateDataSource(exFilePath.toStdString().c_str(), NULL);
    pOutDs = m_pDriverShp->CreateDataSource(m_outDir, NULL);
    pBratIn = pInDs->GetLayerByName(m_qsBratName.toStdString().c_str());
    pDamsIn = pExDS->GetLayerByName(exName.toStdString().c_str());
    pDamsOut = pOutDs->CreateLayer(m_layerName, pBratIn->GetSpatialRef(), wkbPoint, NULL);

    createFields(pDamsOut);

    if (type == 1)
    {
        createDamPoints_Copy(pBratIn, pDamsOut, pDamsIn);
    }
    else if (type == 2)
    {
        createDamPoints_Heights(pBratIn, pDamsOut, pDamsIn);
    }

    OGRDataSource::DestroyDataSource(pInDs);
    OGRDataSource::DestroyDataSource(pOutDs);
}

void DamPoints::compareArea(const char *damsIn, const char *csvOut)
{
    QFileInfo fi2(QString::fromUtf8(damsIn));
    QString exFilePath = fi2.absolutePath();
    QString exName = fi2.baseName();

    OGRDataSource *pModDS, *pCompDS;
    OGRLayer *pModLyr, *pCompLyr;

    pModDS = m_pDriverShp->CreateDataSource(m_outDir, NULL);
    pCompDS = m_pDriverShp->CreateDataSource(exFilePath.toStdString().c_str(), NULL);
    pModLyr = pModDS->GetLayerByName(m_layerName);
    pCompLyr = pCompDS->GetLayerByName(exName.toStdString().c_str());

    QFile file(QString::fromUtf8(csvOut));
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outStream(&file);
    outStream << "fid,area,lo,mid,hi,dist,slope\n";
    //outStream << "fid,type,area,lo,hi\n";

    int nFeats = pCompLyr->GetFeatureCount();

    OGRFeature *pModFeat, *pCompFeat;
    double area, lo, mid, hi, dist, slope;

    for (int i=0; i<nFeats; i++)
    {
        pCompFeat = pCompLyr->GetFeature(i);
        pModFeat = pModLyr->GetFeature(i);
        OGRPoint *point1, *point2;
        OGRGeometry *pGeom = pCompFeat->GetGeometryRef();
        point1 = (OGRPoint*) pGeom;
        pGeom = pModFeat->GetGeometryRef();
        point2 = (OGRPoint*) pGeom;
        dist = Geometry::distance_point(point1->getX(), point1->getY(), point2->getX(), point2->getY());
        area = pCompFeat->GetFieldAsDouble("Shape_Area");
        slope = pCompFeat->GetFieldAsDouble("iGeo_Slope");
        lo = pModFeat->GetFieldAsDouble("area_lo");
        mid = pModFeat->GetFieldAsDouble("area_mid");
        hi = pModFeat->GetFieldAsDouble("area_hi");
        outStream <<i<<","<<area<<","<<lo<<","<<mid<<","<<hi<<","<<dist<<","<<slope<<"\n";
        //outStream <<i<<","<<"Observed"<<","<<area<<"\n";
        //outStream <<i<<","<<"Modeled"<<","<<mid<<","<<lo<<","<<hi<<"\n";
    }

    file.close();
    OGRFeature::DestroyFeature(pModFeat);
    OGRFeature::DestroyFeature(pCompFeat);
    OGRDataSource::DestroyDataSource(pModDS);
    OGRDataSource::DestroyDataSource(pCompDS);
}

void DamPoints::createDamPoints_BRAT(OGRLayer *pBratLyr, OGRLayer *pDamsLyr)
{
    const char *slopeField = "iGeo_Slope";
    const char *densField = "oCC_EX";
    double sampleDist = 50.0;
    int nDams = 0;
    long lastFID = -9999;

    OGRFeature *pBratFeat;
    OGRFeature *pDamFeat = OGRFeature::CreateFeature(pDamsLyr->GetLayerDefn());
    Raster raster_dem;
    int nFeatures = pBratLyr->GetFeatureCount();

    for (int i=0; i<nFeatures; i++)
    {
        //qDebug()<<"BRAT reach "<<i<<" of "<<nFeatures;
        int nDamCount = 0, nErrCount = 0;
        double length, damDens, slope, spacing, elev, azimuthStart, endx, endy, end_elev;
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
        qDebug()<<nDams<<damDens/1000.0<<nDamCount<<m_modCap<<length;

        if (nDamCount > 0)
        {
            spacing = length / (nDamCount*1.0);
        }
        else
        {
            spacing = 0.0;
        }
        azimuthStart = Geometry::calcAzimuth(point2.getX(), point2.getY(), point1.getX(), point1.getY());
        //end_elev = raster_dem.sampleAlongLine_LowVal(m_demPath, pBratLine->getX(0), pBratLine->getY(0), azimuthStart, sampleDist, endx, endy);
        end_elev = raster_dem.sampleAlongLine_RasterVal(m_demPath, m_facPath, pBratLine->getX(0), pBratLine->getY(0), azimuthStart, sampleDist, endx, endy);
        nDams += nDamCount;

        for (int j=0; j<nDamCount; j++)
        {
            OGRPoint damPoint;
            double pointDist = length - (spacing * (j * 1.0));
            Statistics lognormal(Random::randomSeries(1000, RDT_lnorm, -0.09, 0.42), RDT_lnorm);
            lognormal.calcCredibleInterval(CI_95);
            pBratLine->Value(pointDist, &damPoint);
            double x = damPoint.getX();
            double y = damPoint.getY();
            elev = raster_dem.sampleAlongLine_LowVal(m_demPath, damPoint.getX(), damPoint.getY(), azimuthStart, sampleDist, x, y);
            //elev = raster_dem.sampleAlongLine_RasterVal(m_demPath, m_facPath, damPoint.getX(), damPoint.getY(), azimuthStart, sampleDist, x, y);
            damPoint.setX(x);
            damPoint.setY(y);
            setFieldValues(pDamFeat, i, elev, slope, Geometry::calcAzimuth(damPoint.getX(), damPoint.getY(), endx, endy), x, y);
            setDamHeights(pDamFeat, lognormal.getQuantile(0.025), lognormal.getQuantile(0.5), lognormal.getQuantile(0.975), VectorOps::max(lognormal.getData()));

            pDamFeat->SetGeometry(&damPoint);
            qDebug()<<"Field Values "<<i<<elev<<slope<<Geometry::calcAzimuth(damPoint.getX(), damPoint.getY(), endx, endy)<<x<<y;
            if (elev > 0.0)
            {
                if (pDamsLyr->CreateFeature(pDamFeat) != OGRERR_NONE)
                {
                    nErrCount++;
                    qDebug()<<"error writing dam point "<<nErrCount;
                }
            }
        }
    }
    OGRFeature::DestroyFeature(pDamFeat);
    OGRFeature::DestroyFeature(pBratFeat);
}

void DamPoints::createDamPoints_Copy(OGRLayer *pBratLyr, OGRLayer *pDamsLyr, OGRLayer *pExLyr)
{
    const char *slopeField = "iGeo_Slope";
    double sampleDist = 50.0;
    OGRFeature *pBratFeat, *pOldFeat;
    OGRFeature *pDamFeat = OGRFeature::CreateFeature(pDamsLyr->GetLayerDefn());
    Raster raster_dem;
    int nPrimary = 0, nSecondary = 0;
    int nFeatures = pExLyr->GetFeatureCount();

    for (int i=0; i<nFeatures; i++)
    {
        qDebug()<<"creating dam"<<i+1<<"of"<<nFeatures;
        pOldFeat = pExLyr->GetFeature(i);
        int nBratFID = pOldFeat->GetFieldAsInteger("ID");
        pBratFeat = pBratLyr->GetFeature(nBratFID);
        double slope = pBratFeat->GetFieldAsDouble(slopeField);
        OGRGeometry *pGeom = pBratFeat->GetGeometryRef();
        OGRLineString *pBratLine = (OGRLineString*) pGeom;
        pGeom = pOldFeat->GetGeometryRef();
        OGRPoint *pOldDam = (OGRPoint*) pGeom;
        double az = Geometry::calcAzimuth(pBratLine->getX(pBratLine->getNumPoints()-1), pBratLine->getY(pBratLine->getNumPoints()-1), pBratLine->getX(0), pBratLine->getY(0));
        double rnum = ((double) rand() / (RAND_MAX));
        double dht;
        Statistics normDist(Random::randomSeries(1000, RDT_norm, 0.93, 0.17), RDT_norm);
        dht = Random::random_normal(0.92, 0.17);
        if (rnum <= 0.15)
        {
            normDist.setSample(Random::randomSeries(1000, RDT_norm, 1.16, 0.20));
            dht = Random::random_normal(1.16, 0.2);
            nPrimary++;
            //qDebug()<<"primary"<<rnum;
        }
        else
        {
            nSecondary++;
            //qDebug()<<"secondary"<<rnum;
        }

        double x = pOldDam->getX();
        double y = pOldDam->getY();
        //double elev = raster_dem.sampleAlongLine_LowVal(m_demPath, pOldDam->getX(), pOldDam->getY(), az, sampleDist, x, y);
        double elev = raster_dem.sampleAlongLine_RasterVal(m_demPath, m_facPath, pOldDam->getX(), pOldDam->getY(), az, sampleDist, x, y);
        pOldDam->setX(x);
        pOldDam->setY(y);
        setFieldValues(pDamFeat, i, elev, slope, Geometry::calcAzimuth(pOldDam->getX(), pOldDam->getY(), pBratLine->getX(0), pBratLine->getY(0)), x, y);
        double loHt, midHt, hiHt, maxHt;
        loHt = normDist.getQuantile(0.025), midHt = normDist.getQuantile(0.5), hiHt = normDist.getQuantile(0.975), maxHt = VectorOps::max(normDist.getData());
        //setDamHeights(pDamFeat, dht*dht, dht*dht, dht*dht, dht*dht);
        setDamHeights(pDamFeat, loHt*loHt, midHt*midHt, hiHt*hiHt, maxHt*maxHt);
        //setDamHeights(pDamFeat, normDist.getQuantile(0.025), normDist.getQuantile(0.5), normDist.getQuantile(0.975), VectorOps::max(normDist.getData()));
        //qDebug()<<loHt*loHt<<midHt*midHt<<hiHt*hiHt;
        //qDebug()<<"brat id"<<nBratFID;

        pDamFeat->SetGeometry(pOldDam);
        if (elev > 0.0)
        {
            pDamsLyr->CreateFeature(pDamFeat);
        }
    }
    qDebug()<<"primary"<<nPrimary<<"secondary"<<nSecondary;
    OGRFeature::DestroyFeature(pDamFeat);
    OGRFeature::DestroyFeature(pBratFeat);
    OGRFeature::DestroyFeature(pOldFeat);
}

void DamPoints::createDamPoints_Heights(OGRLayer *pBratLyr, OGRLayer *pDamsLyr, OGRLayer *pExLyr)
{
    const char *slopeField = "iGeo_Slope";
    double sampleDist = 50.0;
    OGRFeature *pBratFeat, *pOldFeat;
    OGRFeature *pDamFeat = OGRFeature::CreateFeature(pDamsLyr->GetLayerDefn());
    Raster raster_dem;
    int nFeatures = pExLyr->GetFeatureCount();

    for (int i=0; i<nFeatures; i++)
    {
        pOldFeat = pExLyr->GetFeature(i);
        int nBratFID = pOldFeat->GetFieldAsInteger("ID");
        double damHeight = pOldFeat->GetFieldAsDouble("Dam_Height");
        pBratFeat = pBratLyr->GetFeature(nBratFID);
        double slope = pBratFeat->GetFieldAsDouble(slopeField);
        OGRGeometry *pGeom = pBratFeat->GetGeometryRef();
        OGRLineString *pBratLine = (OGRLineString*) pGeom;
        pGeom = pOldFeat->GetGeometryRef();
        OGRPoint *pOldDam = (OGRPoint*) pGeom;
        double az = Geometry::calcAzimuth(pBratLine->getX(pBratLine->getNumPoints()-1), pBratLine->getY(pBratLine->getNumPoints()-1), pBratLine->getX(0), pBratLine->getY(0));
        double x = pOldDam->getX();
        double y = pOldDam->getY();
        //double elev = raster_dem.sampleAlongLine_LowVal(m_demPath, pOldDam->getX(), pOldDam->getY(), az, sampleDist, x, y);
        double elev = raster_dem.sampleAlongLine_RasterVal(m_demPath, m_facPath, pOldDam->getX(), pOldDam->getY(), az, sampleDist, x, y);
        pOldDam->setX(x);
        pOldDam->setY(y);
        setFieldValues(pDamFeat, i, elev, slope, Geometry::calcAzimuth(pOldDam->getX(), pOldDam->getY(), pBratLine->getX(0), pBratLine->getY(0)), x, y);
        setDamHeights(pDamFeat, damHeight, damHeight, damHeight, damHeight);

        pDamFeat->SetGeometry(pOldDam);
        if (elev > 0.0)
        {
            pDamsLyr->CreateFeature(pDamFeat);
        }
    }

    OGRFeature::DestroyFeature(pDamFeat);
    OGRFeature::DestroyFeature(pBratFeat);
    OGRFeature::DestroyFeature(pOldFeat);
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
    field.SetName("ht_max");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("ht_lo");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("ht_mid");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("ht_hi");
    field.SetType(OFTReal);
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

QString DamPoints::getBratDir()
{
    return m_qsBratDir;
}

QString DamPoints::getBratName()
{
    return m_qsBratName;
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

void DamPoints::setBratPath(QString path)
{
    m_qsBratDir = path;
}

void DamPoints::setBratName(QString name)
{
    m_qsBratName = name;
}

void DamPoints::setDemPath(const char *demPath)
{
    m_demPath = demPath;
}

void DamPoints::setDamHeights(OGRFeature *pFeat, double low, double mid, double high, double max)
{
    pFeat->SetField("ht_lo", low);
    pFeat->SetField("ht_mid", mid);
    pFeat->SetField("ht_hi", high);
    pFeat->SetField("ht_max", max);
}

void DamPoints::setFacPath(const char *facPath)
{
    m_facPath = facPath;
}

void DamPoints::setFieldValues(OGRFeature *pFeat, int bratID, double groundElev, double slope, double azimuth, double ptX, double ptY)
{
    pFeat->SetField("brat_ID", bratID);
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

void DamPoints::setPondAttributes(OGRFeature *pFeat, double lowarea, double midarea, double hiarea, double lowvol, double midvol, double hivol)
{
    pFeat->SetField("area_lo", lowarea);
    pFeat->SetField("area_mid", midarea);
    pFeat->SetField("area_hi", hiarea);
    pFeat->SetField("vol_lo", lowvol);
    pFeat->SetField("vol_mid", midvol);
    pFeat->SetField("vol_hi", hivol);
}

