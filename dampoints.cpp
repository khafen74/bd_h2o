#include "dampoints.h"

//model at percentage (proportion 0-1.0) of existing brat capacity
DamPoints::DamPoints(const char *demPath, const char *bratPath, const char *facPath, const char *outDirPath, double modCap)
{
    m_statPath = "";
    setDemPath(demPath);
    setFacPath(facPath);
    setOutDir(outDirPath);
    setBratCapacity(modCap);
    qDebug()<<"init points";
    init(bratPath);
}

DamPoints::DamPoints(const char *demPath, const char *bratPath, const char *facPath, const char *statPath, const char *outDirPath, double modCap)
{
    m_statPath = statPath;
    m_nDamPlace = 1;
    setDemPath(demPath);
    setFacPath(facPath);
    setOutDir(outDirPath);
    setBratCapacity(modCap);
    qDebug()<<"init points";
    init(bratPath);
}

DamPoints::DamPoints(const char *demPath, const char *bratPath, const char *facPath, const char *statPath, const char *outDirPath, double modCap, int nDamPlaceType)
{
    m_statPath = statPath;
    m_nDamPlace = nDamPlaceType;
    setDemPath(demPath);
    setFacPath(facPath);
    setOutDir(outDirPath);
    setBratCapacity(modCap);
    qDebug()<<"init points";
    init(bratPath);
}

//model locations of existing dams (exPath), type determines if dams are moved to flow accumulation (1) or not moved (2)
DamPoints::DamPoints(const char *demPath, const char *bratPath, const char *facPath, const char *statPath, const char *outDirPath, double modCap, const char *exPath, int type)
{
    qDebug()<<"running from points";
    m_statPath = statPath;
    setDemPath(demPath);
    setFacPath(facPath);
    setOutDir(outDirPath);
    setBratCapacity(modCap);
    init(bratPath, exPath, type);
}

void DamPoints::init(const char *bratPath)
{
    qDebug()<<bratPath;
    qDebug()<<m_outDir;
    m_layerName = "ModeledDamPoints";
    loadDriver();

    //path and name of BRAT shapefile
    QFileInfo fi(QString::fromUtf8(bratPath));
    setBratPath(fi.absolutePath());
    setBratName(fi.baseName());

    OGRDataSource *pInDs, *pOutDs;
    OGRLayer *pBratIn, *pDamsOut;

    //load BRAT shapefile and create output shapefile for dams
    qDebug()<<"creating in dir"<<m_qsBratDir;
    pInDs = m_pDriverShp->CreateDataSource(m_qsBratDir.toStdString().c_str(), NULL);
    qDebug()<<"creating out dir"<<m_outDir;
    pOutDs = m_pDriverShp->CreateDataSource(m_outDir, NULL);
    pBratIn = pInDs->GetLayerByName(m_qsBratName.toStdString().c_str());
    pDamsOut = pOutDs->CreateLayer(m_layerName, pBratIn->GetSpatialRef(), wkbPoint, NULL);

    //create fields for modeled dam points
    createFields(pDamsOut);
    //create modeled dam points based on BRAT estimates
    if (m_nDamPlace == 1)
    {
        createDamPoints_BRAT(pBratIn, pDamsOut);
    }
    else if (m_nDamPlace == 2)
    {
        sortByCapacity(pBratIn);
        createDamPoints_BRATcomplex(pBratIn, pDamsOut);
    }
    else if (m_nDamPlace == 3)
    {
        sortByCapacity(pBratIn);
        createDamPoints_BRATcomplex100(pBratIn, pDamsOut);
    }

    OGRDataSource::DestroyDataSource(pInDs);
    OGRDataSource::DestroyDataSource(pOutDs);
}

void DamPoints::init(const char *bratPath, const char *exPath, int type)
{
    qDebug()<<"initializing";
    m_layerName = "ModeledDamPoints";
    loadDriver();

    //path and name of BRAT linework
    QFileInfo fi(QString::fromUtf8(bratPath));
    setBratPath(fi.absolutePath());
    setBratName(fi.baseName());
    //path and name for existing dam loations
    QFileInfo fi2(QString::fromUtf8(exPath));
    QString exFilePath = fi2.absolutePath();
    QString exName = fi2.baseName();
    qDebug()<<exFilePath;
    qDebug()<<exName;

    OGRDataSource *pInDs, *pOutDs, *pExDS;
    OGRLayer *pBratIn, *pDamsOut, *pDamsIn;

    //input shapefiles (located in BRAT directory)
    pInDs = OGRSFDriverRegistrar::Open(m_qsBratDir.toStdString().c_str(), 1);
    //input shapefiles (located in input dams directory)
    pExDS = m_pDriverShp->CreateDataSource(exFilePath.toStdString().c_str(), NULL);
    //output shapefiles (directory specified as input
    pOutDs = m_pDriverShp->CreateDataSource(m_outDir, NULL);
    //get existing BRAT and dams layers
    pBratIn = pInDs->GetLayerByName(m_qsBratName.toStdString().c_str());
    pDamsIn = pExDS->GetLayerByName(exName.toStdString().c_str());
    //create output shapefile of modeled dams
    pDamsOut = pOutDs->CreateLayer(m_layerName, pBratIn->GetSpatialRef(), wkbPoint, NULL);

    createFields(pDamsOut);

    if (type == 1)
    {
        //use existing dam locations (copy created)
        qDebug()<<"creating dam points";
        createDamPoints_Copy(pBratIn, pDamsOut, pDamsIn);
    }
    else if (type == 2)
    {
        //use existing dam points, maintain locations (do not move to flow accumulation)
        createDamPoints_CopyLoc(pBratIn, pDamsOut, pDamsIn);
    }
    else if (type == 3)
    {
        //use existing dam points with heights
        createDamPoints_Heights(pBratIn, pDamsOut, pDamsIn);
    }
    else if (type == 4)
    {
        //use existing dam points with heights, maintain locations (do not move to flow accumulation)
        createDamPoints_HeightsLoc(pBratIn, pDamsOut, pDamsIn);
    }

    OGRDataSource::DestroyDataSource(pInDs);
    OGRDataSource::DestroyDataSource(pOutDs);
}

//writes csv of modeled area, currently not used, from previous workflow
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
    qDebug()<<"starting dam points from BRAT";
    const char *slopeField = "iGeo_Slope";
    const char *densField = "oCC_EX";
    double sampleDist = 50.0;
    int nDams = 0, nPrimary = 0, nSecondary = 0;
    long lastFID = -9999;

    //BRAT line segment
    OGRFeature *pBratFeat;
    //feature for modeled dams layer
    qDebug()<<"loading ogr feature";
    OGRFeature *pDamFeat = OGRFeature::CreateFeature(pDamsLyr->GetLayerDefn());
    qDebug()<<"done";
    Raster raster_dem;
    qDebug()<<"counting features";
    int nFeatures = pBratLyr->GetFeatureCount();
    qDebug()<<"features"<<nFeatures;

    //loop through all features in BRAT layer
    for (int i=0; i<nFeatures; i++)
    {
        //qDebug()<<"BRAT reach "<<i<<" of "<<nFeatures;
        int nDamCount = 0, nErrCount = 0;
        double length, damDens, slope, spacing, elev, azimuthStart, endx, endy, end_elev;
        OGRPoint point1, point2;
        //get BRAT feature
        pBratFeat = pBratLyr->GetFeature(i);;
        OGRGeometry *pGeom = pBratFeat->GetGeometryRef();
        //geometry of BRAT feature
        OGRLineString *pBratLine = (OGRLineString*) pGeom;
        int nPoints = pBratLine->getNumPoints();

        //start and end points of BRAT segment, used to calculate azimuth
        point1.setX(pBratLine->getX(0));
        point1.setY(pBratLine->getY(0));
        point2.setX(pBratLine->getX(nPoints-1));
        point2.setY(pBratLine->getY(nPoints-1));

        //length of BRAT segment
        length = pBratLine->get_Length();
        //BRAT existing capacity (dams/km)
        damDens = pBratFeat->GetFieldAsDouble(densField);
        //slope of BRAT segment
        slope = pBratFeat->GetFieldAsDouble(slopeField);

        //calculate number of dams from dam density, segment, length, and percent of capacity (m_modCap)
        nDamCount = ceil(length * (damDens/1000.0) * m_modCap);
        //qDebug()<<nDams<<damDens/1000.0<<nDamCount<<m_modCap<<length;

        if (nDamCount > 0)
        {
            //distance between dams
            spacing = length / (nDamCount*1.0);
        }
        else
        {
            spacing = 0.0;
        }
        //calculate azimuth (from start to end) of BRAT segment
        azimuthStart = Geometry::calcAzimuth(point2.getX(), point2.getY(), point1.getX(), point1.getY());
        //end_elev = raster_dem.sampleAlongLine_LowVal(m_demPath, pBratLine->getX(0), pBratLine->getY(0), azimuthStart, sampleDist, endx, endy);
        //find elevation on stream network closes to dam point (on line perpendicular to BRAT segment)
        end_elev = raster_dem.sampleAlongLine_RasterVal(m_demPath, m_facPath, pBratLine->getX(0), pBratLine->getY(0), azimuthStart, sampleDist, endx, endy);
        nDams += nDamCount;

        //create a point for each dam to be modeled on BRAT segment
        for (int j=0; j<nDamCount; j++)
        {
            //location of modeled dam
            OGRPoint damPoint;
            //location on BRAT segment
            double pointDist = length - (spacing * (j * 1.0));
            //random sample of 1000 heights from dam height distribtuion
            double rnum = ((double) rand() / (RAND_MAX));
            double dht;
            Statistics normDist(Random::randomSeries(1000, RDT_norm, 0.92, 0.17), RDT_norm);
            dht = Random::random_normal(0.92, 0.17);
            if (rnum <= 0.15)
            {
                normDist.setSample(Random::randomSeries(1000, RDT_norm, 1.14, 0.19));
                dht = Random::random_normal(1.14, 0.19);
                nPrimary++;
            }
            else
            {
                nSecondary++;
            }
            //set point location on BRAT segment
            pBratLine->Value(pointDist, &damPoint);
            //location of dam point
            double x = damPoint.getX();
            double y = damPoint.getY();
            //move dam point to lowest cross sectional elevation
            //elev = raster_dem.sampleAlongLine_LowVal(m_demPath, damPoint.getX(), damPoint.getY(), azimuthStart, sampleDist, x, y);
            //move dam point to stream network raster
            elev = raster_dem.sampleAlongLine_RasterVal(m_demPath, m_facPath, damPoint.getX(), damPoint.getY(), azimuthStart, sampleDist, x, y);
            //set location of point
            damPoint.setX(x);
            damPoint.setY(y);
            //set field values and dam heights
            setFieldValues(pDamFeat, i, elev, slope, Geometry::calcAzimuth(damPoint.getX(), damPoint.getY(), endx, endy), x, y);
            double loHt, midHt, hiHt, maxHt;
            loHt = normDist.getQuantile(0.025), midHt = normDist.getQuantile(0.5), hiHt = normDist.getQuantile(0.975), maxHt = VectorOps::max(normDist.getData());
            setDamHeights(pDamFeat, loHt*loHt, midHt*midHt, hiHt*hiHt, maxHt*maxHt);

            pDamFeat->SetGeometry(&damPoint);
            //qDebug()<<"Field Values "<<i<<elev<<slope<<Geometry::calcAzimuth(damPoint.getX(), damPoint.getY(), endx, endy)<<x<<y;
            //make sure point is located inside model domain (if the elevation value is NoData, outside domain)
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
    qDebug()<<nDams<<"modeled";
    OGRFeature::DestroyFeature(pDamFeat);
    OGRFeature::DestroyFeature(pBratFeat);
}

void DamPoints::createDamPoints_BRATcomplex(OGRLayer *pBratLyr, OGRLayer *pDamsLyr)
{
    qDebug()<<"starting dam points from BRAT. . . complex";
    const char *slopeField = "iGeo_Slope";
    const char *densField = "oCC_EX";
    double sampleDist = 50.0;
    long lastFID = -9999;
    double maxCap = VectorOps::sum(m_qvMaxDams);
    double modCap = maxCap * m_modCap;
    int nTotalDams = 0, nPrimary = 0, nSecondary = 0;

    //BRAT line segment
    OGRFeature *pBratFeat;
    //feature for modeled dams layer
    qDebug()<<"loading ogr feature";
    OGRFeature *pDamFeat = OGRFeature::CreateFeature(pDamsLyr->GetLayerDefn());
    qDebug()<<"done";
    Raster raster_dem;
    qDebug()<<"counting features";
    int nFeatures = pBratLyr->GetFeatureCount();
    qDebug()<<"features"<<nFeatures;

    //loop through all features in BRAT layer
    int i=0;
    bool stop = false;
    while (i<nFeatures && !stop)
    {
        //qDebug()<<"BRAT reach "<<i<<" of "<<nFeatures;
        int nDamCount = 0, nErrCount = 0;
        double length, slope, spacing, elev, azimuthStart, endx, endy, end_elev;
        OGRPoint point1, point2;
        //get BRAT feature
        pBratFeat = pBratLyr->GetFeature(m_qvBratFID[i]);
        OGRGeometry *pGeom = pBratFeat->GetGeometryRef();
        //geometry of BRAT feature
        OGRLineString *pBratLine = (OGRLineString*) pGeom;
        int nPoints = pBratLine->getNumPoints();

        //start and end points of BRAT segment, used to calculate azimuth
        point1.setX(pBratLine->getX(0));
        point1.setY(pBratLine->getY(0));
        point2.setX(pBratLine->getX(nPoints-1));
        point2.setY(pBratLine->getY(nPoints-1));

        //length of BRAT segment
        length = pBratLine->get_Length();

        //slope of BRAT segment
        slope = pBratFeat->GetFieldAsDouble(slopeField);

        //number of dams for a BRAT segment randomly selected from complex size distribution
        nDamCount = ceil(Random::random_lognormal(1.5516125, 0.7239713));
        //qDebug()<<nDams<<damDens/1000.0<<nDamCount<<m_modCap<<length;

        if (nDamCount > 0)
        {
            if (nDamCount > m_qvMaxDams[i])
            {
                nDamCount = m_qvMaxDams[i];
            }
            if ((nTotalDams + nDamCount) > modCap)
            {
                nDamCount = ceil(modCap - nTotalDams);
            }
            //distance between dams
            spacing = length / (nDamCount*1.0);
        }
        else
        {
            spacing = 0.0;
        }
        nTotalDams += nDamCount;
        //calculate azimuth (from start to end) of BRAT segment
        azimuthStart = Geometry::calcAzimuth(point2.getX(), point2.getY(), point1.getX(), point1.getY());
        //find elevation on stream network closest to dam point (on line perpendicular to BRAT segment)
        end_elev = raster_dem.sampleAlongLine_RasterVal(m_demPath, m_facPath, pBratLine->getX(0), pBratLine->getY(0), azimuthStart, sampleDist, endx, endy);

        //create a point for each dam to be modeled on BRAT segment
        for (int j=0; j<nDamCount; j++)
        {
            //location of modeled dam
            OGRPoint damPoint;
            //location on BRAT segment
            double pointDist = length - (spacing * (j * 1.0));
            //Determine if dam is primary or secondary and get dam heights from distribution
            double rnum = ((double) rand() / (RAND_MAX));
            double dht;
            //Heigth distribution for secondary dams
            Statistics normDist(Random::randomSeries(1000, RDT_norm, 0.92, 0.17), RDT_norm);
            dht = Random::random_normal(0.92, 0.17);
            //Primary dam
            if (rnum <= 0.15)
            {
                //Height distribution for primary dams
                normDist.setSample(Random::randomSeries(1000, RDT_norm, 1.14, 0.19));
                dht = Random::random_normal(1.14, 0.19);
                nPrimary++;
            }
            else
            {
                nSecondary++;
            }
            //set point location on BRAT segment
            pBratLine->Value(pointDist, &damPoint);
            //location of dam point
            double x = damPoint.getX();
            double y = damPoint.getY();
            //move dam point to lowest cross sectional elevation
            //elev = raster_dem.sampleAlongLine_LowVal(m_demPath, damPoint.getX(), damPoint.getY(), azimuthStart, sampleDist, x, y);
            //move dam point to stream network raster
            elev = raster_dem.sampleAlongLine_RasterVal(m_demPath, m_facPath, damPoint.getX(), damPoint.getY(), azimuthStart, sampleDist, x, y);
            //set location of point
            damPoint.setX(x);
            damPoint.setY(y);
            //set field values and dam heights
            setFieldValues(pDamFeat, i, elev, slope, Geometry::calcAzimuth(damPoint.getX(), damPoint.getY(), endx, endy), x, y);
            double loHt, midHt, hiHt, maxHt;
            loHt = normDist.getQuantile(0.025), midHt = normDist.getQuantile(0.5), hiHt = normDist.getQuantile(0.975), maxHt = VectorOps::max(normDist.getData());
            setDamHeights(pDamFeat, loHt*loHt, midHt*midHt, hiHt*hiHt, maxHt*maxHt);

            pDamFeat->SetGeometry(&damPoint);
            //qDebug()<<"Field Values "<<i<<elev<<slope<<Geometry::calcAzimuth(damPoint.getX(), damPoint.getY(), endx, endy)<<x<<y;
            //make sure point is located inside model domain (if the elevation value is NoData, outside domain)
            if (elev > 0.0 && end_elev > 0.0)
            {
                if (pDamsLyr->CreateFeature(pDamFeat) != OGRERR_NONE)
                {
                    nErrCount++;
                    qDebug()<<"error writing dam point "<<nErrCount;
                }
            }
        }
        if (nTotalDams >= modCap)
        {
            stop = true;
        }
        i++;
    }
    qDebug()<<"Dams Modeled"<<nTotalDams<<"Max Capacity for Scenario"<<modCap<<"Maximum Capacity"<<maxCap<<"Percent Capacity"<<m_modCap;
    qDebug()<<"primary"<<nPrimary<<"secondary"<<nSecondary;
    OGRFeature::DestroyFeature(pDamFeat);
    OGRFeature::DestroyFeature(pBratFeat);
}

void DamPoints::createDamPoints_BRATcomplex100(OGRLayer *pBratLyr, OGRLayer *pDamsLyr)
{
    qDebug()<<"starting dam points from BRAT. . . complex 100% update";
    OGRFieldDefn field("totdams", OFTInteger);

    //delete dam and complex count fields if they already exist
    if (pBratLyr->FindFieldIndex("totdams", 1)>=0)
    {
        pBratLyr->DeleteField(pBratLyr->FindFieldIndex("totdams", 1));
    }
    if (pBratLyr->FindFieldIndex("totcomp", 1)>=0)
    {
        pBratLyr->DeleteField(pBratLyr->FindFieldIndex("totcomp", 1));
    }

    //create dam and complex count fields
    pBratLyr->CreateField(&field);
    field.SetName("totcomp");
    field.SetType(OFTInteger);
    pBratLyr->CreateField(&field);
    const char *slopeField = "iGeo_Slope";
    double sampleDist = 50.0;
    double maxCap = VectorOps::sum(m_qvMaxDams);
    double modCap = maxCap * m_modCap;
    int nTotalDams = 0, nPrimary = 0, nSecondary = 0;

    //BRAT line segment
    OGRFeature *pBratFeat;
    //feature for modeled dams layer
    qDebug()<<"loading ogr feature";
    OGRFeature *pDamFeat = OGRFeature::CreateFeature(pDamsLyr->GetLayerDefn());
    qDebug()<<"done";
    Raster raster_dem;
    qDebug()<<"counting features";
    int nFeatures = pBratLyr->GetFeatureCount();
    qDebug()<<"features"<<nFeatures;

    while (nTotalDams < ceil(modCap))
    {
        //loop through all features in BRAT layer
        int i=0;
        bool stop = false;

        while (i<nFeatures && !stop)
        {
            int nDamCount = 0;
            pBratFeat = pBratLyr->GetFeature(m_qvBratFID[i]);
            int exDams = pBratFeat->GetFieldAsInteger("totdams");
            int exComp = pBratFeat->GetFieldAsInteger("totcomp");

            //number of dams for a BRAT segment randomly selected from complex size distribution
            nDamCount = ceil(Random::random_lognormal(1.5516125, 0.7239713));

            if (nDamCount > 0)
            {
                if (nDamCount > ceil(m_qvMaxDams[i]))
                {
                    //set equal to max reach capacity if complex size is larger than reach capacity
                    nDamCount = ceil(m_qvMaxDams[i]);
                }
                if ((nTotalDams + nDamCount) > modCap)
                {
                    //if maximum capacity for the area is reached reduce capacity at the final reach
                    nDamCount = ceil(modCap - nTotalDams);
                }

                pBratFeat->SetField("totdams", exDams + nDamCount);
                if(m_qvMaxDams[i] >= 1.0)
                {
                    pBratFeat->SetField("totcomp", exComp + 1);
                }
                pBratLyr->SetFeature(pBratFeat);
            }
            nTotalDams += nDamCount;
            i++;
            if (nTotalDams >= ceil(modCap))
            {
                stop = true;
            }
        }
        qDebug()<<"finished iteration of complex placement";
    }

    //loop through all features in BRAT layer
    int i=0;
    const char *damType;
    while (i<nFeatures)
    {
        //qDebug()<<"BRAT reach "<<i<<" of "<<nFeatures;
        int nErrCount = 0;
        int nDamCount, nCompCount;
        double length, slope, spacing, elev, azimuthStart, endx, endy, end_elev;
        OGRPoint point1, point2;
        //get BRAT feature
        pBratFeat = pBratLyr->GetFeature(m_qvBratFID[i]);
        OGRGeometry *pGeom = pBratFeat->GetGeometryRef();
        //geometry of BRAT feature
        OGRLineString *pBratLine = (OGRLineString*) pGeom;
        int nPoints = pBratLine->getNumPoints();

        //start and end points of BRAT segment, used to calculate azimuth
        point1.setX(pBratLine->getX(0));
        point1.setY(pBratLine->getY(0));
        point2.setX(pBratLine->getX(nPoints-1));
        point2.setY(pBratLine->getY(nPoints-1));

        //length of BRAT segment
        length = pBratLine->get_Length();

        //slope of BRAT segment
        slope = pBratFeat->GetFieldAsDouble(slopeField);

        //number of dams for a BRAT segment randomly selected from complex size distribution
        nDamCount = pBratFeat->GetFieldAsInteger("totdams");
        nCompCount = pBratFeat->GetFieldAsInteger("totcomp");
        //qDebug()<<nDams<<damDens/1000.0<<nDamCount<<m_modCap<<length;

        if (nDamCount > 0)
        {
            spacing = length / (nDamCount*1.0);
        }
        else
        {
            spacing = 0.0;
        }
        //calculate azimuth (from start to end) of BRAT segment
        azimuthStart = Geometry::calcAzimuth(point2.getX(), point2.getY(), point1.getX(), point1.getY());
        //find elevation on stream network closest to dam point (on line perpendicular to BRAT segment)
        end_elev = raster_dem.sampleAlongLine_RasterVal(m_demPath, m_facPath, pBratLine->getX(0), pBratLine->getY(0), azimuthStart, sampleDist, endx, endy);

        //create a point for each dam to be modeled on BRAT segment
        int damReamin = nDamCount;
        int compRemain = nCompCount;

        for (int j=0; j<nDamCount; j++)
        {
            //location of modeled dam
            OGRPoint damPoint;
            //location on BRAT segment
            double pointDist = length - (spacing * (j * 1.0));
            //Determine if dam is primary or secondary and get dam heights from distribution
            double rnum = ((double) rand() / (RAND_MAX));
            double dht;
            //Heigth distribution for secondary dams
            Statistics normDist(Random::randomSeries(1000, RDT_norm, 0.92, 0.17), RDT_norm);
            dht = Random::random_normal(0.92, 0.17);

            //Different dam type determination methods for different sized dam complexes
            if (nDamCount > 2)
            {
                if (rnum <= (compRemain*1.0)/(damReamin*1.0))
                {
                    //Height distribution for primary dams
                    normDist.setSample(Random::randomSeries(1000, RDT_norm, 1.14, 0.19));
                    dht = Random::random_normal(1.14, 0.19);
                    nPrimary++;
                    damType = "primary";
                    compRemain--;
                }
                else
                {
                    nSecondary++;
                    damType = "secondary";
                }
            }
            //For small complexes
            else
            {
                if (rnum <= 0.15)
                {
                    normDist.setSample(Random::randomSeries(1000, RDT_norm, 1.14, 0.19));
                    dht = Random::random_normal(1.14, 0.19);
                    nPrimary++;
                    damType = "primary";
                }
                else
                {
                    nSecondary++;
                    damType = "secondary";
                }
            }
            //set point location on BRAT segment
            pBratLine->Value(pointDist, &damPoint);
            //location of dam point
            double x = damPoint.getX();
            double y = damPoint.getY();
            //move dam point to lowest cross sectional elevation
            //elev = raster_dem.sampleAlongLine_LowVal(m_demPath, damPoint.getX(), damPoint.getY(), azimuthStart, sampleDist, x, y);
            //move dam point to stream network raster
            elev = raster_dem.sampleAlongLine_RasterVal(m_demPath, m_facPath, damPoint.getX(), damPoint.getY(), azimuthStart, sampleDist, x, y);
            //set location of point
            damPoint.setX(x);
            damPoint.setY(y);
            //set field values and dam heights
            setFieldValues(pDamFeat, i, elev, slope, Geometry::calcAzimuth(damPoint.getX(), damPoint.getY(), endx, endy), x, y);
            double loHt, midHt, hiHt, maxHt;
            loHt = normDist.getQuantile(0.025), midHt = normDist.getQuantile(0.5), hiHt = normDist.getQuantile(0.975), maxHt = VectorOps::max(normDist.getData());
            setDamHeights(pDamFeat, loHt*loHt, midHt*midHt, hiHt*hiHt, maxHt*maxHt);
            pDamFeat->SetField("damType", damType);

            pDamFeat->SetGeometry(&damPoint);
            //qDebug()<<"Field Values "<<i<<elev<<slope<<Geometry::calcAzimuth(damPoint.getX(), damPoint.getY(), endx, endy)<<x<<y;
            //make sure point is located inside model domain (if the elevation value is NoData, outside domain)
            if (elev > 0.0 && end_elev > 0.0)
            {
                if (pDamsLyr->CreateFeature(pDamFeat) != OGRERR_NONE)
                {
                    nErrCount++;
                    qDebug()<<"error writing dam point "<<nErrCount;
                }
            }
            damReamin--;
        }
        i++;
    }
    qDebug()<<"Dams Modeled"<<nTotalDams<<"Max Capacity for Scenario"<<modCap<<"Maximum Capacity"<<maxCap<<"Percent Capacity"<<m_modCap;
    qDebug()<<"primary"<<nPrimary<<"secondary"<<nSecondary;
    OGRFeature::DestroyFeature(pDamFeat);
    OGRFeature::DestroyFeature(pBratFeat);
}

//Create dam points to model from known dam locations. Modeled dam points moved to flow accumulation raster
void DamPoints::createDamPoints_Copy(OGRLayer *pBratLyr, OGRLayer *pDamsLyr, OGRLayer *pExLyr)
{
    qDebug()<<"starting dam points from points";
    const char *slopeField = "iGeo_Slope";
    double sampleDist = 100.0;
    OGRFeature *pBratFeat, *pOldFeat;
    qDebug()<<"features declared";
    OGRFeature *pDamFeat = OGRFeature::CreateFeature(pDamsLyr->GetLayerDefn());
    qDebug()<<"feature created";
    Raster raster_dem;
    int nPrimary = 0, nSecondary = 0;
    qDebug()<<"getting features";
    int nFeatures = pExLyr->GetFeatureCount();
    qDebug()<<"starting loop";
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
        //Heigth distribution for secondary dams
        Statistics normDist(Random::randomSeries(1000, RDT_norm, 0.92, 0.17), RDT_norm);
        dht = Random::random_normal(0.92, 0.17);
        //Primary dam
        if (rnum <= 0.15)
        {
            //Height distribution for primary dams
            normDist.setSample(Random::randomSeries(1000, RDT_norm, 1.14, 0.19));
            dht = Random::random_normal(1.14, 0.19);
            nPrimary++;
        }
        else
        {
            nSecondary++;
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
        else
        {
            qDebug()<<"elevation error creating point "<<elev;
        }
    }
    qDebug()<<"primary"<<nPrimary<<"secondary"<<nSecondary;
    OGRFeature::DestroyFeature(pDamFeat);
    OGRFeature::DestroyFeature(pBratFeat);
    OGRFeature::DestroyFeature(pOldFeat);
}

//Create dam points to model at same location as input dam points. If dam extents created wtih polygons output will have errors
void DamPoints::createDamPoints_CopyLoc(OGRLayer *pBratLyr, OGRLayer *pDamsLyr, OGRLayer *pExLyr)
{
    qDebug()<<"starting dam points from points copy location";
    const char *slopeField = "iGeo_Slope";
    //double sampleDist = 50.0;
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
        //double az = Geometry::calcAzimuth(pBratLine->getX(pBratLine->getNumPoints()-1), pBratLine->getY(pBratLine->getNumPoints()-1), pBratLine->getX(0), pBratLine->getY(0));
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

        //double x = pOldDam->getX();
        //double y = pOldDam->getY();
        //double elev = raster_dem.sampleAlongLine_LowVal(m_demPath, pOldDam->getX(), pOldDam->getY(), az, sampleDist, x, y);
        //double elev = raster_dem.sampleAlongLine_RasterVal(m_demPath, m_facPath, pOldDam->getX(), pOldDam->getY(), az, sampleDist, x, y);
        //pOldDam->setX(x);
        //pOldDam->setY(y);
        double elev = raster_dem.valueAtPoint(m_demPath, pOldDam->getX(), pOldDam->getY());
        setFieldValues(pDamFeat, i, elev, slope, Geometry::calcAzimuth(pOldDam->getX(), pOldDam->getY(), pBratLine->getX(0), pBratLine->getY(0)), pOldDam->getX(), pOldDam->getY());
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

//Create dam points to model from dams with measured heights. Modeled dam points moved to flow accumulation raster
void DamPoints::createDamPoints_Heights(OGRLayer *pBratLyr, OGRLayer *pDamsLyr, OGRLayer *pExLyr)
{
    qDebug()<<"starting dam points height";
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
        double damHeight = pOldFeat->GetFieldAsDouble("DamHt_m"); //Field name for Temple Fork data is "Dam_Height", for all others "DamHt_m", for volume comparison "d_ht_m"
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

//Create dam points to model with measure heights at same location as input dam points. If dam extents created wtih polygons output will have errors
void DamPoints::createDamPoints_HeightsLoc(OGRLayer *pBratLyr, OGRLayer *pDamsLyr, OGRLayer *pExLyr)
{
    qDebug()<<"Dam points with heights, preserve original location";
    const char *slopeField = "iGeo_Slope";
    //double sampleDist = 50.0;
    OGRFeature *pBratFeat, *pOldFeat;
    OGRFeature *pDamFeat = OGRFeature::CreateFeature(pDamsLyr->GetLayerDefn());
    Raster raster_dem;
    int nFeatures = pExLyr->GetFeatureCount();

    for (int i=0; i<nFeatures; i++)
    {
        pOldFeat = pExLyr->GetFeature(i);
        int nBratFID = pOldFeat->GetFieldAsInteger("ID");
        double damHeight = pOldFeat->GetFieldAsDouble("DamHt_m");
        pBratFeat = pBratLyr->GetFeature(nBratFID);
        double slope = pBratFeat->GetFieldAsDouble(slopeField);
        OGRGeometry *pGeom = pBratFeat->GetGeometryRef();
        OGRLineString *pBratLine = (OGRLineString*) pGeom;
        pGeom = pOldFeat->GetGeometryRef();
        OGRPoint *pOldDam = (OGRPoint*) pGeom;
        //double az = Geometry::calcAzimuth(pBratLine->getX(pBratLine->getNumPoints()-1), pBratLine->getY(pBratLine->getNumPoints()-1), pBratLine->getX(0), pBratLine->getY(0));
        //double x = pOldDam->getX();
        //double y = pOldDam->getY();
        //double elev = raster_dem.sampleAlongLine_LowVal(m_demPath, pOldDam->getX(), pOldDam->getY(), az, sampleDist, x, y);
        //double elev = raster_dem.sampleAlongLine_RasterVal(m_demPath, m_facPath, pOldDam->getX(), pOldDam->getY(), az, sampleDist, x, y);
        //pOldDam->setX(x);
        //pOldDam->setY(y);
        double elev = raster_dem.valueAtPoint(m_demPath, pOldDam->getX(), pOldDam->getY());
        setFieldValues(pDamFeat, i, elev, slope, Geometry::calcAzimuth(pOldDam->getX(), pOldDam->getY(), pBratLine->getX(0), pBratLine->getY(0)), pOldDam->getX(), pOldDam->getY());
        setDamHeights(pDamFeat, damHeight, damHeight, damHeight, damHeight);
        //qDebug()<<i<< elev<< slope<< Geometry::calcAzimuth(pOldDam->getX(), pOldDam->getY(), pBratLine->getX(0), pBratLine->getY(0))<< pOldDam->getX()<< pOldDam->getY();

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
    field.SetName("slope");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("ht_max");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("damType");
    field.SetType(OFTString);
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
    field.SetName("ht_lo_mod");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("ht_mid_mod");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("ht_hi_mod");
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
    field.SetName("vol_lo_lp");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_mid_lp");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_hi_lp");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_lo_mp");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_mid_mp");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_hi_mp");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_lo_up");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_mid_up");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("vol_hi_up");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("diff_lo");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("diff_mid");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("diff_hi");
    field.SetType(OFTReal);
    pLayer->CreateField(&field);
    field.SetName("type");
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
    double slp = pFeat->GetFieldAsDouble("slope");
    pFeat->SetField("ht_lo", low);
    pFeat->SetField("ht_mid", mid);
    pFeat->SetField("ht_hi", high);
    pFeat->SetField("ht_lo_mod", low);
    pFeat->SetField("ht_mid_mod", mid);
    pFeat->SetField("ht_hi_mod", high);
    pFeat->SetField("ht_max", max);

    if (QString::fromUtf8(m_statPath) != "")
    {
        Raster raster;
        int x_lo = low/0.05 - 2, x_mid = mid/0.05 - 2, x_hi = high/0.05 - 2, y = slp/0.005;
        qDebug()<<high<<x_lo<<x_mid<<x_hi<<y;
//        if (x_lo == 0)
//        {
//            x_lo = 1;
//        }
//        if (x_mid == 0)
//        {
//            x_mid = 1;
//        }
//        if (x_hi == 0)
//        {
//            x_hi = 1;
//        }
//        if (y == 0)
//        {
//            y = 1;
//        }
        pFeat->SetField("vol_lo_lp", raster.value(m_statPath, y, x_lo, 2));
        pFeat->SetField("vol_mid_lp", raster.value(m_statPath, y, x_mid, 2));
        pFeat->SetField("vol_hi_lp", raster.value(m_statPath, y, x_hi, 2));
        pFeat->SetField("vol_lo_mp", raster.value(m_statPath, y, x_lo, 1));
        pFeat->SetField("vol_mid_mp", raster.value(m_statPath, y, x_mid, 1));
        pFeat->SetField("vol_hi_mp", raster.value(m_statPath, y, x_hi, 1));
        pFeat->SetField("vol_lo_up", raster.value(m_statPath, y, x_lo, 3));
        pFeat->SetField("vol_mid_up", raster.value(m_statPath, y, x_mid, 3));
        pFeat->SetField("vol_hi_up", raster.value(m_statPath, y, x_hi, 3));
    }
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

void DamPoints::sortByCapacity(OGRLayer *pBratLyr)
{
    qDebug()<<"Sorting BRAT by capacity";
    m_qvBratFID.clear();
    m_qvCapacityRank.clear();
    //qDebug()<<"Cleared";
    //BRAT line segment
    OGRFeature *pBratFeat;

    int nFeatures = pBratLyr->GetFeatureCount();
    double capacity;
    int nFid;

    //loop through all features in BRAT layer
    for (int i=0; i<nFeatures; i++)
    {
        //get BRAT feature
        pBratFeat = pBratLyr->GetFeature(i);
        capacity = pBratFeat->GetFieldAsDouble("oCC_EX");
        nFid = pBratFeat->GetFID();
        OGRGeometry *pGeom = pBratFeat->GetGeometryRef();
        //geometry of BRAT feature
        OGRLineString *pBratLine = (OGRLineString*) pGeom;

        //length of BRAT segment
        double length = pBratLine->get_Length();
        double maxDams = ceil(length * (capacity/1000.0));

        if (i > 0)
        {
            if (capacity >= m_qvCapacityRank[0])
            {
                //qDebug()<<"prepending";
                m_qvCapacityRank.prepend(capacity);
                m_qvBratFID.prepend(nFid);
                m_qvMaxDams.prepend(maxDams);
                //qDebug()<<"prepending done";
            }
            else if (capacity <= m_qvCapacityRank.last())
            {
                //qDebug()<<"appending";
                m_qvCapacityRank.append(capacity);
                m_qvBratFID.append(nFid);
                m_qvMaxDams.append(maxDams);
                //qDebug()<<"appending done";
            }
            else
            {
                bool stop = false;
                int j=0;
                while (j<m_qvCapacityRank.length() && !stop)
                {
                    if (capacity >= m_qvCapacityRank[j])
                    {
                        //qDebug()<<"inserting"<<stop<<nFeatures<<i;
                        m_qvCapacityRank.insert(j, capacity);
                        m_qvBratFID.insert(j, nFid);
                        m_qvMaxDams.insert(j, maxDams);
                        //qDebug()<<"inserting done";
                        stop = true;
                    }
                    j++;
                }
            }
        }
        else
        {
            qDebug()<<"adding first obs";
            m_qvBratFID.append(nFid);
            m_qvCapacityRank.append(capacity);
            m_qvMaxDams.append(maxDams);
            qDebug()<<"first added";
        }
    }
    OGRFeature::DestroyFeature(pBratFeat);
    //qDebug()<<"These should be the same"<<m_qvCapacityRank.length()<<m_qvBratFID.length()<<m_qvMaxDams.length();

//    qDebug()<<"BRAT segment RANK TESTING";
//    for (int i=0; i< m_qvCapacityRank.length(); i++)
//    {
//        qDebug()<<m_qvBratFID[i]<<m_qvCapacityRank[i];
//    }
}

bool DamPoints::setPondAttributes(OGRFeature *pFeat, double lowarea, double midarea, double hiarea, double lowvol, double midvol, double hivol)
{
    bool adjusted = true;
    QVector<QString> qsLevels;
    QVector<double> qvArea, qvVol;
    qvArea.append(lowarea), qvArea.append(midarea), qvArea.append(hiarea);
    qvVol.append(lowvol), qvVol.append(midvol), qvVol.append(hivol);
    qsLevels.append("lo"), qsLevels.append("mid"), qsLevels.append("hi");
    double ht, htMod, pred, lwr, upr, vol, diff;
    // amount by which to increment dam height
    double htInc = 0.1;

    for (int i=0; i<qsLevels.length(); i++)
    {
        ht = pFeat->GetFieldAsDouble(QString("ht_"+qsLevels[i]).toStdString().c_str());
        htMod = pFeat->GetFieldAsDouble(QString("ht_"+qsLevels[i]+"_mod").toStdString().c_str());
        pred = pFeat->GetFieldAsDouble(QString("vol_"+qsLevels[i]+"_mp").toStdString().c_str());
        lwr = pFeat->GetFieldAsDouble(QString("vol_"+qsLevels[i]+"_lp").toStdString().c_str());
        upr = pFeat->GetFieldAsDouble(QString("vol_"+qsLevels[i]+"_up").toStdString().c_str());
        vol = pFeat->GetFieldAsDouble(QString("vol_"+qsLevels[i]).toStdString().c_str());
        diff = pFeat->GetFieldAsDouble(QString("diff_"+qsLevels[i]).toStdString().c_str());

        if (pred <= 0.0 && lwr <= 0.0 && upr <= 0.0)
        {
            pFeat->SetField("type", 1.0);
            pFeat->SetField(QString("area_"+qsLevels[i]).toStdString().c_str(), qvArea[i]);
            pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), qvVol[i]);
        }
        else if (vol == 0.0)
        {
            if ((qvVol[i] >= lwr && qvVol[i] <= upr) || vol == -3.0)
            {
                pFeat->SetField("type", 2.0);
                pFeat->SetField(QString("area_"+qsLevels[i]).toStdString().c_str(), qvArea[i]);
                pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), qvVol[i]+0.1);
            }
            else if (qvVol[i] > upr)
            {
                if (htMod > htInc)
                {
                    pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), -1.0);
                    pFeat->SetField(QString("ht_"+qsLevels[i]+"_mod").toStdString().c_str(), htMod - htInc);
                    pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
                    adjusted = false;
                }
                else
                {
                    pFeat->SetField("type", 3.0);
                    pFeat->SetField(QString("area_"+qsLevels[i]).toStdString().c_str(), qvArea[i]);
                    pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), qvVol[i]+0.1);
                    pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
                }
            }
            else if (qvVol[i] < lwr)
            {
                if (htMod < 4.0)
                {
                    pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), -2.0);
                    pFeat->SetField(QString("ht_"+qsLevels[i]+"_mod").toStdString().c_str(), htMod + htInc);
                    pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
                    adjusted = false;
                }
                else
                {
                    //qDebug()<<"IN ZERO: MODELED HEIGHT ABOVE THRESHOLD!!!!"<<htMod<<vol;
                    pFeat->SetField("type", 4.0);
                    pFeat->SetField(QString("area_"+qsLevels[i]).toStdString().c_str(), qvArea[i]);
                    pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), qvVol[i]+0.1);
                    pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
                }
            }
        }
        else if (vol < 0.0)
        {
            if (fabs(qvVol[i] - pred) < 5.0)
            {
                pFeat->SetField("type", 5.0);
                pFeat->SetField(QString("area_"+qsLevels[i]).toStdString().c_str(), qvArea[i]+0.01);
                pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), qvVol[i]+0.01);
                pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
            }
            else if (vol == -1.0)
            {
                if ((qvVol[i] - pred) > 0.0)
                {
                    if (htMod > htInc)
                    {
                        //qDebug()<<"NEGATIVE ONE: HEIGHT CHANGED"<<htMod - 0.05;
                        pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), -1.0);
                        pFeat->SetField(QString("ht_"+qsLevels[i]+"_mod").toStdString().c_str(), htMod - htInc);
                        pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
                        adjusted = false;
                    }
                    else
                    {
                        pFeat->SetField("type", 6.0);
                        pFeat->SetField(QString("area_"+qsLevels[i]).toStdString().c_str(), qvArea[i]);
                        pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), qvVol[i]);
                        pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
                    }
                }
                else
                {
                    if(fabs(diff) > fabs(pred - qvVol[i]))
                    {
                        pFeat->SetField("type", 7.0);
                        pFeat->SetField(QString("area_"+qsLevels[i]).toStdString().c_str(), qvArea[i]);
                        pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), qvVol[i]);
                        pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
                    }
                    else
                    {
                        pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), -3.0);
                        pFeat->SetField(QString("ht_"+qsLevels[i]+"_mod").toStdString().c_str(), htMod + htInc);
                        pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
                        adjusted = false;
                    }
                }
            }
            else if (vol == -2.0)
            {
                if ((qvVol[i] - pred) < 0.0)
                {
                    if (htMod < 4.0)
                    {
                        pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), -2.0);
                        pFeat->SetField(QString("ht_"+qsLevels[i]+"_mod").toStdString().c_str(), htMod + htInc);
                        pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
                        adjusted = false;
                    }
                    else
                    {
                        pFeat->SetField("type", 8.0);
                        pFeat->SetField(QString("area_"+qsLevels[i]).toStdString().c_str(), qvArea[i]);
                        pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), qvVol[i]);
                        pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
                    }
                }
                else
                {
                    if(fabs(diff) >= fabs(pred - qvVol[i]))
                    {
                        pFeat->SetField("type", 9.0);
                        pFeat->SetField(QString("area_"+qsLevels[i]).toStdString().c_str(), qvArea[i]);
                        pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), qvVol[i]);
                        pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
                    }
                    else
                    {
                        pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), -3.0);
                        pFeat->SetField(QString("ht_"+qsLevels[i]+"_mod").toStdString().c_str(), htMod - htInc);
                        pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
                        adjusted = false;
                    }
                }
            }
            else if (vol == -3.0)
            {
                pFeat->SetField("type", 10.0);
                pFeat->SetField(QString("area_"+qsLevels[i]).toStdString().c_str(), qvArea[i]);
                pFeat->SetField(QString("vol_"+qsLevels[i]).toStdString().c_str(), qvVol[i]);
                pFeat->SetField(QString("diff_"+qsLevels[i]).toStdString().c_str(), qvVol[i] - pred);
            }
        }
    }

    return adjusted;
}

