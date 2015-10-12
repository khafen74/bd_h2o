#include <QCoreApplication>
#include <QDebug>
#include "gdal.h"
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "ogr_core.h"
#include "ogr_api.h"

const double PI = 3.14159265;
const double ANGLE_OFFSET[5] = {-90.0, -45.0, 0.0, 45.0, 90.0};

double addDegrees(double base, double addValue);
double angleBetweenLines(double x1, double y1, double x2, double y2, double x3, double y3);
double calcAzimuth(double startX, double startY, double endX, double endY);
int calcCoords(double startX, double startY, double azimuth, double distance, double &newX, double &newY);
int createDamPoints(const char *demPath, const char *inputFeaturePath, int nFeatureType, const char *outputFeaturePath = 0);
int createRasterFromPoint(const char *rasterPath, const char *pointPath, int rows, int cols, double transform[]);
int createSearchPolygons(const char *outputFeaturePath);
int getRasterCol(double transform[6], double xCoord);
int getRasterRow(double transform[6], double yCoord);
double getRasterValueAtPoint(const char *rasterPath, double xCoord, double yCoord);
int pointsInPolygon(const char *pointsPath, const char *polygonPath);
double sampleRasterAlongLine_LowVal(const char * rasterPath, double startX, double startY, double azimuth, double distance, double &x, double &y);

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    GDALAllRegister();
    OGRRegisterAll();

    const char *shpIn = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/01_shpIn";
    const char *shpOut = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut";
    const char *demIn = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/templefk_10m_ws.tif";
    const char *depOut = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/ponddepth_10m.tif";

    GDALDataset *pDem = (GDALDataset*) GDALOpen(demIn, GA_ReadOnly);
    int rows = pDem->GetRasterYSize();
    int cols = pDem->GetRasterXSize();
    double transform[6];
    pDem->GetGeoTransform(transform);
    GDALClose(pDem);

    //qDebug()<<angleBetweenLines(0.75,1,1,0.75,0,0);
    pointsInPolygon(shpIn, shpOut);
    //createDamPoints(demIn, shpIn, 1, shpOut);
    //createSearchPolygons(shpOut);
    //createRasterFromPoint(depOut, shpOut, rows, cols, transform);

    qDebug()<<"done";

    return a.exec();
}

double addDegrees(double base, double addValue)
{
    double value, remainder;

    value = base + addValue;

    if (value > 360.0)
    {
        remainder = value - 360.0;
        value = remainder;
    }
    else if (value < 0.0)
    {
        remainder = fabs(value);
        value = 360.0 - remainder;
    }

    return value;
}

double angleBetweenLines(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double dx1, dy1, dx2, dy2, d, l2, m1, m2, angle;

    dx1 = x3-x1;
    dy1 = y3-y1;
    dx2 = x3-x2;
    dy2 = y3-y2;

    m1 = dy1/dx1;
    m2 = dy2/dx2;
    //qDebug()<<"slopes"<<m1<<m2;

    //qDebug()<<"result"<<(fabs((m1-m2)/(1+m1*m2)));
    angle = tan(fabs((m1-m2)/(1+m1*m2)));
    if ((angle*180/PI) > 180.0)
    {
        qDebug()<<"error angle too big"<<angle*180/PI;
    }
    //qDebug()<<"angle rad "<<angle<<" angle deg "<< angle * 180.0/PI;

    return angle * 180.0/PI;
}

double calcAzimuth(double startX, double startY, double endX, double endY)
{
    double theta, azimuth;

    if (startX > endX && startY > endY)
    {
        theta = atan2((startY-endY),(startX-endX)) * 180.0 / PI;
        azimuth = 180.0 + theta;
    }
    else if (startX < endX && startY > endY)
    {
        theta = atan2((startY-endY),(endX-startX)) * 180.0 / PI;
        azimuth = 360.0 - theta;
    }
    else if (startX < endX && startY < endY)
    {
        theta = atan2((endY-startY),(endX-startX)) * 180.0 / PI;
        azimuth = theta;
    }
    else if (startX > endX && startY < endY)
    {
        theta = atan2((startX-endX),(endY-startY)) * 180.0 / PI;
        azimuth = 90.0 + theta;
    }
    else
    {
        qDebug()<<"AZIMUTH ERROR: may be a straight line";
    }

    return azimuth;
}

int calcCoords(double startX, double startY, double azimuth, double distance, double &newX, double &newY)
{
    double deltaX, deltaY, theta;

    if (azimuth > 0.0 && azimuth < 90.0)
    {
        theta = azimuth;
        deltaY = sin(theta*(PI/180.0))*distance;
        deltaX = cos(theta*(PI/180.0))*distance;
    }
    else if (azimuth > 90.0 && azimuth < 180.0)
    {
        theta = azimuth - 90.0;
        deltaX = sin(theta*(PI/180.0))*distance*(-1.0);
        deltaY = cos(theta*(PI/180.0))*distance;
    }
    else if (azimuth > 180.0 && azimuth < 270.0)
    {
        theta = azimuth - 180.0;
        deltaY = sin(theta*(PI/180.0))*distance*(-1.0);
        deltaX = cos(theta*(PI/180.0))*distance*(-1.0);
    }
    else if (azimuth > 270.0 && azimuth < 360.0)
    {
        theta = 360.0 - azimuth;
        deltaY = sin(theta*(PI/180.0))*distance*(-1.0);
        deltaX = cos(theta*(PI/180.0))*distance;
    }
    else
    {

    }

    newX = startX + deltaX;
    newY = startY + deltaY;

    return 0;
}

int createDamPoints(const char *demPath, const char *inputFeaturePath, int nFeatureType, const char *outputFeaturePath)
{
    //GDALDataset *pDem;
    OGRDataSource *pInDs, *pOutDs;
    OGRLayer *pDamsIn, *pDamsOut, *pBratIn;

    GDALDriver *pDriverTiff;
    pDriverTiff = GetGDALDriverManager()->GetDriverByName("GTiff");
    qDebug()<<"gdal driver created";

    OGRSFDriver *pDriverShp;
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    pDriverShp = registrar->GetDriverByName("ESRI Shapefile");
    qDebug()<<"ogr driver created";

    pInDs = pDriverShp->CreateDataSource(inputFeaturePath, NULL);
    pOutDs = pDriverShp->CreateDataSource(outputFeaturePath, NULL);

    pDamsIn = pInDs->GetLayerByName("Dams_BRAT_join5_UTM12N");
    pBratIn = pInDs->GetLayerByName("BRAT_TempleFk_WS");

    pDamsOut = pOutDs->CreateLayer("ModeledDamPoints", pBratIn->GetSpatialRef(), wkbPoint, NULL);
    qDebug()<<"creating fields";

    OGRFieldDefn field("endx", OFTReal);
    pDamsOut->CreateField(&field);
    field.SetName("endy");
    field.SetType(OFTReal);
    pDamsOut->CreateField(&field);
    field.SetName("az_us");
    field.SetType(OFTReal);
    pDamsOut->CreateField(&field);
    field.SetName("g_elev");
    field.SetType(OFTReal);
    pDamsOut->CreateField(&field);
    field.SetName("d_elev");
    field.SetType(OFTReal);
    pDamsOut->CreateField(&field);
    field.SetName("slope");
    field.SetType(OFTReal);
    pDamsOut->CreateField(&field);
    qDebug()<<"done";

    double x, y, az, value, sampleDist, slope;
    double damHeight = 1.0;
    OGRFeature *pDamFeatureOld, *pDamFeatureNew, *pBratFeature;
    OGRGeometry *pGeom;
    OGRPoint *pDamOld, pDamNew;
    OGRLineString *pLineString;

    int nDams, nBratFID;
    nDams = pDamsIn->GetFeatureCount();
    pDamFeatureNew = OGRFeature::CreateFeature(pDamsOut->GetLayerDefn());
    qDebug()<<"starting loop";

    for (int i=0; i<nDams; i++)
    {
        pDamFeatureOld = pDamsIn->GetFeature(i);
        nBratFID = pDamFeatureOld->GetFieldAsInteger("ID");
        qDebug()<<nBratFID;
        pBratFeature = pBratIn->GetFeature(nBratFID);
        pDamFeatureOld = pDamsIn->GetFeature(i);

        pGeom = pBratFeature->GetGeometryRef();
        pLineString = (OGRLineString*) pGeom;
        pGeom = pDamFeatureOld->GetGeometryRef();
        pDamOld = (OGRPoint*) pGeom;
        az = calcAzimuth(pDamOld->getX(), pDamOld->getY(), pLineString->getX(0), pLineString->getY(0));

        //qDebug()<<az<<pDamOld->getX()<<pDamOld->getY()<<pLineString->getX(0)<<pLineString->getY(0);

        x = pDamOld->getX(), y = pDamOld->getY();

        slope = pBratFeature->GetFieldAsDouble("iGeo_Slope");
        sampleDist = damHeight/slope;
        value = sampleRasterAlongLine_LowVal(demPath, pDamOld->getX(), pDamOld->getY(), az, sampleDist, x, y);
        //qDebug()<<"setting fields";
        pDamNew.setX(x);
        pDamNew.setY(y);
        //qDebug()<<"xy set";
        pDamFeatureNew->SetField("g_elev", value);
        pDamFeatureNew->SetField("d_elev", value+damHeight);
        pDamFeatureNew->SetField("slope", slope);
        //qDebug()<<"elevs set";
        value = sampleRasterAlongLine_LowVal(demPath, pLineString->getX(0), pLineString->getY(0), az, sampleDist, x, y);
        //qDebug()<<"endpoint found";
        pDamFeatureNew->SetField("endx", x);
        pDamFeatureNew->SetField("endy", y);
        //qDebug()<<"enpoint added";
        pDamFeatureNew->SetField("az_us", calcAzimuth(pDamNew.getX(), pDamNew.getY(), x, y));
        //qDebug()<<"azimuth added";
        pDamFeatureNew->SetGeometry(&pDamNew);
        pDamsOut->CreateFeature(pDamFeatureNew);

        qDebug()<<x<<y<<value<<getRasterValueAtPoint(demPath, pLineString->getX(0), pLineString->getY(0));
    }
    OGRFeature::DestroyFeature(pDamFeatureNew);
    OGRFeature::DestroyFeature(pDamFeatureOld);
    OGRFeature::DestroyFeature(pBratFeature);

    OGRDataSource::DestroyDataSource(pInDs);
    OGRDataSource::DestroyDataSource(pOutDs);
}

int createRasterFromPoint(const char *rasterPath, const char *pointPath, int rows, int cols, double transform[6])
{
    GDALDataset *pRaster;
    GDALDriver *pDriver;
    OGRDataSource *pPoints;

    double noData = -9999;

    pDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    pRaster = pDriver->Create(rasterPath,cols, rows, 1, GDT_Float32, NULL);
    pRaster->SetGeoTransform(transform);
    pRaster->GetRasterBand(1)->Fill(noData);
    pRaster->GetRasterBand(1)->SetNoDataValue(noData);

    pPoints = OGRSFDriverRegistrar::Open(pointPath);
    OGRLayer *layer = pPoints->GetLayerByName("PondPts");
    int nFeatures = layer->GetFeatureCount();

    QString damElevName = "d_elev";
    int nDamEl, nDemEl, row, col;
    QString demElevName = "GRID_CODE";
    QString fieldName;
    double damElev, demElev, watDep;

    OGRFeatureDefn *featDfn = layer->GetLayerDefn();

    for (int i=0; i<featDfn->GetFieldCount(); i++)
    {
        fieldName = QString::fromUtf8(featDfn->GetFieldDefn(i)->GetNameRef());
        if (fieldName == damElevName)
        {
            nDamEl = i;
        }
        else if (fieldName == demElevName)
        {
            nDemEl = i;
        }
    }

    OGRFeature *feature;

    float *depth = (float*) CPLMalloc(sizeof(float)*1);

    layer->ResetReading();
    while((feature = layer->GetNextFeature()) != NULL)
    {
        damElev = feature->GetFieldAsDouble(nDamEl);
        demElev = feature->GetFieldAsDouble(nDemEl);
        watDep = damElev - demElev;

        OGRPoint *point = (OGRPoint*) feature->GetGeometryRef();

        if (watDep > 0.0)
        {
            col = getRasterCol(transform, point->getX());
            row = getRasterRow(transform, point->getY());
            *depth = watDep;

            pRaster->GetRasterBand(1)->RasterIO(GF_Write, col, row, 1, 1, depth, 1, 1, GDT_Float32, 0,0);
        }
    }

    OGRFeature::DestroyFeature(feature);
    OGRDataSource::DestroyDataSource(pPoints);

    CPLFree(depth);

    GDALClose(pRaster);
}

int createSearchPolygons(const char *outputFeaturePath)
{
    OGRSFDriver *pDriverShp;
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    pDriverShp = registrar->GetDriverByName("ESRI Shapefile");

    OGRDataSource *pOutDs = pDriverShp->CreateDataSource(outputFeaturePath, NULL);
    OGRLayer *pDamPoints, *pSearchPolygons;
    pDamPoints = pOutDs->GetLayerByName("ModeledDamPoints");
    pSearchPolygons = pOutDs->CreateLayer("DamSearchPolygons", pDamPoints->GetSpatialRef(), wkbPolygon, NULL);
    OGRFieldDefn field("d_elev", OFTReal);
    pSearchPolygons->CreateField(&field);

    double xCoords[5], yCoords[5];
    double azimuthCurrent, azimuthStart, distance, slope;
    int nDams = pDamPoints->GetFeatureCount();
    qDebug()<<nDams;

    for (int i=0; i<nDams; i++)
    {
        OGRFeature *pDamFeature;
        OGRPolygon polygon;
        OGRLinearRing ring;
        OGRPoint point, *damPoint;
        OGRFeature *pPolyFeature = OGRFeature::CreateFeature(pSearchPolygons->GetLayerDefn());

        pDamFeature = pDamPoints->GetFeature(i);
        pPolyFeature->SetField("d_elev", pDamFeature->GetFieldAsDouble("d_elev"));
        azimuthStart = pDamFeature->GetFieldAsDouble("az_us");
        slope = pDamFeature->GetFieldAsDouble("slope");
        distance = (pDamFeature->GetFieldAsDouble("d_elev")-pDamFeature->GetFieldAsDouble("g_elev"))/slope;
        qDebug()<<pDamFeature->GetFieldAsDouble("slope")<<distance;
        damPoint = (OGRPoint*) pDamFeature->GetGeometryRef();

        for (int j=0; j<5; j++)
        {
            azimuthCurrent = addDegrees(azimuthStart, ANGLE_OFFSET[j]);
            calcCoords(damPoint->getX(), damPoint->getY(), azimuthCurrent, distance, xCoords[j], yCoords[j]);

            point.setX(xCoords[j]);
            point.setY(yCoords[j]);

            ring.addPoint(&point);
        }

        ring.addPoint(xCoords[0], yCoords[0]);
        polygon.addRing(&ring);
        pPolyFeature->SetGeometry(&polygon);
        pSearchPolygons->CreateFeature(pPolyFeature);
        OGRFeature::DestroyFeature(pPolyFeature);
        OGRFeature::DestroyFeature(pDamFeature);
    }

    OGRDataSource::DestroyDataSource(pOutDs);
}

int getRasterCol(double transform[6], double xCoord)
{
    int col;
    double xOffset, xDiv;

    xOffset = xCoord - transform[0];
    xDiv = xOffset/transform[1];
    col = floor(xDiv);

    return col;
}

int getRasterRow(double transform[6], double yCoord)
{
    int row;

    double yOffset, yDiv;

    yOffset = transform[3] - yCoord;
    yDiv = yOffset/fabs(transform[5]);
    row = floor(yDiv);

    return row;
}

double getRasterValueAtPoint(const char *rasterPath, double xCoord, double yCoord)
{
    GDALDataset *pRaster;
    double transform[6];
    double value, xOffset, yOffset, xDiv, yDiv;
    int row, col;

    pRaster = (GDALDataset*) GDALOpen(rasterPath, GA_ReadOnly);

    pRaster->GetGeoTransform(transform);

    xOffset = xCoord - transform[0];
    yOffset = transform[3] - yCoord;

    xDiv = xOffset/transform[1];
    yDiv = yOffset/transform[1];

    row = floor(yDiv);
    col = floor(xDiv);

    float *rasVal = (float*) CPLMalloc(sizeof(float)*1);

    pRaster->GetRasterBand(1)->RasterIO(GF_Read, col, row, 1, 1, rasVal, 1, 1, GDT_Float32, 0, 0);

    value = *rasVal;

    GDALClose(pRaster);
    CPLFree(rasVal);

    return value;
}

int pointsInPolygon(const char *pointsPath, const char *polygonPath)
{
    qDebug()<<"starting points to poly";
    OGRDataSource *pPointDS, *pPolyDS;
    OGRSFDriver *pDriverShp;
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    pDriverShp = registrar->GetDriverByName("ESRI Shapefile");

    //qDebug()<<"loading data";
    pPointDS = pDriverShp->CreateDataSource(pointsPath);
    //qDebug()<<"points loaded";
    pPolyDS = pDriverShp->CreateDataSource(polygonPath);
    //qDebug()<<"polys loaded";
    OGRLayer *pPointsLayer = pPointDS->GetLayerByName("dempoints_10m");
    //qDebug()<<"dem points loaded"<<pPointsLayer->GetFeatureCount();
    OGRLayer *pPolyLayer = pPolyDS->GetLayerByName("DamSearchPolygons");
    //qDebug()<<"pond polys loaded";
    OGRLayer *pDamPointLayer = pPolyDS->CreateLayer("PondPtsTest", pPointsLayer->GetSpatialRef(), wkbPoint, NULL);
    //qDebug()<<"new points created";
    OGRFieldDefn field("endx", OFTReal);
    pDamPointLayer->CreateField(&field);

    qDebug()<<"starting loop";
    int nPolyCount = pPolyLayer->GetFeatureCount();
    double angleDeg;
    qDebug()<<"polygons"<<nPolyCount;
    for (int i=0; i<nPolyCount; i++)
    {
        OGRFeature *pPolyFeat = pPolyLayer->GetFeature(i);
        OGRPolygon *pPoly = (OGRPolygon*) pPolyFeat->GetGeometryRef();
        OGRLinearRing *pRing = pPoly->getExteriorRing();
        qDebug()<<pPointsLayer->GetFeatureCount()<<pRing->getNumPoints();

        for (int j=0; j<pPointsLayer->GetFeatureCount(); j++)
        {
            OGRFeature *pPointFeat = pPointsLayer->GetFeature(j);
            OGRPoint *pPoint = (OGRPoint*) pPointFeat->GetGeometryRef();
            angleDeg = 0.0;
            //qDebug()<<"calcing intersect";
            for (int k=0; k<pRing->getNumPoints()-1; k++)
            {
                //qDebug()<<pRing->getX(k)<< pRing->getY(k)<< pRing->getX(k+1)<< pRing->getY(k+1)<< pPoint->getX()<< pPoint->getY();
                angleDeg += angleBetweenLines(pRing->getX(k), pRing->getY(k), pRing->getX(k+1), pRing->getY(k+1), pPoint->getX(), pPoint->getY());
                //qDebug()<<angleDeg;
                //qDebug()<<angleDeg;

            }
            //qDebug()<<"final angle"<<angleDeg;
            if (angleDeg < 180)
            {
                //qDebug()<<angleDeg;
            }
            else
            {
                qDebug()<<"add point "<<angleDeg;
                OGRFeature *newFeature = OGRFeature::CreateFeature(pDamPointLayer->GetLayerDefn());
                qDebug()<<"feature created";
                OGRPoint newPoint;
                qDebug()<<"point created";
                newPoint.setX(pPoint->getX());
                qDebug()<<"xset";
                newPoint.setY(pPoint->getY());
                qDebug()<<"yset";
                newFeature->SetGeometry(&newPoint);
                qDebug()<<"geom set";
                pDamPointLayer->SetFeature(newFeature);
                qDebug()<<"feature set";
                OGRFeature::DestroyFeature(newFeature);
            }
            //system("pause");

            OGRFeature::DestroyFeature(pPointFeat);
        }

        OGRFeature::DestroyFeature(pPolyFeat);
        qDebug()<<"finished feature"<<i;
    }

    OGRDataSource::DestroyDataSource(pPointDS);
    OGRDataSource::DestroyDataSource(pPolyDS);

    return 0;
}

double sampleRasterAlongLine_LowVal(const char * rasterPath, double startX, double startY, double azimuth, double distance, double &x, double &y)
{
    //qDebug()<<"sampling";
    double az1, az2, interval;
    double transform[6];
    double newX, newY, rasValue, lowValue;
    int nSamples;
    az1 = addDegrees(azimuth, 90.0);
    az2 = addDegrees(azimuth, -90.0);

    //qDebug()<<"opening raster";
    GDALDataset *pRas;
    pRas = (GDALDataset*) GDALOpen(rasterPath, GA_ReadOnly);
    pRas->GetGeoTransform(transform);
    GDALClose(pRas);
    //qDebug()<<"raster closed";

    interval = transform[1];
    nSamples = ceil(distance/interval);

    //qDebug()<<"starting loop";
    for (int i=0; i<nSamples; i++)
    {
        calcCoords(startX, startY, az1, interval*(i+1), newX, newY);
        rasValue = getRasterValueAtPoint(rasterPath, newX, newY);
        if (i==0)
        {
            lowValue = rasValue;
            x = newX, y = newY;
        }
        else
        {
            if (rasValue < lowValue)
            {
                lowValue = rasValue;
                x = newX, y = newY;
            }
        }
        calcCoords(startX, startY, az2, interval*(i+1), newX, newY);
        rasValue = getRasterValueAtPoint(rasterPath, newX, newY);
        if (rasValue < lowValue)
        {
            lowValue = rasValue;
            x = newX, y = newY;
        }
    }
    //qDebug()<<"done sampling";

    return lowValue;
}
