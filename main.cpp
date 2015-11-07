#include <QCoreApplication>
#include <QtCore>
#include <QDebug>
#include "gdal.h"
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "ogr_core.h"
#include "ogr_api.h"
#include "raster.h"

const double ANGLE_OFFSET[5] = {-90.0, -45.0, 0.0, 45.0, 90.0};

int run();
int testRegions();
int testRandom();
double addDegrees(double base, double addValue);
double angleBetweenLines(double x1, double y1, double x2, double y2, double x3, double y3);
double calcAzimuth(double startX, double startY, double endX, double endY);
int calcCoords(double startX, double startY, double azimuth, double distance, double &newX, double &newY);
int cleanup(const char *cleanDir);
int cleanInundationRaster(const char *rasterPath);
double createDamPoints(const char *demPath, const char *inputFeaturePath, const char *outputLayerName, const char *layerName);
int createInundationRaster(const char *rasterPath, int rows, int cols, double transform[]);
int createRasterFromPoint(const char *rasterPath, const char *pointPath, int rows, int cols, double transform[]);
int createSearchPolygons(const char *outputFeaturePath);
double getDamHeight();
double getDamHeightLnorm();
int getRasterCol(double transform[6], double xCoord);
int getRasterRow(double transform[6], double yCoord);
double getRasterX(double transform[], int col);
double getRasterY(double transform[], int row);
double getRasterValueAtPoint(const char *rasterPath, double xCoord, double yCoord);
double getWetArea(const char *rasterPath);
int pointsInPolygon2(const char *rasterPath, const char *rasterOut, const char *polygonPath);
int pointsInPolygon(const char *pointsPath, const char *polygonPath, const char *pointLayerName);
int regionsToDepth(const char *regRas, const char *depRas, const char *demRas, int nRegions);
double sampleRasterAlongLine_LowVal(const char * rasterPath, double startX, double startY, double azimuth, double distance, double &x, double &y);
int summarizeInundationRaster(const char *rasterPath, const char *outputCsv, int nValues, QVector<int> &thresholds, QVector<double> &areas);
int updateInundationRaster(const char *updatePath, const char *inputPath);
double uniformRandom();

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QDateTime startTime = QDateTime::currentDateTime();

    run();
    //testRegions();
    //testRandom();

    QDateTime endTime = QDateTime::currentDateTime();

    qDebug()<<"done"<<startTime.secsTo(endTime)<<startTime.secsTo(endTime)/60.0;

    return a.exec();
}

int run()
{
    GDALAllRegister();
    OGRRegisterAll();

//    const char *shpIn = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/01_shpIn";
//    const char *shpOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut";

    const char *shpIn = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/01_shpIn";
    const char *shpOut = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/03_shpOut";

//    const char *demIn = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/templefk_10m_ws.tif";
//    const char *depOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/ponddepth_10m.tif";
//    const char *freqOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_10m.tif";
//    const char *csv = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_10m.csv";

    const char *demIn = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/templefk_10m_ws.tif";
    const char *depOut = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/ponddepth_10m.tif";
    const char *freqOut = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_10m.tif";
    const char *csv = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_10m.csv";

//    const char *demIn = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/fme450000.tif";
//    const char *depOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/ponddepth_1m.tif";
//    const char *freqOut = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_1m.tif";
//    const char *csv = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_1m.csv";


//    const char *binRas = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/bin.tif";
//    const char *regRas = "E:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/reg.tif";

    const char *binRas = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/bin.tif";
    const char *regRas = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/reg.tif";

    const char *damLayerName = "Dams_BRAT_join5_UTM12N";
    //const char *damLayerName = "Dams_ESRI";

    GDALDataset *pDem = (GDALDataset*) GDALOpen(demIn, GA_ReadOnly);
    int rows = pDem->GetRasterYSize();
    int cols = pDem->GetRasterXSize();
    double transform[6];
    QVector<int> qvWetThresh;
    QVector<double> qvWetArea;
    pDem->GetGeoTransform(transform);
    GDALClose(pDem);

    createInundationRaster(freqOut, rows, cols, transform);
    int nIterations = 1000;
    double damHeightSum = 0.0, areaSum = 0.0;

    srand(time(NULL));
    for (int i=0; i<nIterations; i++)
    {
        cleanup(shpOut);
        damHeightSum = createDamPoints(demIn, shpIn, shpOut, damLayerName);
        createSearchPolygons(shpOut);
        pointsInPolygon2(demIn, depOut,shpOut);
        updateInundationRaster(freqOut, depOut);
        areaSum = getWetArea(depOut);
        qDebug()<<"Finished"<<i+1<<" of "<<nIterations<< ": Area = "<<areaSum<<" : Mean Dam Height = "<<damHeightSum;
    }

    cleanInundationRaster(freqOut);
    summarizeInundationRaster(freqOut, csv, nIterations, qvWetThresh, qvWetArea);

    Raster raster;
    raster.greaterThan(freqOut, binRas, 650.0);
    int regions = raster.regions(binRas, regRas);
    regionsToDepth(regRas, depOut, demIn, regions);

    return 0;
}

int testRegions()
{
    const char *freqOut = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/freqwet_10m.tif";
    const char *binRas = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/bin.tif";
    const char *regRas = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/reg.tif";
    const char *demIn = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/02_rasIn/templefk_10m_ws.tif";
    const char *depOut = "C:/etal/Projects/NonLoc/BeaverModeling/02_Data/z_TestRuns/04_rasOut/ponddepth_10m.tif";

    Raster raster;
    raster.greaterThan(freqOut, binRas, 65.0);
    int count = raster.regions(binRas, regRas);
    qDebug()<<"Regions"<<count;
    regionsToDepth(regRas, depOut, demIn, count);

    return 0;
}

int testRandom()
{
    for (int i=0; i<1000; i++)
    {
        qDebug()<<getDamHeight();
        //qDebug()<<uniformRandom();
    }
    return 0;
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
    double angle;

    angle = atan2(y1-y3, x1-x3) - atan2(y2-y3, x2-x3);

    while (angle < -PI)
    {
        angle += 2*PI;
    }
    while (angle > PI)
    {
        angle -= 2*PI;
    }

    return angle;
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

int cleanup(const char *cleanDir)
{
    QString path = QString::fromUtf8(cleanDir);
    QDir dir(path);
    dir.setNameFilters(QStringList() << "*.*");

    foreach (QString dirFile, dir.entryList())
    {
        dir.remove(dirFile);
    }

    return 0;
}

int cleanInundationRaster(const char *rasterPath)
{
    GDALDataset *pRaster;
    pRaster = (GDALDataset*) GDALOpen(rasterPath, GA_Update);

    int nCols = pRaster->GetRasterXSize();
    int nRows = pRaster->GetRasterYSize();

    float *rasVal = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pRaster->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, rasVal, nCols, 1, GDT_Float32, 0, 0);

        for (int j=0; j<nCols; j++)
        {
            if (rasVal[j] <= 0.0)
            {
                rasVal[j] = -9999;
            }
        }

        pRaster->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, rasVal, nCols, 1, GDT_Float32, 0, 0);
    }

    pRaster->GetRasterBand(1)->SetNoDataValue(-9999);

    CPLFree(rasVal);

    GDALClose(pRaster);

    return 0;
}

double createDamPoints(const char *demPath, const char *inputFeaturePath, const char *outputFeaturePath, const char *layerName)
{
    OGRDataSource *pInDs, *pOutDs;
    OGRLayer *pDamsIn, *pDamsOut, *pBratIn;

    GDALDriver *pDriverTiff;
    pDriverTiff = GetGDALDriverManager()->GetDriverByName("GTiff");

    OGRSFDriver *pDriverShp;
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    pDriverShp = registrar->GetDriverByName("ESRI Shapefile");

    pInDs = pDriverShp->CreateDataSource(inputFeaturePath, NULL);
    pOutDs = pDriverShp->CreateDataSource(outputFeaturePath, NULL);

    pDamsIn = pInDs->GetLayerByName(layerName);
    pBratIn = pInDs->GetLayerByName("BRAT_TempleFk_WS");

    pDamsOut = pOutDs->CreateLayer("ModeledDamPoints", pBratIn->GetSpatialRef(), wkbPoint, NULL);

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

    double x, y, az, value, sampleDist, slope;
    double damHeight;
    OGRFeature *pDamFeatureOld, *pDamFeatureNew, *pBratFeature;
    OGRGeometry *pGeom;
    OGRPoint *pDamOld, pDamNew;
    OGRLineString *pLineString;

    int nDams, nBratFID;
    double damSum = 0.0;
    nDams = pDamsIn->GetFeatureCount();
    pDamFeatureNew = OGRFeature::CreateFeature(pDamsOut->GetLayerDefn());

    for (int i=0; i<nDams; i++)
    {
        damHeight = getDamHeight();
        damSum += damHeight;
        if (damHeight < 0.3 || damHeight > 2.5)
        {
            //qDebug()<<"Dam Height Error: "<<damHeight;
        }
        else
        {
            //qDebug()<<damHeight;
        }
        pDamFeatureOld = pDamsIn->GetFeature(i);
        nBratFID = pDamFeatureOld->GetFieldAsInteger("ID");
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
    }
    OGRFeature::DestroyFeature(pDamFeatureNew);
    OGRFeature::DestroyFeature(pDamFeatureOld);
    OGRFeature::DestroyFeature(pBratFeature);

    OGRDataSource::DestroyDataSource(pInDs);
    OGRDataSource::DestroyDataSource(pOutDs);

    return damSum/(nDams*1.0);
}

int createInundationRaster(const char *rasterPath, int rows, int cols, double transform[])
{
    GDALDataset *pRaster;
    GDALDriver *pDriverTiff = GetGDALDriverManager()->GetDriverByName("GTiff");

    pRaster = pDriverTiff->Create(rasterPath, cols, rows, 1, GDT_Float32, NULL);
    pRaster->SetGeoTransform(transform);
    pRaster->GetRasterBand(1)->Fill(0.0);
    pRaster->GetRasterBand(1)->SetNoDataValue(-9999);

    GDALClose(pRaster);

    return 0;
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

    QString damElevName = "dam_elev";
    int nDamEl, nDemEl, row, col;
    QString demElevName = "elev";
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

double getDamHeight()
{
//    int randVal, high, low;
//    double returnVal;

//    high = 22;
//    low = 3;

//    randVal = (qrand() % ((high+1)-low)+low);
//    returnVal = randVal/10.0;

//    return returnVal;

//    double u1 = uniformRandom();
//    double u2 = uniformRandom();
//    return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));

    double x1, x2, w, y1, y2;

    do
    {
        x1 = 2.0 * uniformRandom() - 1.0;
        x2 = 2.0 * uniformRandom() - 1.0;
        w = x1 * x1 + x2 * x2;
    }
    while (w >= 1.0);

    //qDebug()<<"w"<<w<<"logw"<<log(w);
    w = sqrt((-2.0 * log(w)) / w);

    y1 = x1 * w;
    y2 = x2 * w;

    return exp(-0.0999+0.42*y1);
}

double getDamHeightLnorm()
{
    double returnVal;

    return returnVal;
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

double getRasterX(double transform[6], int col)
{
    double x;

    x = transform[0] + (col*transform[1] + (0.5*transform[1]));

    return x;
}

double getRasterY(double transform[6], int row)
{
    double y;

    y = transform[3] - ((row*fabs(transform[5])+ (0.5*fabs(transform[5]))));

    return y;
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

double getWetArea(const char *rasterPath)
{
    GDALDataset *pRaster;
    pRaster = (GDALDataset*) GDALOpen(rasterPath, GA_ReadOnly);

    int nCols = pRaster->GetRasterXSize();
    int nRows = pRaster->GetRasterYSize();
    double transform[6];
    double area;
    int count = 0;
    pRaster->GetGeoTransform(transform);

    float *value = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pRaster->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, value, nCols, 1, GDT_Float32, 0, 0);

        for (int j=0; j<nCols; j++)
        {
            if (value[j] > 0.0)
            {
                count++;
            }
        }
    }

    area = fabs(transform[1]*transform[5]) * (count*1.0);

    CPLFree(value);

    GDALClose(pRaster);

    return area;
}

int pointsInPolygon(const char *pointsPath, const char *polygonPath, const char *pointLayerName)
{
    OGRDataSource *pPointDS, *pPolyDS;
    OGRSFDriver *pDriverShp;
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    pDriverShp = registrar->GetDriverByName("ESRI Shapefile");

    pPointDS = pDriverShp->CreateDataSource(pointsPath);
    pPolyDS = pDriverShp->CreateDataSource(polygonPath);
    OGRLayer *pPointsLayer = pPointDS->GetLayerByName(pointLayerName);
    OGRLayer *pPolyLayer = pPolyDS->GetLayerByName("DamSearchPolygons");
    OGRLayer *pDamPointLayer = pPolyDS->CreateLayer("PondPts", pPointsLayer->GetSpatialRef(), wkbPoint, NULL);
    OGRFieldDefn field("dam_elev", OFTReal);
    pDamPointLayer->CreateField(&field);
    field.SetName("elev");
    field.SetType(OFTReal);
    pDamPointLayer->CreateField(&field);

    int nPolyCount = pPolyLayer->GetFeatureCount();
    double angle;
    for (int i=0; i<nPolyCount; i++)
    {
        OGRFeature *pPolyFeat = pPolyLayer->GetFeature(i);
        OGRPolygon *pPoly = (OGRPolygon*) pPolyFeat->GetGeometryRef();
        OGRLinearRing *pRing = pPoly->getExteriorRing();

        for (int j=0; j<pPointsLayer->GetFeatureCount(); j++)
        {
            OGRFeature *pPointFeat = pPointsLayer->GetFeature(j);
            OGRPoint *pPoint = (OGRPoint*) pPointFeat->GetGeometryRef();
            angle = 0.0;
            double px = pPoint->getX();
            double py = pPoint->getY();
            for (int k=0; k<pRing->getNumPoints()-1; k++)
            {
                angle += angleBetweenLines(pRing->getX(k), pRing->getY(k), pRing->getX(k+1), pRing->getY(k+1), px, py);
            }
            if (angle <PI)
            {
                //point not in polygon
            }
            else
            {
                double delev = pPolyFeat->GetFieldAsDouble("d_elev");
                double elev = pPointFeat->GetFieldAsDouble("GRID_CODE");
                OGRFeature *newFeature = OGRFeature::CreateFeature(pDamPointLayer->GetLayerDefn());
                newFeature->SetField("dam_elev", delev);
                newFeature->SetField("elev", elev);
                OGRPoint newPoint;
                newPoint.setX(pPoint->getX());
                newPoint.setY(pPoint->getY());
                newFeature->SetGeometry(&newPoint);
                pDamPointLayer->CreateFeature(newFeature);
                OGRFeature::DestroyFeature(newFeature);
            }

            OGRFeature::DestroyFeature(pPointFeat);
        }

        OGRFeature::DestroyFeature(pPolyFeat);
        //qDebug()<<"finished feature"<<i;
    }

    OGRDataSource::DestroyDataSource(pPointDS);
    OGRDataSource::DestroyDataSource(pPolyDS);

    return 0;
}

int pointsInPolygon2(const char *rasterPath, const char *rasterOut, const char *polygonPath)
{
    OGRDataSource *pPolyDS;
    OGRSFDriver *pDriverShp;
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    pDriverShp = registrar->GetDriverByName("ESRI Shapefile");

    GDALDataset *pRaster = (GDALDataset*) GDALOpen(rasterPath, GA_ReadOnly);
    double geot[6];
    pRaster->GetGeoTransform(geot);

    double noData = -9999;

    GDALDriver *pDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    GDALDataset *pRasterOut = pDriver->Create(rasterOut,pRaster->GetRasterXSize(), pRaster->GetRasterYSize(), 1, GDT_Float32, NULL);
    pRasterOut->SetGeoTransform(geot);
    pRasterOut->GetRasterBand(1)->Fill(noData);
    pRasterOut->GetRasterBand(1)->SetNoDataValue(noData);

    pPolyDS = pDriverShp->CreateDataSource(polygonPath);
    OGRLayer *pPolyLayer = pPolyDS->GetLayerByName("DamSearchPolygons");

    float *val = (float*) CPLMalloc(sizeof(float));

    int nPolyCount = pPolyLayer->GetFeatureCount();
    double angle;
    for (int i=0; i<nPolyCount; i++)
    {
        OGRFeature *pPolyFeat = pPolyLayer->GetFeature(i);
        OGRPolygon *pPoly = (OGRPolygon*) pPolyFeat->GetGeometryRef();
        OGRLinearRing *pRing = pPoly->getExteriorRing();

        OGREnvelope ringBound;
        pRing->getEnvelope(&ringBound);

        int left = getRasterCol(geot, ringBound.MinX)-1;
        int right = getRasterCol(geot, ringBound.MaxX)+1;
        int top = getRasterRow(geot, ringBound.MaxY)-1;
        int bottom = getRasterRow(geot, ringBound.MinY)+1;

        for (int i=top; i<bottom; i++)
        {
            for (int j=left; j<right; j++)
            {
                int x = getRasterX(geot, j);
                int y = getRasterY(geot, i);

                angle = 0.0;
                for (int k=0; k<pRing->getNumPoints()-1; k++)
                {
                    angle += angleBetweenLines(pRing->getX(k), pRing->getY(k), pRing->getX(k+1), pRing->getY(k+1), x, y);
                }
                if (angle <PI)
                {
                    //point not in polygon
                }
                else
                {
                    double delev = pPolyFeat->GetFieldAsDouble("d_elev");
                    pRaster->GetRasterBand(1)->RasterIO(GF_Read, j, i, 1, 1, val, 1, 1, GDT_Float32, 0, 0);
                    *val = delev - *val;
                    pRasterOut->GetRasterBand(1)->RasterIO(GF_Write, j, i, 1, 1, val, 1, 1, GDT_Float32, 0, 0);
                }
            }
        }

        OGRFeature::DestroyFeature(pPolyFeat);
        //qDebug()<<"finished feature"<<i;
    }

    OGRDataSource::DestroyDataSource(pPolyDS);
    CPLFree(val);
    GDALClose(pRaster);
    GDALClose(pRasterOut);

    return 0;
}

int regionsToDepth(const char *regRas, const char *depRas, const char *demRas, int nRegions)
{
    GDALDataset *pRegionsRaster, *pDepthRaster, *pDem;
    pRegionsRaster = (GDALDataset*) GDALOpen(regRas, GA_ReadOnly);
    pDepthRaster = (GDALDataset*) GDALOpen(depRas, GA_Update);
    pDem = (GDALDataset*) GDALOpen(demRas, GA_ReadOnly);
    pDepthRaster->GetRasterBand(1)->Fill(-9999);
    pDepthRaster->GetRasterBand(1)->SetNoDataValue(-9999);

    int nCols = pRegionsRaster->GetRasterXSize();
    int nRows = pRegionsRaster->GetRasterYSize();

    QVector<double> maxElev(nRegions, 0.0);

    float *regVal = (float*) CPLMalloc(sizeof(float)*nCols);
    float *depVal = (float*) CPLMalloc(sizeof(float)*nCols);
    float *demVal = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pRegionsRaster->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, regVal, nCols, 1, GDT_Float32, 0, 0);
        pDem->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, demVal, nCols, 1, GDT_Float32, 0, 0);
        for (int j=0; j<nCols; j++)
        {
            if (regVal[j] > 0)
            {
                int index = regVal[j] - 1;
                if (demVal[j] > maxElev[index])
                {
                    maxElev[index] = demVal[j];
                }
            }
        }
    }

    for (int i=0; i<nRows; i++)
    {
        pRegionsRaster->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, regVal, nCols, 1, GDT_Float32, 0, 0);
        pDem->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, demVal, nCols, 1, GDT_Float32, 0, 0);
        for (int j=0; j<nCols; j++)
        {
            if (regVal[j] > 0)
            {
                int index = regVal[j] - 1;
                depVal[j] = (maxElev[index] + 0.01) - demVal[j];
            }
            else
            {
                depVal[j] = -9999;
            }
        }
        pDepthRaster->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, depVal, nCols, 1, GDT_Float32, 0, 0);
    }

    CPLFree(regVal);
    CPLFree(depVal);
    CPLFree(demVal);

    GDALClose(pRegionsRaster);
    GDALClose(pDepthRaster);
    GDALClose(pDem);
}

double sampleRasterAlongLine_LowVal(const char * rasterPath, double startX, double startY, double azimuth, double distance, double &x, double &y)
{
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

    return lowValue;
}

int summarizeInundationRaster(const char *rasterPath, const char *outputCsv, int nValues, QVector<int> &thresholds, QVector<double> &areas)
{
    GDALDataset *pRaster;
    pRaster = (GDALDataset*) GDALOpen(rasterPath, GA_ReadOnly);
    double transform[6];
    pRaster->GetGeoTransform(transform);
    double cellArea = fabs(transform[1]*transform[5]);
    double area;
    int nRows = pRaster->GetRasterYSize();
    int nCols = pRaster->GetRasterXSize();
    float *value = (float*) CPLMalloc(sizeof(float)*nCols);
    QString csvPath = QString::fromUtf8(outputCsv);
    QFile file(csvPath);
    if (!file.open(QIODevice::ReadWrite))
    {
        qDebug()<<"error opening summary csv";
        return 1;
    }
    QTextStream stream(&file);
    stream<<"Value,Count,Area\n";

    for (int i=0; i<nValues; i++)
    {
        int cellCount = 0;

        for (int j=0; j<nRows; j++)
        {
            pRaster->GetRasterBand(1)->RasterIO(GF_Read, 0, j, nCols, 1, value, nCols, 1, GDT_Float32, 0, 0);
            for (int k=0; k<nCols; k++)
            {
                if (value[k] != -9999)
                {
                    if (value[k] > (i*1.0))
                    {
                        cellCount++;
                    }
                }
            }
        }
        area = (cellCount*1.0)*cellArea;
        stream<<QString::number(i+1)+","+QString::number(cellCount)+","+QString::number(area,'g',10);
        thresholds.append(i+1);
        areas.append(area);
        stream<<endl;
    }

    file.close();
    CPLFree(value);
    GDALClose(pRaster);

    return 0;
}

int updateInundationRaster(const char *updatePath, const char *inputPath)
{
    GDALDataset *pUpdate, *pInput;
    pUpdate = (GDALDataset*) GDALOpen(updatePath, GA_Update);
    pInput = (GDALDataset*) GDALOpen(inputPath, GA_ReadOnly);

    int nCols = pUpdate->GetRasterXSize();
    int nRows = pUpdate->GetRasterYSize();

    float *inputVal = (float*) CPLMalloc(sizeof(float)*nCols);
    float *updateVal = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pInput->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, inputVal, nCols, 1, GDT_Float32, 0, 0);
        pUpdate->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, updateVal, nCols, 1, GDT_Float32, 0, 0);

        for (int j=0; j<nCols; j++)
        {
            if (inputVal[j] > 0.0)
            {
                updateVal[j] += 1.0;
            }
        }

        pUpdate->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, updateVal, nCols, 1, GDT_Float32, 0, 0);
    }

    CPLFree(inputVal);
    CPLFree(updateVal);

    GDALClose(pInput);
    GDALClose(pUpdate);

    return 0;
}

double uniformRandom()
{
    return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}
