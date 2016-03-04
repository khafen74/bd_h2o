#include "raster.h"

Raster::Raster()
{

}

void Raster::add(const char *addPath)
{
    GDALDataset *pSourceDS, *pAddDS;

    pSourceDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_Update);
    pAddDS = (GDALDataset*) GDALOpen(addPath, GA_ReadOnly);

    float *srcRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *addRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *newRow = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pSourceDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, srcRow, nCols, 1, GDT_Float32, 0, 0);
        pAddDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, addRow, nCols, 1, GDT_Float32, 0, 0);

        for (int j=0; j<nCols; j++)
        {
            if (srcRow[j] == noData)
            {
                newRow[j] = noData;
            }
            else if (addRow[j] == noData)
            {
                newRow[j] = srcRow[j];
            }
            else
            {
                newRow[j] = srcRow[j] + addRow[j];
            }
        }

        pSourceDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, newRow, nCols, 1, GDT_Float32, 0, 0);
    }

    GDALClose(pSourceDS);
    GDALClose(pAddDS);

    CPLFree(srcRow);
    CPLFree(addRow);
    CPLFree(newRow);
}

void Raster::add(const char *sourcePath, const char *addPath)
{
    setProperties(sourcePath);
    add(addPath);
}

double Raster::area()
{
    GDALDataset *pRaster;

    pRaster  = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);

    float *value = (float*) CPLMalloc(sizeof(float)*nCols);

    int nCount = 0.0;

    for (int i=0; i<nRows; i++)
    {
        pRaster->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, value, nCols, 1, GDT_Float32, 0, 0);

        for (int j=0; j<nCols; j++)
        {
            if (value[j] != noData)
            {
                nCount++;
            }
        }
    }

    double area = nCount * fabs(transform[1]*transform[5]);

    CPLFree(value);
    GDALClose(pRaster);

    return area;
}

double Raster::area(const char *sourcePath)
{
    setProperties(sourcePath);
    return area();
}

void Raster::aspect(const char *aspectPath)
{
    GDALDataset *pSourceDS, *pAspectDS;

    pSourceDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);

    pAspectDS = pDriverTiff->Create(aspectPath, nCols, nRows, 1, GDT_Float32, NULL);
    pAspectDS->SetGeoTransform(transform);
    pAspectDS->GetRasterBand(1)->SetNoDataValue(noData);
    pAspectDS->GetRasterBand(1)->Fill(noData);

    bool calc;
    double dx, dy;

    float *eVals = (float*) CPLMalloc(sizeof(float)*9);
    float *aVal = (float*) CPLMalloc(sizeof(float)*1);

    for (int i=1; i<nRows-1; i++)
    {
        for (int j=1; j<nCols-1; j++)
        {
            pSourceDS->GetRasterBand(1)->RasterIO(GF_Read, j-1, i-1, 3, 3, eVals, 3, 3, GDT_Float32, 0, 0);

            calc = true;
            for (int k=0; k<9; k++)
            {
                if (eVals[k] == noData)
                {
                    *aVal = noData;
                    calc = false;
                    break;
                }
            }

            if (calc)
            {
                dx = ((eVals[2] + eVals[5] + eVals[5] + eVals[8]) - (eVals[0] + eVals[3] + eVals[3] + eVals[6]));
                dy = ((eVals[6] + eVals[7] + eVals[7] + eVals[8]) - (eVals[0] + eVals[1] + eVals[1] + eVals[2]));

                *aVal = atan2(dy/8.0, ((-1.0)*dx/8.0)) * (180.0/PI);

                if (dx == 0.0)
                {
                    if (*aVal < 0.0)
                    {
                        *aVal = 90.0 - *aVal;
                    }
                    else if (*aVal > 90.0)
                    {
                        *aVal = 360.0 - *aVal + 90.0;
                    }
                    else
                    {
                        *aVal = 90.0 - *aVal;
                    }
                }
                else
                {
                    if (*aVal > 90.0)
                    {
                        *aVal = 450.0 - *aVal;
                    }
                    else
                    {
                        *aVal = 90.0 - *aVal;
                    }
                }

                if (*aVal == 0.0)
                {
                    *aVal = 0.001;
                }
            }

            pAspectDS->GetRasterBand(1)->RasterIO(GF_Write, j, i, 1, 1, aVal, 1, 1, GDT_Float32, 0, 0);
        }
    }

    GDALClose(pSourceDS);
    GDALClose(pAspectDS);

    CPLFree(eVals);
    CPLFree(aVal);
}

void Raster::aspect(const char *sourcePath, const char *aspectPath)
{
    setProperties(sourcePath);

    aspect(aspectPath);
}

void Raster::extractByMask_CellCenters(const char *rasterOut, const char *polygonPath)
{
    OGRDataSource *pPolyDS;
    OGRSFDriver *pDriverShp;
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    pDriverShp = registrar->GetDriverByName("ESRI Shapefile");

    GDALDataset *pRaster = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
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

        int left = getCol(ringBound.MinX)-1;
        int right = getCol(ringBound.MaxX)+1;
        int top = getRow(ringBound.MaxY)-1;
        int bottom = getRow(ringBound.MinY)+1;

        for (int i=top; i<bottom; i++)
        {
            for (int j=left; j<right; j++)
            {
                int x = xCoordinate(j);
                int y = yCoordinate(i);

                angle = 0.0;
                for (int k=0; k<pRing->getNumPoints()-1; k++)
                {
                    angle += Geometry::angleBetweenLines(pRing->getX(k), pRing->getY(k), pRing->getX(k+1), pRing->getY(k+1), x, y);
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
}

void Raster::extractByMask_CellCenters(const char *rasterPath, const char *rasterOut, const char *polygonPath)
{
    setProperties(rasterPath);

    extractByMask_CellCenters(rasterOut, polygonPath);
}

void Raster::filterLowPass(const char *filterRaster)
{
    GDALDataset *pSource, *pFilter;

    pSource = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pFilter = pDriverTiff->Create(filterRaster, nCols, nRows, 1, GDT_Float32, NULL);

    pFilter->SetGeoTransform(transform);
    pFilter->GetRasterBand(1)->Fill(noData);
    pFilter->GetRasterBand(1)->SetNoDataValue(noData);

    int count;
    double average;

    float *value = (float*) CPLMalloc(sizeof(float)*1);
    float *window = (float*) CPLMalloc(sizeof(float)*9);

    for (int i=1; i<nRows-1; i++)
    {
        for (int j=1; j<nCols-1; j++)
        {
            count = 0;
            average = 0.0;

            pSource->GetRasterBand(1)->RasterIO(GF_Read, j-1, i-1, 3, 3, window, 3, 3, GDT_Float32, 0, 0);

            if (window[4] == noData)
            {
                *value = noData;
            }
            else
            {
                for (int k=0; k<9; k++)
                {
                    if (window[k] != noData)
                    {
                        count++;
                        average += window[k];
                    }
                }

                *value = average / (count*1.0);

            }

            pFilter->GetRasterBand(1)->RasterIO(GF_Write, j, i, 1, 1, value, 1, 1, GDT_Float32, 0, 0);

        }
    }

    GDALClose(pSource);
    GDALClose(pFilter);

    CPLFree(window);
    CPLFree(value);
}

void Raster::filterLowPass(const char *sourceRaster, const char *filterRaster)
{
    setProperties(sourceRaster);
    filterLowPass(filterRaster);
}

double Raster::findMax(const char *rasterPath)
{
    setProperties(rasterPath);

    GDALDataset *pRaster = (GDALDataset*) GDALOpen(rasterPath, GA_ReadOnly);

    double max = noData;
    bool firstVal = false;

    float *row = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pRaster->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, row, nCols, 1, GDT_Float32, 0, 0);

        for (int j=0; j<nCols; j++)
        {
            if (row[j] != noData)
            {
                if (firstVal)
                {
                    if (row[j] > max)
                    {
                        max = row[j];
                    }
                }
                else
                {
                    if (row[j] != noData)
                    {
                        max = row[j];
                        firstVal = true;
                    }
                }
            }
        }
    }

    GDALClose(pRaster);

    CPLFree(row);

    return max;
}

void Raster::fromXYZ(const char *rasterPath, const char *xyzPath, int cols, int rows, double noDataValue, double inTransform[6], int headerRows)
{
    int i, j;
    double x, y, yTLCenter, xTLCenter;
    float z;
    QString qsDummy, qsX, qsY, qsZ;

    GDALDataset *pDatasetNew;
    pDatasetNew = pDriverTiff->Create(rasterPath, cols, rows, 1, GDT_Float32, NULL);
    pDatasetNew->SetGeoTransform(inTransform);
    pDatasetNew->GetRasterBand(1)->Fill(noDataValue);
    pDatasetNew->GetRasterBand(1)->SetNoDataValue(noDataValue);

    xTLCenter = inTransform[0] + (inTransform[1] / 2.0);
    yTLCenter = inTransform[3] - (inTransform[1] / 2.0);

    float *rasVal = (float*) CPLMalloc(sizeof(float)*1);

    QFile inFile (QString::fromUtf8(xyzPath));

    if (inFile.open(QIODevice::ReadOnly | QIODevice::Text))
    {

        QTextStream stream(&inFile);

        int count = 0;

        while(!stream.atEnd())
        {
            if (count < headerRows)
            {
                qsDummy = stream.readLine();
                count++;
            }
            else
            {
                stream >> qsX;
                x = qsX.toDouble();
                stream >> qsY;
                y = qsY.toDouble();
                stream >> qsZ;
                z = qsZ.toDouble();

                *rasVal = z;

                i = (yTLCenter - y) / inTransform[1];
                j = (x - xTLCenter) / inTransform[1];

                if ((i>=0 && i<rows) && (j>=0 && j<cols))
                {
                    pDatasetNew->GetRasterBand(1)->RasterIO(GF_Write, j, i, 1, 1, rasVal, 1, 1, GDT_Float32, 0, 0);
                }

                count++;

            }
        }
    }

    setProperties(rasterPath);

    GDALClose(pDatasetNew);
    CPLFree(rasVal);
}

int Raster::getCol(double xCoord)
{
    int col;
    double xOffset, xDiv;

    xOffset = xCoord - transform[0];
    xDiv = xOffset/transform[1];
    col = floor(xDiv);

    return col;
}

int Raster::getCols()
{
    return nCols;
}

const char *Raster::getPath()
{
    return m_rasterPath.toStdString().c_str();
}

double Raster::getRow(double yCoord)
{
    int row;

    double yOffset, yDiv;

    yOffset = transform[3] - yCoord;
    yDiv = yOffset/fabs(transform[5]);
    row = floor(yDiv);

    return row;
}

int Raster::getRows()
{
    return nRows;
}

void Raster::greaterThan(const char *outPath, double value)
{
    GDALDataset *pSourceDS, *pGreaterDS;

    pSourceDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pGreaterDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);
    pGreaterDS->SetGeoTransform(transform);
    pGreaterDS->GetRasterBand(1)->SetNoDataValue(noData);
    pGreaterDS->GetRasterBand(1)->Fill(noData);

    float *oldVal = (float*) CPLMalloc(sizeof(float)*nCols);
    float *newVal = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pSourceDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, oldVal, nCols, 1, GDT_Float32, 0, 0);
        for (int j=0; j<nCols; j++)
        {
            if (oldVal[j] > value)
            {
                newVal[j] = 0.0;
            }
            else
            {
                newVal[j] = noData;
            }
        }
        pGreaterDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, newVal, nCols, 1, GDT_Float32, 0, 0);
    }

    CPLFree(oldVal);
    CPLFree(newVal);

    GDALClose(pSourceDS);
    GDALClose(pGreaterDS);
}

void Raster::greaterThan(const char *inPath, const char *outPath, double value)
{
    setProperties(inPath);
    greaterThan(outPath, value);
}

void Raster::hillshade(const char *hlsdPath)
{
    GDALDataset *pSourceDS, *pHlsdDS;

    pSourceDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pHlsdDS = pDriverTiff->Create(hlsdPath, nCols, nRows, 1, GDT_Byte, NULL);
    pHlsdDS->SetGeoTransform(transform);
    pHlsdDS->GetRasterBand(1)->SetNoDataValue(0);
    pHlsdDS->GetRasterBand(1)->Fill(0);

    double altDeg = 45.0, azimuth = 315.0, zFactor = 1.0,
            zenDeg = 90.0 - altDeg, zenRad = zenDeg *PI / 180.0,
            azimuthMath = 360.0 - azimuth + 90.0,
            azimuthRad, dzdx, dzdy, slopeRad, aspectRad, hlsdByte;

    if (azimuthMath > 360.0)
    {
        azimuthMath = azimuthMath - 360.0;
    }

    azimuthRad = azimuthMath * PI / 180.0;

    float *elevWin = (float*) CPLMalloc(sizeof(float)*9);
    unsigned char *hlsdRow = (unsigned char*) CPLMalloc(sizeof(int)*nCols);

    bool calculate;

    for (int i=1; i<nRows-1; i++)
    {
        hlsdRow[0] = 0, hlsdRow[nCols-1] = 0;

        for (int j=1; j<nCols-1; j++)
        {
            pSourceDS->GetRasterBand(1)->RasterIO(GF_Read, j-1, i-1, 3, 3, elevWin, 3, 3, GDT_Float32, 0, 0);

            calculate = true;

            for (int k=0; k<9; k++)
            {
                if (elevWin[k] == noData)
                {
                    calculate = false;
                }
            }

            if (calculate)
            {
                dzdx = ((elevWin[2]+(2*elevWin[5])+elevWin[8]) - (elevWin[0]+(2*elevWin[3])+elevWin[6])) / (8*transform[1]);
                dzdy = ((elevWin[6]+(2*elevWin[7])+elevWin[8]) - (elevWin[0]+(2*elevWin[1])+elevWin[2])) / (8*transform[1]);
                slopeRad = atan(zFactor * sqrt(pow(dzdx,2)+pow(dzdy,2)));

                if (dzdx != 0.0)
                {
                    aspectRad = atan2(dzdy, (dzdx * (-1.0)));

                    if ( aspectRad < 0.0)
                    {
                        aspectRad = 2.0 * PI + aspectRad;
                    }
                }
                else
                {
                    if (dzdy > 0.0)
                    {
                        aspectRad = PI / 2.0;
                    }
                    else if (dzdy < 0.0)
                    {
                        aspectRad = (2.0 * PI) - (PI / 2.0);
                    }
                    else
                    {
                        aspectRad = aspectRad;
                    }
                }

                hlsdByte = round(254 * ((cos(zenRad) * cos(slopeRad)) + (sin(zenRad) * sin(slopeRad) * cos(azimuthRad - aspectRad)))) + 1.0;
                hlsdRow[j] = hlsdByte;

            }
            else
            {
                hlsdRow[j] = 0;
            }
        }

        pHlsdDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, hlsdRow, nCols, 1, GDT_Byte, 0, 0);

    }

    GDALClose(pHlsdDS);
    GDALClose(pSourceDS);

    CPLFree(elevWin);
    CPLFree(hlsdRow);
}

void Raster::hillshade(const char *rasterPath, const char *hlsdPath)
{
    setProperties(rasterPath);
    hillshade(hlsdPath);
}

int Raster::regions(const char *regionsRaster)
{
    GDALDataset *pInputRaster, *pRegionsRaster;

    pInputRaster = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);

    pRegionsRaster = pDriverTiff->CreateCopy(regionsRaster, pInputRaster, FALSE, NULL, NULL, NULL);

    float *inpVal = (float*) CPLMalloc(sizeof(float)*1);
    float *regVal = (float*) CPLMalloc(sizeof(float)*1);
    float *inpWin = (float*) CPLMalloc(sizeof(float)*9);
    float *regWin = (float*) CPLMalloc(sizeof(float)*9);

    int regionCount, changedCount, regionValue, valueCount, totalCount;
    bool regionDone, rasterDone;

    rasterDone = false;
    regionValue = 0;
    totalCount = 0;

    while (!rasterDone && totalCount < 4)
    {
        changedCount = 0;

        for (int i=1; i<nRows-1; i++)
        {
            for (int j=1; j<nCols-1; j++)
            {
                pInputRaster->GetRasterBand(1)->RasterIO(GF_Read, j, i, 1, 1, inpVal, 1, 1, GDT_Float32, 0, 0);
                pRegionsRaster->GetRasterBand(1)->RasterIO(GF_Read, j, i, 1, 1, regVal, 1, 1, GDT_Float32, 0, 0);

                if (*inpVal == 0 && *regVal == 0)
                {
                    regionValue++;
                    changedCount++;
                    *regVal = regionValue;
                    pRegionsRaster->GetRasterBand(1)->RasterIO(GF_Write, j, i, 1, 1, regVal, 1, 1, GDT_Float32, 0, 0);
                    regionDone = false;

                    valueCount = 0;

                    while (!regionDone)
                    {
                        regionCount = 0;

                        for (int k=1; k<nRows-1; k++)
                        {
                            for (int l=1; l<nCols-1; l++)
                            {
                                pRegionsRaster->GetRasterBand(1)->RasterIO(GF_Read, l-1, k-1, 3, 3, regWin, 3, 3, GDT_Float32, 0, 0);
                                pInputRaster->GetRasterBand(1)->RasterIO(GF_Read, l-1, k-1, 3, 3, inpWin, 3, 3, GDT_Float32, 0, 0);

                                if (inpWin[4] == 0.0 && regWin[4] == *regVal)
                                {
                                    for (int m=0; m<9; m++)
                                    {
                                        if (regWin[m] == 0.0)
                                        {
                                            pRegionsRaster->GetRasterBand(1)->RasterIO(GF_Write, l+COL_OFFSET[m], k+ROW_OFFSET[m], 1, 1, regVal, 1, 1, GDT_Float32, 0, 0);
                                            regionCount++;
                                            valueCount++;
                                        }
                                    }
                                }
                            }
                        }

                        if (regionCount == 0)
                        {
                            regionDone = true;
                        }
                    }
                }
            }
        }

        if (changedCount == 0)
        {
            rasterDone = true;
        }

        totalCount++;
    }

    GDALClose(pInputRaster);
    GDALClose(pRegionsRaster);

    CPLFree(inpVal);
    CPLFree(regVal);
    CPLFree(inpWin);
    CPLFree(regWin);

    return regionValue;
}

int Raster::regions(const char *inputRaster, const char *regionsRaster)
{
    int regionsValue;
    setProperties(inputRaster);
    regionsValue = regions(regionsRaster);

    return regionsValue;
}

double Raster::sampleAlongLine_LowVal(double startX, double startY, double azimuth, double distance, double &x, double &y)
{
    double az1, az2, interval;
    double transform[6];
    double newX, newY, rasValue, lowValue;
    int nSamples;
    az1 = Geometry::addDegrees(azimuth, 90.0);
    az2 = Geometry::addDegrees(azimuth, -90.0);

    //qDebug()<<"opening raster";
    GDALDataset *pRas;
    pRas = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pRas->GetGeoTransform(transform);
    GDALClose(pRas);
    //qDebug()<<"raster closed";

    interval = transform[1]/5.0;
    nSamples = ceil(distance/interval);

    //qDebug()<<"starting loop";
    for (int i=0; i<nSamples; i++)
    {
        Geometry::calcCoords(startX, startY, az1, interval*(i+1), newX, newY);
        rasValue = valueAtPoint(newX, newY);
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
        Geometry::calcCoords(startX, startY, az2, interval*(i+1), newX, newY);
        rasValue = valueAtPoint(newX, newY);
        if (rasValue < lowValue)
        {
            lowValue = rasValue;
            x = newX, y = newY;
        }
    }

    return lowValue;
}

double Raster::sampleAlongLine_LowVal(const char *rasterPath, double startX, double startY, double azimuth, double distance, double &x, double &y)
{
    setProperties(rasterPath);

    double value = sampleAlongLine_LowVal(startX, startY, azimuth, distance, x, y);

    return value;
}

void Raster::setProperties(const char *rasterPath)
{
    loadDrivers();

    m_rasterPath = QString::fromUtf8(rasterPath);
    GDALDataset *pRaster = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);

    nRows = pRaster->GetRasterBand(1)->GetYSize();
    nCols = pRaster->GetRasterBand(1)->GetXSize();
    pRaster->GetGeoTransform(transform);
    noData = pRaster->GetRasterBand(1)->GetNoDataValue();

    GDALClose(pRaster);
}

void Raster::slopeTOF(const char *slopePath)
{
    GDALDataset *pSourceDS, *pSlopeDS;

    pSourceDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pSlopeDS = pDriverTiff->Create(slopePath, nCols, nRows, 1, GDT_Float32, NULL);
    pSlopeDS->SetGeoTransform(transform);
    pSlopeDS->GetRasterBand(1)->SetNoDataValue(noData);
    pSlopeDS->GetRasterBand(1)->Fill(noData);

    double xslope, yslope, xyslope, xypow;

    float* eVals = (float*) CPLMalloc(sizeof(float)*9);
    float* sVals = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows-1; i++)
    {
        for (int j=0; j<nCols-1; j++)
        {
            if (i == 0 || i == nRows || j == 0 || j == nCols)
            {
                sVals[j] = noData;
            }
            else
            {
                pSourceDS->GetRasterBand(1)->RasterIO(GF_Read,j-1,i-1,3,3,eVals,3,3,GDT_Float32,0,0);
                if (eVals[4] == noData || eVals[0] == noData || eVals[1] == noData || eVals[2] == noData || eVals[3] == noData || eVals[5] == noData || eVals[6] == noData || eVals[7] == noData || eVals[8] == noData)
                {
                    sVals[j] = noData;
                }
                else
                {
                    xslope = ((eVals[2]-eVals[0]) + ((2*eVals[5])-(2*eVals[3])) + (eVals[8]-eVals[6])) / (8*transform[1]);
                    yslope = ((eVals[0]-eVals[6]) + ((2*eVals[1])-(2*eVals[7])) + (eVals[2]-eVals[8]))/(8*transform[1]);
                    xyslope = pow(xslope,2.0) + pow(yslope,2.0);
                    xypow = pow(xyslope,0.5);
                    sVals[j] = (atan(xypow)*180.0/PI);
                }
            }
        }
        pSlopeDS->GetRasterBand(1)->RasterIO(GF_Write,0,i,nCols,1,sVals,nCols,1,GDT_Float32,0,0);
    }

    GDALClose(pSourceDS);
    GDALClose(pSlopeDS);

    CPLFree(eVals);
    CPLFree(sVals);
}

void Raster::slopeTOF(const char *sourcePath, const char *slopePath)
{
    setProperties(sourcePath);

    slopeTOF(slopePath);
}

void Raster::subtract(const char *subtractPath)
{
    GDALDataset *pSourceDS, *pSubtractDS;

    loadDrivers();

    pSourceDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_Update);
    pSubtractDS = (GDALDataset*) GDALOpen(subtractPath, GA_ReadOnly);

    float *srcRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *subRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *newRow = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pSourceDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, srcRow, nCols, 1, GDT_Float32, 0, 0);
        pSubtractDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, subRow, nCols, 1, GDT_Float32, 0, 0);

        for (int j=0; j<nCols; j++)
        {
            if (srcRow[j] == noData || subRow[j] == noData)
            {
                newRow[j] = noData;
            }
            else if (subRow[j] == noData)
            {
                newRow[j] = noData;
            }
            else
            {
                newRow[j] = srcRow[j] - subRow[j];
            }
        }

        pSourceDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, newRow, nCols, 1, GDT_Float32, 0, 0);
    }

    GDALClose(pSourceDS);
    GDALClose(pSubtractDS);

    CPLFree(srcRow);
    CPLFree(subRow);
    CPLFree(newRow);
}

void Raster::subtract(const char *sourcePath, const char *subtractPath)
{
    setProperties(sourcePath);
    subtract(subtractPath);
}

void Raster::subtract(const char *sourcePath, const char *subtractPath, const char *outputPath)
{
    setProperties(sourcePath);

    GDALDataset *pSourceDS, *pSubtractDS, *pOutDS;

    pSourceDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_Update);
    pSubtractDS = (GDALDataset*) GDALOpen(subtractPath, GA_ReadOnly);
    pOutDS = pDriverTiff->Create(outputPath, nCols, nRows, 1, GDT_Float32, NULL);

    pOutDS->GetRasterBand(1)->Fill(noData);
    pOutDS->GetRasterBand(1)->SetNoDataValue(noData);
    pOutDS->SetGeoTransform(transform);

    float *srcRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *subRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *newRow = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pSourceDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, srcRow, nCols, 1, GDT_Float32, 0, 0);
        pSubtractDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, subRow, nCols, 1, GDT_Float32, 0, 0);

        for (int j=0; j<nCols; j++)
        {
            if (srcRow[j] == noData)
            {
                newRow[j] = noData;
            }
            else if (subRow[j] == noData)
            {
                newRow[j] = noData;
            }
            else
            {
                newRow[j] = srcRow[j] - subRow[j];
                if (newRow[j] < -9990 || newRow[j] > 9990)
                {
                    newRow[j] = noData;
                }
            }
        }

        pOutDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, newRow, nCols, 1, GDT_Float32, 0, 0);
    }

    GDALClose(pSourceDS);
    GDALClose(pSubtractDS);
    GDALClose(pOutDS);

    CPLFree(srcRow);
    CPLFree(subRow);
    CPLFree(newRow);
}

double Raster::sum()
{
    GDALDataset *pRaster;

    pRaster = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);

    double dSum = 0.0;

    float *row = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pRaster->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, row, nCols, 1, GDT_Float32, 0, 0);

        for (int j=0; j<nCols; j++)
        {
            if (row[j] != noData)
            {
                dSum += row[j];
            }
        }
    }

    GDALClose(pRaster);

    CPLFree(row);

    return dSum;
}

double Raster::sum(const char *rasterPath)
{
    setProperties(rasterPath);

    double dSum = 0.0;
    dSum = sum();

    return dSum;
}

double Raster::valueAtPoint(double xCoord, double yCoord)
{
    GDALDataset *pRaster;
    double transform[6];
    double value, xOffset, yOffset, xDiv, yDiv;
    int row, col;

    pRaster = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);

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

double Raster::valueAtPoint(const char *rasterPath, double xCoord, double yCoord)
{
    setProperties(rasterPath);

    double value = valueAtPoint(xCoord, yCoord);

    return value;
}

void Raster::writeCellValue(double xCoord, double yCoord, double value)
{
    float *newVal = (float*) CPLMalloc(sizeof(float)*1);
    *newVal = value;
    int row = getRow(yCoord);
    int col = getCol(xCoord);
    GDALDataset *pRaster;
    pRaster = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_Update);
    pRaster->GetRasterBand(1)->RasterIO(GF_Write, col, row, 1, 1, newVal, 1, 1, GDT_Float32, 0, 0);
    CPLFree(newVal);
    GDALClose(pRaster);
}

void Raster::writeCellValue(const char *rasterPath, double xCoord, double yCoord, double value)
{
    setProperties(rasterPath);
    writeCellValue(xCoord, yCoord, value);
}

double Raster::xCoordinate(int col)
{
    double x;

    x = transform[0] + (col*transform[1] + (0.5*transform[1]));

    return x;
}

double Raster::yCoordinate(int row)
{
    double y;

    y = transform[3] - ((row*fabs(transform[5])+ (0.5*fabs(transform[5]))));

    return y;
}

void Raster::zeroToNoData(const char *sourcePath, double noDataValue)
{
    setProperties(sourcePath);
    noData = noDataValue;

    GDALDataset *pSourceRaster;

    pSourceRaster = (GDALDataset*) GDALOpen(sourcePath, GA_Update);

    float *oldRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *newRow = (float*) CPLMalloc(sizeof(float)*nCols);

    GDALClose(pSourceRaster);

    for (int i=0; i<nRows; i++)
    {
        pSourceRaster->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, oldRow, nCols, 1, GDT_Float32, 0, 0);

        for (int j=0; j<nCols; j++)
        {
            if (oldRow[j] == 0.0 || oldRow[j] == noData)
            {
                newRow[j] = noData;
            }
            else
            {
                newRow[j] = oldRow[j];
            }
        }

        pSourceRaster->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, newRow, nCols, 1, GDT_Float32, 0, 0);
    }

    GDALClose(pSourceRaster);

    CPLFree(oldRow);
    CPLFree(newRow);
}

void Raster::loadDrivers()
{
    GDALAllRegister();

    pDriverTiff = GetGDALDriverManager()->GetDriverByName("GTiff");
}
