#include "raster.h"

Raster::Raster()
{

}

void Raster::adjustSoil(const char *mukeyPath, const char *varPath, const char *outPath)
{
    setProperties(mukeyPath);
    GDALDataset *mukDS, *varDS, *outDS;
    double transform[6];

    mukDS = (GDALDataset*) GDALOpen(mukeyPath, GA_ReadOnly);
    varDS = (GDALDataset*) GDALOpen(varPath, GA_ReadOnly);

    outDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);

    float *mukVal = (float*) CPLMalloc(sizeof(float)*nCols);
    float *varVal = (float*) CPLMalloc(sizeof(float)*nCols);
    float *outVal = (float*) CPLMalloc(sizeof(float)*nCols);

    double mukNoData = mukDS->GetRasterBand(1)->GetNoDataValue();
    double varNoData = varDS->GetRasterBand(1)->GetNoDataValue();
    qDebug()<<mukNoData<<varNoData;
    qDebug()<<nCols<<nRows;

    mukDS->GetGeoTransform(transform);
    outDS->SetGeoTransform(transform);
    outDS->SetProjection(mukDS->GetProjectionRef());
    outDS->GetRasterBand(1)->SetNoDataValue(-9999.0);
    qDebug()<<"starting loop";
    for (int i=0; i<nRows; i++)
    {
        mukDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, mukVal, nCols, 1, GDT_Float32, 0, 0);
        varDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, varVal, nCols, 1, GDT_Float32, 0, 0);
        for (int j=0; j<nCols; j++)
        {
            if (mukVal[j]>=2147483648.0 || mukVal<0)
            {
                outVal[j] = -9999.0;
            }
            else if (varVal[j]==varNoData || varVal[j] <= 0.0)
            {
                outVal[j] = 0.0;
            }
            else
            {
                outVal[j] = varVal[j];
                //qDebug()<<outVal[j];
            }
            if (varVal[j] > 0.0 && varVal[j] < 200.0)
            {
                //qDebug()<<"outval"<<outVal[j];
            }
        }
        outDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, outVal, nCols, 1, GDT_Float32, 0, 0);
        if (i%500 == 0)
        {
            qDebug()<<"row "<<i+1<<" of "<<nRows;
        }
    }
    qDebug()<<"done";

    CPLFree(mukVal);
    CPLFree(varVal);
    CPLFree(outVal);

    GDALClose(outDS);
    GDALClose(varDS);
    GDALClose(mukDS);
}
void Raster::add(const char *addPath, const char *outPath)
{
    GDALDataset *pSourceDS, *pAddDs, *pOutDS;
    //qDebug()<<m_rasterPath<<addPath<<outPath;
    pSourceDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pAddDs = (GDALDataset*) GDALOpen(addPath, GA_ReadOnly);
    pOutDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);
    pOutDS->SetGeoTransform(transform);
    pOutDS->GetRasterBand(1)->SetNoDataValue(noData);
    pOutDS->GetRasterBand(1)->Fill(noData);

    float *srcRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *addRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *outRow = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pSourceDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, srcRow, nCols, 1, GDT_Float32, 0, 0);
        pAddDs->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, addRow, nCols, 1, GDT_Float32, 0, 0);

        for (int j=0; j<nCols; j++)
        {
            if (srcRow[j] == noData)
            {
                outRow[j] = noData;
            }
            else if (addRow[j] == noData)
            {
                outRow[j] = srcRow[j];
            }
            else
            {
                outRow[j] = srcRow[j] + addRow[j];
            }
        }

        pOutDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, outRow, nCols, 1, GDT_Float32, 0, 0);
    }

    CPLFree(srcRow);
    CPLFree(addRow);
    CPLFree(outRow);

    GDALClose(pSourceDS);
    GDALClose(pAddDs);
    GDALClose(pOutDS);
}

void Raster::add(const char *sourcePath, const char *addPath, const char *outPath)
{
    setProperties(sourcePath);
    add(addPath, outPath);
}

void Raster::addTo(const char *addPath)
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

void Raster::addTo(const char *sourcePath, const char *addPath)
{
    setProperties(sourcePath);
    addTo(addPath);
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

int Raster::checkRowCol(int row, int col)
{
    if (row < 1 || row >= nRows-1)
    {
        //qDebug()<<"row value out of range "<<row;
        return 0;
    }
    if (col < 0 || col >= nCols-1)
    {
        //qDebug()<<"column out of range "<<col;
        return 0;
    }
    return 1;
}

bool Raster::drainsToMe(int index, int fdir)
{
    if (index == 4)
    {
        return false;
    }
    else if (index == 0 && fdir == FLOW_DIR[8])
    {
        return true;
    }
    else if (index == 1 && fdir == FLOW_DIR[7])
    {
        return true;
    }
    else if (index == 2 && fdir == FLOW_DIR[6])
    {
        return true;
    }
    else if (index == 3 && fdir == FLOW_DIR[5])
    {
        return true;
    }
    else if (index == 5 && fdir == FLOW_DIR[3])
    {
        return true;
    }
    else if (index == 6 && fdir == FLOW_DIR[2])
    {
        return true;
    }
    else if (index == 7 && fdir == FLOW_DIR[1])
    {
        return true;
    }
    else if (index == 8 && fdir == FLOW_DIR[0])
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Raster::extractByMask_CellCenters(const char *rasterOut, const char *polygonPath, const char *lyrName)
{
    OGRDataSource *pPolyDS;
    OGRSFDriver *pDriverShp;
    OGRRegisterAll();
    OGRSFDriverRegistrar *registrar = OGRSFDriverRegistrar::GetRegistrar();
    pDriverShp = registrar->GetDriverByName("ESRI Shapefile");

    GDALDataset *pRaster = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    double geot[6];
    pRaster->GetGeoTransform(geot);

    double noData = -9999;

    GDALDataset *pRasterOut = pDriverTiff->Create(rasterOut,pRaster->GetRasterXSize(), pRaster->GetRasterYSize(), 1, GDT_Float32, NULL);
    pRasterOut->SetGeoTransform(geot);
    pRasterOut->GetRasterBand(1)->Fill(noData);
    pRasterOut->GetRasterBand(1)->SetNoDataValue(noData);

    pPolyDS = pDriverShp->CreateDataSource(polygonPath);
    OGRLayer *pPolyLayer = pPolyDS->GetLayerByName(lyrName);

    float *val = (float*) CPLMalloc(sizeof(float)*1);

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

void Raster::extractByMask_CellCenters(const char *rasterPath, const char *rasterOut, const char *polygonPath, const char *lyrName)
{
    setProperties(rasterPath);

    extractByMask_CellCenters(rasterOut, polygonPath, lyrName);
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

void Raster::fromXYZ(const char *rasterPath, const char *xyzPath, int cols, int rows, double noDataValue, int headerRows)
{
    loadDrivers();
    double x, y;
    float pred, lwr, upr;
    double inTransform[6] = {1000.0, 1.0, 0.0, 1000.0, 0.0, -1.0};
    QString qsDummy, qsX, qsY, qsPred, qsLwr, qsUpr;

    GDALDataset *pDatasetNew;
    //3 band raster
    pDatasetNew = pDriverTiff->Create(rasterPath, cols, rows, 3, GDT_Float32, NULL);
    pDatasetNew->SetGeoTransform(inTransform);
    pDatasetNew->GetRasterBand(1)->Fill(noDataValue);
    pDatasetNew->GetRasterBand(1)->SetNoDataValue(noDataValue);
    pDatasetNew->GetRasterBand(2)->Fill(noDataValue);
    pDatasetNew->GetRasterBand(2)->SetNoDataValue(noDataValue);
    pDatasetNew->GetRasterBand(3)->Fill(noDataValue);
    pDatasetNew->GetRasterBand(3)->SetNoDataValue(noDataValue);



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
                x = qsX.toInt();
                stream >> qsY;
                y = qsY.toInt();
                stream >> qsPred;
                pred = qsPred.toDouble();
                stream >> qsLwr;
                lwr = qsLwr.toDouble();
                stream >> qsUpr;
                upr = qsUpr.toDouble();

                *rasVal = pred;
                qDebug()<<x<<y;

                if ((y>=0 && y<rows) && (x>=0 && x<cols))
                {
                    pDatasetNew->GetRasterBand(1)->RasterIO(GF_Write, x-1, y-1, 1, 1, rasVal, 1, 1, GDT_Float32, 0, 0);
                    *rasVal = lwr;
                    pDatasetNew->GetRasterBand(2)->RasterIO(GF_Write, x-1, y-1, 1, 1, rasVal, 1, 1, GDT_Float32, 0, 0);
                    *rasVal = upr;
                    pDatasetNew->GetRasterBand(3)->RasterIO(GF_Write, x-1, y-1, 1, 1, rasVal, 1, 1, GDT_Float32, 0, 0);
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

int Raster::getD8Index(int nFdir)
{
    if (nFdir == 32)
    {
        return 0;
    }
    else if (nFdir == 64)
    {
        return 1;
    }
    else if (nFdir == 128)
    {
        return 2;
    }
    else if (nFdir == 16)
    {
        return 3;
    }
    else if (nFdir == 1)
    {
        return 5;
    }
    else if (nFdir == 8)
    {
        return 6;
    }
    else if (nFdir == 4)
    {
        return 7;
    }
    else if (nFdir == 2)
    {
        return 8;
    }
    else
    {
        return -1;
    }
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

void Raster::heightAboveNetwork(const char *fdirPath, const char *facPath, const char *outPath)
{
    GDALDataset *pDemDS, *pFdirDS, *pFacDS, *pOutDS;

    pDemDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pFdirDS = (GDALDataset*) GDALOpen(fdirPath, GA_ReadOnly);
    pFacDS = (GDALDataset*) GDALOpen(facPath, GA_ReadOnly);
    pOutDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);
    pOutDS->SetGeoTransform(transform);
    pOutDS->GetRasterBand(1)->SetNoDataValue(noData);
    pOutDS->GetRasterBand(1)->Fill(noData);
    //qDebug()<<"no data value "<<noData;

    int nIndex, startRow, startCol, newRow, newCol;
    QVector<QString> indices;
    bool done, write;
    unsigned char *fdirWin = (unsigned char*) CPLMalloc(sizeof(unsigned char)*9);
    signed long int *facWin = (signed long int*) CPLMalloc(sizeof(signed long int)*9);
    float *elevValStart = (float*) CPLMalloc(sizeof(float)*1);
    float *elevVal = (float*) CPLMalloc(sizeof(float)*1);
    //qDebug()<<"starting HAND loop";

    for (int i=1; i<nRows-1; i++)
    {
        for (int j=1; j<nCols-1; j++)
        {
            pDemDS->GetRasterBand(1)->RasterIO(GF_Read, j, i, 1, 1, elevValStart, 1, 1, GDT_Float32, 0, 0);
            startRow = i, startCol = j, newRow = i, newCol = j;

            done = false, write = false;
            int nCount = 0;
            while (!done && nCount<2500)
            {
                //qDebug()<<newRow<<newCol;
                indices.clear();
                pFdirDS->GetRasterBand(1)->RasterIO(GF_Read, newCol-1, newRow-1, 3, 3, fdirWin, 3, 3, GDT_Byte, 0, 0);
                pFacDS->GetRasterBand(1)->RasterIO(GF_Read, newCol-1, newRow-1, 3, 3, facWin, 3, 3, GDT_Int32, 0, 0);

                if (facWin[4] > 0)
                {
                    //qDebug()<<"fac value > 0 "<<facWin[4];
                    write = true;
                    nIndex = 4;
                }
                else if (fdirWin[4] > 0)
                {
                    //qDebug()<<"flow dir "<<fdirWin[4];
                    nIndex = getD8Index(fdirWin[4]);

                    if (checkRowCol(newRow+ROW_OFFSET[nIndex], newCol+COL_OFFSET[nIndex]))
                    {
                        indices.append(QString::number(newRow) + " "+ QString::number(newCol));
                        newRow += ROW_OFFSET[nIndex], newCol += COL_OFFSET[nIndex];
                        if (indices.indexOf(QString::number(newRow) + " " + QString::number(newCol)) != -1)
                        {
                            qDebug()<<"matching index";
                            done = true;
                        }
                        if (facWin[nIndex] != noData)
                        {
                            write = true;
                        }
                    }
                    else
                    {
                        done = true;
                    }
                }
                else
                {
                    done = true;
                }

                if (write)
                {
                    pDemDS->GetRasterBand(1)->RasterIO(GF_Read, newCol, newRow, 1, 1, elevVal, 1, 1, GDT_Float32, 0, 0);
                    *elevVal = *elevValStart - *elevVal;
                    pOutDS->GetRasterBand(1)->RasterIO(GF_Write, startCol, startRow, 1, 1, elevVal, 1, 1, GDT_Float32, 0, 0);
                    done = true;
                }
                nCount++;
            }
            //qDebug()<<"col "<<j+1<<" of "<<nRows<<" completed";
        }
        //qDebug()<<"row "<<i+1<<" of "<<nRows<<" completed";
    }


    CPLFree(fdirWin);
    CPLFree(facWin);
    CPLFree(elevValStart);
    CPLFree(elevVal);

    GDALClose(pDemDS);
    GDALClose(pFdirDS);
    GDALClose(pFacDS);
    GDALClose(pOutDS);
}

void Raster::heightAboveNetwork(const char *demPath, const char *fdirPath, const char *facPath, const char *outPath)
{
    setProperties(demPath);
    heightAboveNetwork(fdirPath, facPath, outPath);
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
            if (rasValue < lowValue && rasValue != noData)
            {
                lowValue = rasValue;
                x = newX, y = newY;
            }
        }
        Geometry::calcCoords(startX, startY, az2, interval*(i+1), newX, newY);
        rasValue = valueAtPoint(newX, newY);
        if (rasValue < lowValue && rasValue != noData)
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

double Raster::sampleAlongLine_RasterVal(const char *checkRasPath, double startX, double startY, double azimuth, double distance, double &x, double &y)
{
    double az1, az2, interval;
    double transform[6];
    double newX, newY, rasValue, rasValue2, lowValue, reValue;
    bool found;
    int nSamples;
    az1 = Geometry::addDegrees(azimuth, 90.0);
    az2 = Geometry::addDegrees(azimuth, -90.0);

    //qDebug()<<"opening raster";
    GDALDataset *pRas;
    pRas = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pRas->GetGeoTransform(transform);
    GDALClose(pRas);
    //qDebug()<<"raster closed";

    interval = transform[1]/10.0;
    nSamples = ceil(distance/interval);

    //qDebug()<<"starting loop";
    found = false;
    for (int i=0; i<nSamples; i++)
    {
        Geometry::calcCoords(startX, startY, az1, interval*(i+1), newX, newY);
        rasValue = valueAtPoint(newX, newY);
        rasValue2 = rasterValueAtPoint(checkRasPath, newX, newY);
        if (rasValue2 > 0.5)
        {
            reValue = rasValue;
            found = true;
            i = nSamples;
            x = newX, y = newY;
        }
        if (i==0)
        {
            lowValue = rasValue;
            x = newX, y = newY;
        }
        else
        {
            if (rasValue < lowValue && rasValue != noData)
            {
                lowValue = rasValue;
                x = newX, y = newY;
            }
        }
        Geometry::calcCoords(startX, startY, az2, interval*(i+1), newX, newY);
        rasValue = valueAtPoint(newX, newY);
        rasValue2 = rasterValueAtPoint(checkRasPath, newX, newY);
        if (rasValue2 > 0.5 && rasValue != noData)
        {
            reValue = rasValue;
            found = true;
            i = nSamples;
            x = newX, y = newY;
        }
        if (rasValue < lowValue && rasValue != noData)
        {
            lowValue = rasValue;
            x = newX, y = newY;
        }
    }

    if (!found)
    {
        reValue = lowValue;
    }

    return reValue;
}

double Raster::sampleAlongLine_RasterVal(const char *rasterPath, const char *checkRasPath, double startX, double startY, double azimuth, double distance, double &x, double &y)
{
    setProperties(rasterPath);

    double value = sampleAlongLine_RasterVal(checkRasPath, startX, startY, azimuth, distance, x, y);

    return value;
}

void Raster::setNoData(double noDataValue, double minDataValue, double maxDataValue)
{
    GDALDataset *pSourceDS;
    pSourceDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_Update);
    pSourceDS->GetRasterBand(1)->SetNoDataValue(noDataValue);

    float *readRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *writeRow = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pSourceDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, readRow, nCols, 1, GDT_Float32, 0, 0);
        for (int j=0; j<nCols; j++)
        {
            if (readRow[j] < minDataValue || readRow[j] > maxDataValue)
            {
                writeRow[j] = noDataValue;
            }
            else
            {
                writeRow[j] = readRow[j];
            }
        }
        pSourceDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, writeRow, nCols, 1, GDT_Float32, 0, 0);
    }

    CPLFree(readRow);
    CPLFree(writeRow);

    GDALClose(pSourceDS);
}

void Raster::setNoData(const char *rasterPath, double noDataValue, double minDataValue, double maxDataValue)
{
    setProperties(rasterPath);
    setNoData(noDataValue, minDataValue, maxDataValue);
    setProperties(rasterPath);
//    for (int i=0; i<6; i++)
//    {
//        qDebug()<<transform[i];
//    }
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

    pSourceDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
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

    CPLFree(srcRow);
    CPLFree(subRow);
    CPLFree(newRow);

    GDALClose(pSourceDS);
    GDALClose(pSubtractDS);
    GDALClose(pOutDS);
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

void Raster::toXYZ(const char *rasterPath, const char *xyzPath)
{
    setProperties(rasterPath);

    GDALDataset *pRaster;
    pRaster = (GDALDataset*) GDALOpen(rasterPath, GA_ReadOnly);

    double xtlCenter, ytlCenter, xCenter, yCenter;

    xtlCenter = transform[0] + (transform[1]/2.0);
    ytlCenter = transform[3] + (fabs(transform[5])/2.0);

    float *read = (float*) CPLMalloc(sizeof(float)*1);

    QFile fout(xyzPath);
    fout.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&fout);
    out.setRealNumberNotation(QTextStream::FixedNotation);

    for (int i=0; i<nRows; i++)
    {
        yCenter = ytlCenter + (i*fabs(transform[5]));
        for (int j=0; j<nCols; j++)
        {
            pRaster->GetRasterBand(1)->RasterIO(GF_Read,j,i,1,1,read,1,1,GDT_Float32,0,0);
            if (*read != noData)
            {
                xCenter = xtlCenter + (j*transform[1]);
                out.setRealNumberPrecision(5);
                out << xCenter << "\t" <<yCenter << "\t" << *read <<"\n";
            }
        }
    }
    fout.close();
    CPLFree(read);
    GDALClose(pRaster);
}

void Raster::translateToGeoTIFF(const char *inPath, const char *outPath)
{
    setProperties(inPath);

    GDALDataset *inDS, *outDS;
    inDS = (GDALDataset*) GDALOpen(inPath, GA_ReadOnly);
    outDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);

    float *inval = (float*) CPLMalloc(sizeof(float)*nCols);
    float *outval = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=500; i<502; i++)
    {
        inDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, inval, nCols, 1, GDT_Float32, 0, 0);
        for (int j=0; j<nCols; j++)
        {
            outval[j] = inval[j];
            qDebug()<<inval[j];
        }
        outDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, outval, nCols, 1, GDT_Float32, 0, 0);
    }
}

double Raster::value(const char *rasterPath, int row, int col, int band)
{
    setProperties(rasterPath);
    double value;
    GDALDataset *pRaster;
    pRaster = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);

    float *rasVal = (float*) CPLMalloc(sizeof(float)*1);

    pRaster->GetRasterBand(band)->RasterIO(GF_Read, col, row, 1, 1, rasVal, 1, 1, GDT_Float32, 0, 0);

    value = *rasVal;

    GDALClose(pRaster);
    CPLFree(rasVal);

    return value;
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

    CPLFree(rasVal);

    GDALClose(pRaster);

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

double Raster::rasterValueAtPoint(const char *rasterPath, double xCoord, double yCoord)
{
    Raster raster;
    raster.setProperties(rasterPath);
    double value = raster.valueAtPoint(xCoord, yCoord);
    return value;
}

void Raster::loadDrivers()
{
    GDALAllRegister();

    pDriverTiff = GetGDALDriverManager()->GetDriverByName("GTiff");
}
