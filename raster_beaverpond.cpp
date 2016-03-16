#include "raster_beaverpond.h"

Raster_BeaverPond::Raster_BeaverPond()
{

}

void Raster_BeaverPond::createHANDInput(const char *pondPath, const char *facPath, const char *outPath)
{
    setProperties(pondPath);
    GDALDataset *pPondDS, *pFacDS, *pOutDS;

    pPondDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_Update);
    pFacDS = (GDALDataset*) GDALOpen(facPath, GA_ReadOnly);
    pOutDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);
    pOutDS->SetGeoTransform(transform);
    pOutDS->GetRasterBand(1)->SetNoDataValue(noData);
    pOutDS->GetRasterBand(1)->Fill(noData);
    float *outRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *facRow = (float *) CPLMalloc(sizeof(float)*nCols);
    float *pndRow = (float *) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pFacDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, facRow, nCols, 1, GDT_Float32, 0, 0);
        pPondDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, pndRow, nCols, 1, GDT_Float32, 0, 0);
        for (int j=0; j<nCols; j++)
        {
            if (pndRow[j] != noData)
            {
                outRow[j] = pndRow[j];
            }
            else if (facRow[j] > 0 && pndRow[j] == noData)
            {
                outRow[j] = -1.0;
            }
            else
            {
                outRow[j] = noData;
            }
        }
        pOutDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, outRow, nCols, 1, GDT_Float32, 0, 0);
    }

    CPLFree(outRow);
    CPLFree(facRow);
    CPLFree(pndRow);

    GDALClose(pPondDS);
    GDALClose(pFacDS);
    GDALClose(pOutDS);
}

void Raster_BeaverPond::heightAboveNetwork(const char *fdirPath, const char *facPath, const char *outPath, const char *outPondID)
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
                indices.clear();
                pFdirDS->GetRasterBand(1)->RasterIO(GF_Read, newCol-1, newRow-1, 3, 3, fdirWin, 3, 3, GDT_Byte, 0, 0);
                pFacDS->GetRasterBand(1)->RasterIO(GF_Read, newCol-1, newRow-1, 3, 3, facWin, 3, 3, GDT_Int32, 0, 0);

                if (facWin[4] > -1)
                {
                    write = true;
                    nIndex = 4;
                }
                else if (fdirWin[4] > 0)
                {
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
                        if (facWin[nIndex] > -1)
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
                    if (*elevVal < 0 || nIndex == 4)
                    {
                        *elevVal = noData;
                    }
                    pOutDS->GetRasterBand(1)->RasterIO(GF_Write, startCol, startRow, 1, 1, elevVal, 1, 1, GDT_Float32, 0, 0);
                    done = true;
                }
                nCount++;
            }
        }
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

void Raster_BeaverPond::heightAboveNetwork(const char *demPath, const char *fdirPath, const char *facPath, const char *outPath, const char *outPondID)
{
    setProperties(demPath);
    heightAboveNetwork(fdirPath, facPath, outPath, outPondID);
}

void Raster_BeaverPond::subtractHAND(const char *endPath, const char *outPath)
{
    GDALDataset *pSourceDS, *pEndDs, *pOutDS;
    pSourceDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pEndDs = (GDALDataset*) GDALOpen(endPath, GA_ReadOnly);
    pOutDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);
    pOutDS->SetGeoTransform(transform);
    pOutDS->GetRasterBand(1)->SetNoDataValue(noData);
    pOutDS->GetRasterBand(1)->Fill(noData);

    float *srcRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *endRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *outRow = (float*) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pSourceDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, srcRow, nCols, 1, GDT_Float32, 0, 0);
        pEndDs->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, endRow, nCols, 1, GDT_Float32, 0, 0);

        for (int j=0; j<nCols; j++)
        {
            if (srcRow[j] == noData || endRow[j] == noData)
            {
                outRow[j] = noData;
            }
            else
            {
                outRow[j] = srcRow[j] - endRow[j];
                if (outRow[j] <= 0)
                {
                    outRow[j] = noData;
                }
            }
        }

        pOutDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, outRow, nCols, 1, GDT_Float32, 0, 0);
    }

    CPLFree(srcRow);
    CPLFree(endRow);
    CPLFree(outRow);

    GDALClose(pSourceDS);
    GDALClose(pEndDs);
    GDALClose(pOutDS);
}

void Raster_BeaverPond::subtractHAND(const char *startPath, const char *endPath, const char *outPath)
{
    setProperties(startPath);
    subtractHAND(endPath, outPath);
}

