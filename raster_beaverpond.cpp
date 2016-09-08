#include "raster_beaverpond.h"

Raster_BeaverPond::Raster_BeaverPond()
{

}

void Raster_BeaverPond::groundwaterDepth(const char *startDepth, const char *newDepth, const char *outPath)
{
    setProperties(startDepth);
    GDALDataset *pStartDS, *pNewDS, *pOutDS;
    pStartDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pNewDS = (GDALDataset*) GDALOpen(newDepth, GA_ReadOnly);
    pOutDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);
    pOutDS->SetGeoTransform(transform);
    pOutDS->GetRasterBand(1)->SetNoDataValue(-9999);
    pOutDS->GetRasterBand(1)->Fill(-9999);
    float *startRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *newRow = (float *) CPLMalloc(sizeof(float)*nCols);
    float *outRow = (float *) CPLMalloc(sizeof(float)*nCols);

    for (int i=0; i<nRows; i++)
    {
        pStartDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, startRow, nCols, 1, GDT_Float32, 0, 0);
        pNewDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, newRow, nCols, 1, GDT_Float32, 0, 0);
        for (int j=0; j<nCols; j++)
        {
            if (startRow[j] < 0.0 && newRow[j] < 0.0)
            {
                outRow[j] = -9999;
            }
            else
            {
                if (newRow[j] < startRow[j] && newRow[j] > 0.0)
                {
                    outRow[j] = newRow[j];
                }
                else
                {
                    outRow[j] = startRow[j];
                }
            }
        }
        pOutDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, outRow, nCols, 1, GDT_Float32, 0, 0);
    }

    CPLFree(startRow);
    CPLFree(newRow);
    CPLFree(outRow);

    GDALClose(pStartDS);
    GDALClose(pNewDS);
    GDALClose(pOutDS);
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
            if (pndRow[j] != noData && facRow[j] > 0) //now only applies pond data to stream raster
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

void Raster_BeaverPond::head(const char *demPath, const char *facPath, const char *outPath)
{
    setProperties(demPath);

    GDALDataset *pDemDS, *pFacDS, *pOutDS;

    pDemDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pFacDS = (GDALDataset*) GDALOpen(facPath, GA_ReadOnly);
    pOutDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);
    pOutDS->SetGeoTransform(transform);
    pOutDS->GetRasterBand(1)->Fill(-9999);
    pOutDS->GetRasterBand(1)->SetNoDataValue(-9999);
    float *eVal = (float*) CPLMalloc(sizeof(float)*nCols);
    float *fVal = (float*) CPLMalloc(sizeof(float)*nCols);
    float *hVal = (float*) CPLMalloc(sizeof(float)*nCols);
    for (int i=0; i<nRows; i++)
    {
        pDemDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, eVal, nCols, 1, GDT_Float32, 0, 0);
        pFacDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, fVal, nCols, 1, GDT_Float32, 0, 0);
        for (int j=0; j<nCols; j++)
        {
            if (fVal[j] >= 0.0)
            {
                hVal[j] = eVal[j];
            }
            else
            {
                hVal[j] = -9999.0;
            }
        }
        pOutDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, hVal, nCols, 1, GDT_Float32, 0, 0);
    }

    CPLFree(eVal);
    CPLFree(fVal);
    CPLFree(hVal);

    GDALClose(pDemDS);
    GDALClose(pFacDS);
    GDALClose(pOutDS);
}

void Raster_BeaverPond::head(const char *demPath, const char *facPath, const char *wetPath, const char *outPath)
{
    setProperties(demPath);

    GDALDataset *pDemDS, *pFacDS,*pWetDS, *pOutDS;
    //qDebug()<<demPath<<facPath<<wetPath<<outPath;
    pDemDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pFacDS = (GDALDataset*) GDALOpen(facPath, GA_ReadOnly);
    pWetDS = (GDALDataset*) GDALOpen(wetPath, GA_ReadOnly);
    pOutDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);
    pOutDS->SetGeoTransform(transform);
    pOutDS->GetRasterBand(1)->Fill(-9999);
    pOutDS->GetRasterBand(1)->SetNoDataValue(-9999);
    float *eVal = (float*) CPLMalloc(sizeof(float)*nCols);
    float *fVal = (float*) CPLMalloc(sizeof(float)*nCols);
    float *wVal = (float*) CPLMalloc(sizeof(float)*nCols);
    float *hVal = (float*) CPLMalloc(sizeof(float)*nCols);
    //qDebug()<<"starting loop";
    for (int i=0; i<nRows; i++)
    {
        //qDebug()<<"reading row"<<i;
        pDemDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, eVal, nCols, 1, GDT_Float32, 0, 0);
        //qDebug()<<"dem";
        pFacDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, fVal, nCols, 1, GDT_Float32, 0, 0);
        //qDebug()<<"fac";
        pWetDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, wVal, nCols, 1, GDT_Float32, 0, 0);
        //qDebug()<<"wet";
        for (int j=0; j<nCols; j++)
        {
            //qDebug()<<wVal[j]<<fVal[j]<<hVal[j];
            if (wVal[j] >= 0.0 || fVal[j] >= 0.0)
            {
                hVal[j] = eVal[j];
            }
            else
            {
                hVal[j] = -9999.0;
            }
        }
        pOutDS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, hVal, nCols, 1, GDT_Float32, 0, 0);
    }

    CPLFree(eVal);
    CPLFree(fVal);
    CPLFree(wVal);
    CPLFree(hVal);

    GDALClose(pDemDS);
    GDALClose(pFacDS);
    GDALClose(pWetDS);
    GDALClose(pOutDS);
}

void Raster_BeaverPond::flowDownstream(const char *fdirPath, const char *facPath, const char *demPath, const char *pidPath, const char *gwePath, const char *outPath)
{
    GDALDataset *pGwDS, *pFdirDS, *pFacDS, *pGweDS, *pOutDS;

    pGwDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_Update);
    pFdirDS = (GDALDataset*) GDALOpen(fdirPath, GA_ReadOnly);
    pFacDS = (GDALDataset*) GDALOpen(facPath, GA_ReadOnly);
    pGweDS = (GDALDataset*) GDALOpen(gwePath, GA_Update);
    pOutDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);
    pOutDS->SetGeoTransform(transform);
    pOutDS->GetRasterBand(1)->Fill(-9999);
    pOutDS->GetRasterBand(1)->SetNoDataValue(-9999);

    unsigned char *fdirVal = (unsigned char*) CPLMalloc(sizeof(unsigned char));
    signed long int *facVal = (signed long int*) CPLMalloc(sizeof(signed long int));
    float *gwVal = (float*) CPLMalloc(sizeof(float)*1);
    float *gweVal = (float*) CPLMalloc(sizeof(float)*1);

    double gweStart;
    int newRow, newCol, nIndex;
    int nChanged = 0;

    for (int i=0; i<nRows; i++)
    {
        for (int j=0; j<nCols; j++)
        {
            pGwDS->GetRasterBand(1)->RasterIO(GF_Read, j, i, 1, 1, gwVal, 1, 1, GDT_Float32, 0, 0);

            if (*gwVal != noData && *gwVal > 0.0)
            {
                newRow = i, newCol = j;
                pGweDS->GetRasterBand(1)->RasterIO(GF_Read, j, i, 1, 1, gweVal, 1, 1, GDT_Float32, 0, 0);
                gweStart = *gweVal;
                pFdirDS->GetRasterBand(1)->RasterIO(GF_Read, j, i, 1, 1, fdirVal, 1, 1, GDT_Byte, 0, 0);
                nIndex = getD8Index(*fdirVal);
                newRow += ROW_OFFSET[nIndex];
                newCol += COL_OFFSET[nIndex];
                pGwDS->GetRasterBand(1)->RasterIO(GF_Read, newCol, newRow, 1, 1, gwVal, 1, 1, GDT_Float32, 0, 0);
                pFacDS->GetRasterBand(1)->RasterIO(GF_Read, newCol, newRow, 1, 1, facVal, 1, 1, GDT_Int32, 0, 0);
                int nCount = 0;

                while (nIndex > -1 && *facVal < 1 && *gwVal <= 0.0 && nCount < 500)
                {

                    pFdirDS->GetRasterBand(1)->RasterIO(GF_Read, newCol, newRow, 1, 1, fdirVal, 1, 1, GDT_Byte, 0, 0);
                    if (*gweVal > gweStart)
                    {
                        *gweVal = gweStart;
                    }
                    pGweDS->GetRasterBand(1)->RasterIO(GF_Write, newCol, newRow, 1, 1, gweVal, 1, 1, GDT_Float32, 0, 0);
                    nIndex = getD8Index(*fdirVal);
                    newRow += ROW_OFFSET[nIndex];
                    newCol += COL_OFFSET[nIndex];
                    pGwDS->GetRasterBand(1)->RasterIO(GF_Read, newCol, newRow, 1, 1, gwVal, 1, 1, GDT_Float32, 0, 0);
                    pFacDS->GetRasterBand(1)->RasterIO(GF_Read, newCol, newRow, 1, 1, facVal, 1, 1, GDT_Int32, 0, 0);
                    nCount++;
                    nChanged++;
                    if (nCount == 499)
                    {
                        qDebug()<<"tracing over 500 cells downstream";
                    }
                }
            }
        }
    }

    qDebug()<<nChanged<<" gw cells updated";

    CPLFree(fdirVal);
    CPLFree(facVal);
    CPLFree(gwVal);
    CPLFree(gweVal);

    GDALClose(pGwDS);
    GDALClose(pFdirDS);
    GDALClose(pFacDS);
    GDALClose(pGweDS);
    GDALClose(pOutDS);
}

void Raster_BeaverPond::flowDownstream(const char *gwPath, const char *fdirPath, const char *facPath, const char *demPath, const char *pidPath, const char *gwePath, const char *outPath)
{
    qDebug()<<"starting flow downstream";
    //qDebug()<<"file paths"<<gwPath<<fdirPath<<facPath<<demPath<<pidPath<<gwePath<<outPath;
    setProperties(gwPath);
    flowDownstream(fdirPath, facPath, demPath, pidPath, gwePath, outPath);
}

void Raster_BeaverPond::heightAboveNetwork(const char *fdirPath, const char *facPath, const char *outPath, const char *outPondID)
{
    GDALDataset *pDemDS, *pFdirDS, *pFacDS, *pOutDS, *pIdDS;

    pDemDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pFdirDS = (GDALDataset*) GDALOpen(fdirPath, GA_ReadOnly);
    pFacDS = (GDALDataset*) GDALOpen(facPath, GA_ReadOnly);
    pOutDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);
    pIdDS = pDriverTiff->Create(outPondID, nCols, nRows, 1, GDT_Int32, NULL);
    pOutDS->SetGeoTransform(transform);
    pOutDS->GetRasterBand(1)->SetNoDataValue(noData);
    pOutDS->GetRasterBand(1)->Fill(noData);
    pIdDS->SetGeoTransform(transform);
    pIdDS->GetRasterBand(1)->SetNoDataValue(noData);
    pIdDS->GetRasterBand(1)->Fill(noData);

    int nIndex, startRow, startCol, newRow, newCol;
    QVector<QString> indices;
    bool done, write;
    unsigned char *fdirWin = (unsigned char*) CPLMalloc(sizeof(unsigned char)*9);
    signed long int *facWin = (signed long int*) CPLMalloc(sizeof(signed long int)*9);
    float *elevValStart = (float*) CPLMalloc(sizeof(float)*1);
    float *elevVal = (float*) CPLMalloc(sizeof(float)*1);
    signed long int *pondVal = (signed long int*) CPLMalloc(sizeof(signed long int));

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

                if (facWin[4] >= -1)
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
                    *pondVal = facWin[nIndex];
                    if (*elevVal < 0 || nIndex == 4)
                    {
                        *elevVal = noData;
                    }
                    pOutDS->GetRasterBand(1)->RasterIO(GF_Write, startCol, startRow, 1, 1, elevVal, 1, 1, GDT_Float32, 0, 0);
                    pIdDS->GetRasterBand(1)->RasterIO(GF_Write, startCol, startRow, 1, 1, pondVal, 1, 1, GDT_Int32, 0, 0);
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
    CPLFree(pondVal);

    GDALClose(pDemDS);
    GDALClose(pFdirDS);
    GDALClose(pFacDS);
    GDALClose(pOutDS);
    GDALClose(pIdDS);
}

void Raster_BeaverPond::heightAboveNetwork(const char *demPath, const char *fdirPath, const char *facPath, const char *outPath, const char *outPondID)
{
    setProperties(demPath);
    heightAboveNetwork(fdirPath, facPath, outPath, outPondID);
}

void Raster_BeaverPond::heightAboveNetwork_ponds(const char *demPath, const char *fdirPath, const char *facPath
                                                 , const char *heightPathLo, const char *heightPathMid, const char *heightPathHi
                                                 , const char *outPath, const char *outPondID
                                                 , const char *outHeightLo, const char *outHeightMid, const char *outHeightHi)
{
    setProperties(demPath);

    GDALDataset *pDemDS, *pFdirDS, *pFacDS, *pOutDS, *pIdDS
            , *pHtOutLoDS, *pHtOutMidDS, *pHtOutHiDS
            , *pHtInLoDS, *pHtInMidDS, *pHtInHiDS;

    pDemDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pFdirDS = (GDALDataset*) GDALOpen(fdirPath, GA_ReadOnly);
    pFacDS = (GDALDataset*) GDALOpen(facPath, GA_ReadOnly);

    pHtInLoDS = (GDALDataset*) GDALOpen(heightPathLo, GA_ReadOnly);
    pHtInMidDS = (GDALDataset*) GDALOpen(heightPathMid, GA_ReadOnly);
    pHtInHiDS = (GDALDataset*) GDALOpen(heightPathHi, GA_ReadOnly);

    pOutDS = pDriverTiff->Create(outPath, nCols, nRows, 1, GDT_Float32, NULL);
    pIdDS = pDriverTiff->Create(outPondID, nCols, nRows, 1, GDT_Float32, NULL);

    pHtOutLoDS = pDriverTiff->Create(outHeightLo, nCols, nRows, 1, GDT_Float32, NULL);
    pHtOutMidDS = pDriverTiff->Create(outHeightMid, nCols, nRows, 1, GDT_Float32, NULL);
    pHtOutHiDS = pDriverTiff->Create(outHeightHi, nCols, nRows, 1, GDT_Float32, NULL);

    pOutDS->SetGeoTransform(transform);
    pOutDS->GetRasterBand(1)->SetNoDataValue(noData);
    pOutDS->GetRasterBand(1)->Fill(noData);
    pIdDS->SetGeoTransform(transform);
    pIdDS->GetRasterBand(1)->SetNoDataValue(noData);
    pIdDS->GetRasterBand(1)->Fill(noData);
    pHtOutLoDS->SetGeoTransform(transform);
    pHtOutLoDS->GetRasterBand(1)->SetNoDataValue(noData);
    pHtOutLoDS->GetRasterBand(1)->Fill(noData);
    pHtOutMidDS->SetGeoTransform(transform);
    pHtOutMidDS->GetRasterBand(1)->SetNoDataValue(noData);
    pHtOutMidDS->GetRasterBand(1)->Fill(noData);
    pHtOutHiDS->SetGeoTransform(transform);
    pHtOutHiDS->GetRasterBand(1)->SetNoDataValue(noData);
    pHtOutHiDS->GetRasterBand(1)->Fill(noData);
    qDebug()<<"datasets all initialized";

    int nIndex, startRow, startCol, newRow, newCol;
    QVector<QString> indices;
    bool done, write;
    unsigned char *fdirWin = (unsigned char*) CPLMalloc(sizeof(unsigned char)*9);
    signed long int *facWin = (signed long int*) CPLMalloc(sizeof(signed long int)*9);
    float *elevValStart = (float*) CPLMalloc(sizeof(float)*1);
    float *elevVal = (float*) CPLMalloc(sizeof(float)*1);
    float *htValLo = (float*) CPLMalloc(sizeof(float)*1);
    float *htValMid = (float*) CPLMalloc(sizeof(float)*1);
    float *htValHi = (float*) CPLMalloc(sizeof(float)*1);
    float *pondVal = (float*) CPLMalloc(sizeof(float));

    //qDebug()<<"starting loop";
    for (int i=1; i<nRows-1; i++)
    {
        if ((i+1)%500 == 0)
        {
            qDebug()<<"row"<<i+1<<"of"<<nRows;
        }
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

                if (facWin[4] >= -1)
                {
                    //make true if you want to inundate the cell(s) marking dam location
                    write = false;
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
                    pHtInLoDS->GetRasterBand(1)->RasterIO(GF_Read, newCol, newRow, 1, 1, htValLo, 1, 1, GDT_Float32, 0, 0);
                    pHtInMidDS->GetRasterBand(1)->RasterIO(GF_Read, newCol, newRow, 1, 1, htValMid, 1, 1, GDT_Float32, 0, 0);
                    pHtInHiDS->GetRasterBand(1)->RasterIO(GF_Read, newCol, newRow, 1, 1, htValHi, 1, 1, GDT_Float32, 0, 0);
                    *elevVal = *elevValStart - *elevVal;
                    *pondVal = facWin[nIndex];
                    if (*elevVal < 0 || nIndex == 4)
                    {
                        *elevVal = noData;
                    }
                    pOutDS->GetRasterBand(1)->RasterIO(GF_Write, startCol, startRow, 1, 1, elevVal, 1, 1, GDT_Float32, 0, 0);
                    pIdDS->GetRasterBand(1)->RasterIO(GF_Write, startCol, startRow, 1, 1, pondVal, 1, 1, GDT_Float32, 0, 0);
                    pHtOutLoDS->GetRasterBand(1)->RasterIO(GF_Write, startCol, startRow, 1, 1, htValLo, 1, 1, GDT_Float32, 0, 0);
                    pHtOutMidDS->GetRasterBand(1)->RasterIO(GF_Write, startCol, startRow, 1, 1, htValMid, 1, 1, GDT_Float32, 0, 0);
                    pHtOutHiDS->GetRasterBand(1)->RasterIO(GF_Write, startCol, startRow, 1, 1, htValHi, 1, 1, GDT_Float32, 0, 0);
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
    CPLFree(pondVal);
    CPLFree(htValLo);
    CPLFree(htValMid);
    CPLFree(htValHi);

    GDALClose(pDemDS);
    GDALClose(pFdirDS);
    GDALClose(pFacDS);
    GDALClose(pOutDS);
    GDALClose(pIdDS);
    GDALClose(pHtInLoDS);
    GDALClose(pHtInMidDS);
    GDALClose(pHtInHiDS);
    GDALClose(pHtOutLoDS);
    GDALClose(pHtOutMidDS);
    GDALClose(pHtOutHiDS);
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

