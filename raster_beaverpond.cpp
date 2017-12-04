#include "raster_beaverpond.h"

Raster_BeaverPond::Raster_BeaverPond()
{

}

//Work the HAND algorithm backward to determine pond depth
void Raster_BeaverPond::backwardHAND(GDALDataset *flowDir, GDALDataset *dem, GDALDataset *idOut, GDALDataset *out, int startX, int startY, double startE, float *pondID, int maxCount, int &count)
{
    //Represents maximum dam height
    double maxHeight = 4.0;

    if (startX > 0 && startY > 0 && m_BHiter < m_maxBHiter)
    {
        for (int i=0; i<9; i++)
        {

            m_BHiter++;
            unsigned char *fdirWin = (unsigned char*) CPLMalloc(sizeof(unsigned char)*9);
            float *demWin = (float*) CPLMalloc(sizeof(float)*9);
            float *htAbove = (float*) CPLMalloc(sizeof(float)*1);
            float *htOld = (float*) CPLMalloc(sizeof(float)*1);

            flowDir->GetRasterBand(1)->RasterIO(GF_Read, startX-1, startY-1, 3, 3, fdirWin, 3, 3, GDT_Byte, 0, 0);
            dem->GetRasterBand(1)->RasterIO(GF_Read, startX-1, startY-1, 3, 3, demWin, 3, 3, GDT_Float32, 0, 0);

            int newX = startX;
            int newY = startY;
            *htAbove = demWin[i] - startE;

            //Index cell must drain to target cell
            if (drainsToMe(i, fdirWin[i]) && *htAbove < maxHeight && *htAbove > -10.0 && count < maxCount)
            {
                newX += COL_OFFSET[i];
                newY += ROW_OFFSET[i];
                out->GetRasterBand(1)->RasterIO(GF_Read, newX, newY, 1, 1, htOld, 1, 1, GDT_Float32, 0, 0);

                //Only write new value if it is less than previous value (deepest pond always wins)
                if (*htOld >= *htAbove || *htOld <= 0.0)
                {
                    idOut->GetRasterBand(1)->RasterIO(GF_Write, newX, newY, 1, 1, pondID, 1, 1, GDT_Float32, 0, 0);
                    out->GetRasterBand(1)->RasterIO(GF_Write, newX, newY, 1, 1, htAbove, 1, 1, GDT_Float32, 0, 0);
                    count++;
                }
                //Free memory before recursive call (should prevent allocating too much memory)
                CPLFree(htOld);
                CPLFree(demWin);
                CPLFree(fdirWin);
                CPLFree(htAbove);

                //If pond is less than max size run algorithm again

                backwardHAND(flowDir, dem, idOut, out, newX, newY, startE, pondID, maxCount, count);
            }
            else
            {
                CPLFree(htOld);
                CPLFree(demWin);
                CPLFree(fdirWin);
                CPLFree(htAbove);
            }
        }
    }
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
    pOutDS->SetProjection(pDemDS->GetProjectionRef());
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
    pOutDS->SetProjection(pDemDS->GetProjectionRef());
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
            if (wVal[j] > 0.0 || fVal[j] >= 0.0)
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

//Determine pond depth with the HAND algorithm, but works backward from dam cells and uses recursive function 'backwardHAND'
void Raster_BeaverPond::pondDepth_backwardHAND(const char *demPath, const char *fdirPath, const char *idPath, const char *outHtPath, const char *outIdPath)
{
    //the maximum allowable pond area in m^2
    double maxPondArea = 20000.0;
    setProperties(demPath);

    GDALDataset *pDem, *pFdir, *pId, *pHeight, *pIdOut;
    pDem = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_ReadOnly);
    pFdir = (GDALDataset*) GDALOpen(fdirPath, GA_ReadOnly);
    pId = (GDALDataset*) GDALOpen(idPath, GA_Update);

    pIdOut = pDriverTiff->Create(outIdPath, nCols, nRows, 1, GDT_Float32, NULL);
    pIdOut->SetGeoTransform(transform);
    pIdOut->GetRasterBand(1)->SetNoDataValue(noData);
    pIdOut->GetRasterBand(1)->Fill(noData);

    qDebug()<<"creating ht above";
    pHeight = pDriverTiff->Create(outHtPath, nCols, nRows, 1, GDT_Float32, NULL);
    pHeight->SetGeoTransform(transform);
    pHeight->GetRasterBand(1)->SetNoDataValue(noData);
    pHeight->GetRasterBand(1)->Fill(noData);

    float* idVal = (float*) CPLMalloc(sizeof(float)*1);
    float* demVal = (float*) CPLMalloc(sizeof(float)*1);
    unsigned char *fdirWin = (unsigned char*) CPLMalloc(sizeof(unsigned char)*9);

    int maxCount = ceil(maxPondArea/abs(transform[1]*transform[5]));
    m_maxBHiter = ceil(1000000.0/abs(transform[1]*transform[5]));
    qDebug()<<"max cell count per pond"<<maxCount;

    for (int i=1; i<nRows-1; i++)
    {
        for (int j=1; j<nCols-1; j++)
        {
            pId->GetRasterBand(1)->RasterIO(GF_Read, j, i, 1, 1, idVal, 1, 1, GDT_Float32, 0, 0);

            if (*idVal >= 0.0)
            {
                pDem->GetRasterBand(1)->RasterIO(GF_Read, j, i, 1, 1, demVal, 1, 1, GDT_Float32, 0, 0);
                pFdir->GetRasterBand(1)->RasterIO(GF_Read, j-1, i-1, 3, 3, fdirWin, 3, 3, GDT_Byte, 0, 0);

                for (int l=0; l<9; l++)
                {
                    int count=0;
                    m_BHiter = 0;
                    int maxCount6 = ceil(maxCount/6.0);
                    double startE = *demVal;

                    //Meat of the algorithm, a recursive function
                    backwardHAND(pFdir, pDem, pIdOut, pHeight, j, i, startE, idVal, maxCount6, count);
                }
            }
        }
        //qDebug()<<"row"<<i<<"done of"<<nRows;
    }

    CPLFree(idVal);
    CPLFree(demVal);
    CPLFree(fdirWin);

    GDALClose(pDem);
    GDALClose(pFdir);
    GDALClose(pId);
    GDALClose(pIdOut);
    GDALClose(pHeight);
    qDebug()<<"all datasets closed";
}

void Raster_BeaverPond::soilRasterCreation(const char *demPath, const char *huc8Path, const char *huc12Path, const char *inSoil, const char *outSoil, double maxVal)
{
    setProperties(demPath);
    GDALDataset *pDemDs, *pOutDs;
    pDemDs = (GDALDataset*) GDALOpen(demPath, GA_ReadOnly);
    pOutDs = pDriverTiff->Create(outSoil, nCols, nRows, 1, GDT_Float32, NULL);
    pOutDs->SetGeoTransform(transform);
    pOutDs->GetRasterBand(1)->SetNoDataValue(noData);
    pOutDs->GetRasterBand(1)->Fill(noData);

    float *demRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *outRow = (float*) CPLMalloc(sizeof(float)*nCols);

    double x, y, val8, val12, soilVal;
    int index8, index12;
    QVector<double> id8, sum8, id12, sum12;
    QVector<int> count8, count12;

    for (int i=0; i<nRows; i++)
    {
        pDemDs->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, demRow, nCols, 1, GDT_Float32, 0, 0);
        //pInDs->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, soilRow, nCols, 1, GDT_Float32, 0, 0);
        y = yCoordinate(i);

        for (int j=0; j<nCols; j++)
        {
            if (demRow[j] > 0.0)
            {
                x = xCoordinate(j);
                soilVal = rasterValueAtPoint(inSoil, x, y);

                if (soilVal > 0.0 && soilVal < maxVal)
                {
                    outRow[j] = soilVal;
                    val8 = rasterValueAtPoint(huc8Path, x, y);
                    val12 = rasterValueAtPoint(huc12Path, x, y);
                    index8 = id8.indexOf(val8);
                    index12 = id12.indexOf(val12);

                    if (index8 > -1)
                    {
                        sum8[index8] += soilVal;
                        count8[index8] += 1;
                    }
                    else if (val8 > 0.0)
                    {
                        id8.append(val8);
                        sum8.append(soilVal);
                        count8.append(1);
                        qDebug()<<"HUC 8 added"<<val8<<index8;
                    }
                    else
                    {
                        qDebug()<<"HUC 8 sampled outside of extent";
                    }

                    if (index12 > -1)
                    {
                        sum12[index12] += soilVal;
                        count12[index12] += 1;
                    }
                    else if (val12 > 0.0)
                    {
                        id12.append(val12);
                        sum12.append(soilVal);
                        count12.append(1);
                        qDebug()<<"HUC 12 added"<<val12<<index12;
                    }
                    else
                    {
                        qDebug()<<"HUC 12 sampled outside of extent";
                    }
                }
                else
                {
                    outRow[j] = noData;
                }
            }
            else
            {
                outRow[j] = noData;
            }
        }

        pOutDs->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, outRow, nCols, 1, GDT_Float32, 0, 0);
        if (i % 1000 == 0)
        {
            qDebug()<<"row"<<i+1<<"of"<<nRows;
        }
    }

    for (int i=0; i<nRows; i++)
    {
        pDemDs->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, demRow, nCols, 1, GDT_Float32, 0, 0);
        pOutDs->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, outRow, nCols, 1, GDT_Float32, 0, 0);
        y = yCoordinate(i);

        for (int j=0; j<nCols; j++)
        {
            x = xCoordinate(j);

            if (demRow[j] > 0.0 && outRow[j] <= 0.0000000001)
            {
                val12 = rasterValueAtPoint(huc12Path, x, y);
                index12 = id12.indexOf(val12);

                if(index12 > -1)
                {
                    outRow[j] = sum12[index12] / (count12[index12] * 1.0);
                    if (outRow[j] <= 0.0 || outRow[j] > maxVal)
                    {
                        qDebug()<<"Error, mean HUC12 value not correct"<<val12<<id12[index12]<<sum12[index12]<<count12[index12];
                    }
                }
                else
                {
                    val8 = rasterValueAtPoint(huc8Path, x, y);
                    index8 = id8.indexOf(val8);

                    if(index8 > -1)
                    {
                        outRow[j] = sum8[index8] / (count8[index8] * 1.0);
                        if (outRow[j] <= 0.0 || outRow[j] > maxVal)
                        {
                            qDebug()<<"Error, mean HUC12 value not correct"<<val8<<id8[index8]<<sum8[index8]<<count8[index8];
                        }
                    }
                    else
                    {
                        qDebug()<<"Error, no data value for soil"<<val12<<val8;
                    }
                }

            }
        }
        pOutDs->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, outRow, nCols, 1, GDT_Float32, 0, 0);
    }

    QFile fout("huc8.csv");
    fout.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&fout);

    for (int i=0; i<id8.length() ; i++)
    {
        out<<id8[i]<<","<<sum8[i]<<","<<count8[i]<<","<<sum8[i] / (count8[i] * 1.0)<<"\n";
        qDebug()<<id8[i]<<","<<sum8[i]<<","<<count8[i]<<","<<sum8[i] / (count8[i] * 1.0);
    }
    fout.close();
    QFile fout2("huc12.csv");
    fout2.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out2(&fout2);

    qDebug()<<"huc12s";
    for (int i=0; i<id12.length() ; i++)
    {
        out2<<id12[i]<<","<<sum12[i]<<","<<count12[i]<<","<<sum12[i] / (count12[i] * 1.0)<<"\n";
        qDebug()<<id12[i]<<","<<sum12[i]<<","<<count12[i]<<","<<sum12[i] / (count12[i] * 1.0);
    }
    fout2.close();
    CPLFree(demRow);
    CPLFree(outRow);

    GDALClose(pDemDs);
    GDALClose(pOutDs);
}

void Raster_BeaverPond::soilRasterCreation_table(const char *demPath, const char *huc8Path, const char *huc12Path, const char *soilPath, double maxVal)
{
    QVector<double> id8, mean8, id12, mean12;
    QFile fileh8 ("G:/01_etal/GIS_Data/USA/Soil/SSURGO/Utah/BearRiver/16010101_UpperBear/huc8_fc.csv");
    if (!fileh8.open(QIODevice::ReadOnly))
    {
        qDebug()<<"error opening file";
    }
    QTextStream streamh8(&fileh8);
    while (!fileh8.atEnd())
    {
        QString line = streamh8.readLine();
        id8.append(line.split(",").first().toDouble());
        mean8.append(line.split(",").last().toDouble());
        qDebug()<<id8.last()<<mean8.last();
    }
    fileh8.close();
    QFile fileh12 ("G:/01_etal/GIS_Data/USA/Soil/SSURGO/Utah/BearRiver/16010101_UpperBear/huc12_fc.csv");
    if (!fileh12.open(QIODevice::ReadOnly))
    {
        qDebug()<<"error opening file";
    }
    QTextStream streamh12 (&fileh12);
    while (!fileh12.atEnd())
    {
        QString line = streamh12.readLine();
        id12.append(line.split(",").first().toDouble());
        mean12.append(line.split(",").last().toDouble());
        qDebug()<<id12.last()<<mean12.last();
    }
    fileh12.close();
    setProperties(demPath);
    GDALDataset *pDemDs, *pOutDs;
    pDemDs = (GDALDataset*) GDALOpen(demPath, GA_ReadOnly);
    pOutDs = (GDALDataset*) GDALOpen(soilPath, GA_Update);
    float *demRow = (float*) CPLMalloc(sizeof(float)*nCols);
    float *outRow = (float*) CPLMalloc(sizeof(float)*nCols);

    double x, y, val8, val12;
    int index8, index12;

    for (int i=0; i<nRows; i++)
    {
        pDemDs->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, demRow, nCols, 1, GDT_Float32, 0, 0);
        pOutDs->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, outRow, nCols, 1, GDT_Float32, 0, 0);
        y = yCoordinate(i);

        for (int j=0; j<nCols; j++)
        {
            if (demRow[j] > 0.0 && outRow[j] <= 0.0000000001)
            {
                x = xCoordinate(j);
                val12 = rasterValueAtPoint(huc12Path, x, y);
                index12 = id12.indexOf(val12);

                if(index12 > -1)
                {
                    outRow[j] = mean12[index12];
                    if (outRow[j] <= 0.0 || outRow[j] > maxVal)
                    {
                        qDebug()<<"Error, mean HUC12 value not correct";
                    }
                }
                else
                {
                    val8 = rasterValueAtPoint(huc8Path, x, y);
                    index8 = id8.indexOf(val8);

                    if(index8 > -1)
                    {
                        outRow[j] = mean8[index8];
                        if (outRow[j] <= 0.0 || outRow[j] > maxVal)
                        {
                            qDebug()<<"Error, mean HUC12 value not correct";
                            outRow[j] = mean8[0];
                        }
                    }
                    else
                    {
                        qDebug()<<"Error, no data value for soil"<<val12<<val8;
                    }
                }

            }
        }
        pOutDs->GetRasterBand(1)->RasterIO(GF_Write, 0, i, nCols, 1, outRow, nCols, 1, GDT_Float32, 0, 0);
    }
    CPLFree(demRow);
    CPLFree(outRow);
    GDALClose(pDemDs);
    GDALClose(pOutDs);
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

