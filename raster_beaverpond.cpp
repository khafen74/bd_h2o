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

