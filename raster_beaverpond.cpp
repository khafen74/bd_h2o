#include "raster_beaverpond.h"

Raster_BeaverPond::Raster_BeaverPond()
{

}

void Raster_BeaverPond::createHANDInput(const char *pondPath, const char *facPath)
{
    setProperties(pondPath);
    GDALDataset *pPondDS, *pFacDS;

    pPondDS = (GDALDataset*) GDALOpen(m_rasterPath.toStdString().c_str(), GA_Update);
    pFacDS = (GDALDataset*) GDALOpen(facPath, GA_ReadOnly);
    float *pndVal = (float*) CPLMalloc(sizeof(float));
    float *facRow = (float *) CPLMalloc(sizeof(float)*nCols);
    *pndVal = -1.0;

    for (int i=0; i<nRows; i++)
    {
        pFacDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i, nCols, 1, facRow, nCols, 1, GDT_Float32, 0, 0);
        for (int j=0; j<nCols; j++)
        {
            if (facRow[j] > 0)
            {
                pPondDS->GetRasterBand(1)->RasterIO(GF_Write, j, i, 1, 1, pndVal, 1, 1, GDT_Float32, 0, 0);
            }
        }
    }

    CPLFree(pndVal);
    CPLFree(facRow);

    GDALClose(pPondDS);
    GDALClose(pFacDS);
}

