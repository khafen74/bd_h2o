#ifndef RASTER_H
#define RASTER_H

#include "gdal_priv.h"
#include <QtCore>

const double PI = 3.14159265;
const int ROW_OFFSET[9] = {-1,-1,-1,0,0,0,1,1,1};
const int COL_OFFSET[9] = {-1,0,1,-1,0,1,-1,0,1};
const int ROW_OFFSET5[25] = {-2,-2,-2,-2,-2,-1,-1,-1,-1,-1,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2};
const int COL_OFFSET5[25] = {-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2};

class Raster
{
public:
    Raster();

    void add(const char *addPath);
    void add(const char *sourcePath, const char *addPath);
    double area();
    double area(const char *sourcePath);
    void aspect(const char *aspectPath);
    void aspect(const char *sourcePath, const char *aspectPath);
    void demOfDifference(const char *oldDem, const char *newDem, const char *dodRaster);
    void filterLowPass(const char *filterRaster);
    void filterLowPass(const char *sourceRaster, const char *filterRaster);
    double findMax(const char *rasterPath);
    void fromXYZ(const char *rasterPath, const char *xyzPath, int cols, int rows, double noDataValue, double inTransform[], int headerRows = 0);
    int getCols();
    const char *getPath();
    int getRows();
    void greaterThan(const char *outPath, double value);
    void greaterThan(const char *inPath, const char *outPath, double value);
    void hillshade(const char *hlsdPath);
    void hillshade(const char *rasterPath, const char *hlsdPath);
    int regions(const char *regionsRaster);
    int regions(const char *inputRaster, const char *regionsRaster);
    void setProperties(const char *rasterPath);
    void slopeTOF(const char *slopePath);
    void slopeTOF(const char *sourcePath, const char *slopePath);
    void slopeDeg(const char *slopePath);
    void slopeDeg(const char *sourcePath, const char *slopePath);
    void subtract(const char *subtractPath);
    void subtract(const char *sourcePath, const char *subtractPath);
    void subtract(const char *sourcePath, const char *subtractPath, const char *outputPath);
    double sum();
    double sum(const char *rasterPath);
    void zeroToNoData(const char *sourcePath, double noDataValue);
protected:

private:
    int nRows, nCols;
    double  transform[6], noData;
    const char *m_rasterPath;
    GDALDriver *pDriverTiff;

    void loadDrivers();
};

#endif // RASTER_H
