#ifndef RASTER_H
#define RASTER_H

#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "ogr_core.h"
#include "ogr_api.h"
#include <QtCore>
#include "geometry.h"

const int ROW_OFFSET[9] = {-1,-1,-1,0,0,0,1,1,1};
const int COL_OFFSET[9] = {-1,0,1,-1,0,1,-1,0,1};
const int FLOW_DIR[9] = {32,64,128,16,0,1,8,4,2};
const int ROW_OFFSET5[25] = {-2,-2,-2,-2,-2,-1,-1,-1,-1,-1,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2};
const int COL_OFFSET5[25] = {-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2};

class Raster
{
public:
    Raster();

    void add(const char *addPath, const char *outPath);
    void add(const char *sourcePath, const char *addPath, const char *outPath);
    void addTo(const char *addPath);
    void addTo(const char *sourcePath, const char *addPath);
    double area();
    double area(const char *sourcePath);
    void aspect(const char *aspectPath);
    void aspect(const char *sourcePath, const char *aspectPath);
    int checkRowCol(int row, int col);
    void demOfDifference(const char *oldDem, const char *newDem, const char *dodRaster);
    void extractByMask_CellCenters(const char *rasterOut, const char *polygonPath, const char *lyrName);
    void extractByMask_CellCenters(const char *rasterPath, const char *rasterOut, const char *polygonPath, const char *lyrName);
    void filterLowPass(const char *filterRaster);
    void filterLowPass(const char *sourceRaster, const char *filterRaster);
    double findMax(const char *rasterPath);
    void fromXYZ(const char *rasterPath, const char *xyzPath, int cols, int rows, double noDataValue, double inTransform[], int headerRows = 0);
    int getCol(double xCoord);
    int getCols();
    int getD8Index(int nFdir);
    const char *getPath();
    double getRow(double yCoord);
    int getRows();
    void greaterThan(const char *outPath, double value);
    void greaterThan(const char *inPath, const char *outPath, double value);
    void heightAboveNetwork(const char *fdirPath, const char *facPath, const char *outPath);
    void heightAboveNetwork(const char *demPath, const char *fdirPath, const char *facPath, const char *outPath);
    void hillshade(const char *hlsdPath);
    void hillshade(const char *rasterPath, const char *hlsdPath);
    int regions(const char *regionsRaster);
    int regions(const char *inputRaster, const char *regionsRaster);
    double sampleAlongLine_LowVal(double startX, double startY, double azimuth, double distance, double &x, double &y);
    double sampleAlongLine_LowVal(const char * rasterPath, double startX, double startY, double azimuth, double distance, double &x, double &y);
    void setNoData(double noDataValue, double minDataValue, double maxDataValue);
    void setNoData(const char *rasterPath, double noDataValue, double minDataValue, double maxDataValue);
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
    double valueAtPoint(double xCoord, double yCoord);
    double valueAtPoint(const char *rasterPath, double xCoord, double yCoord);
    void writeCellValue(double xCoord, double yCoord, double value);
    void writeCellValue(const char *rasterPath, double xCoord, double yCoord, double value);
    double xCoordinate(int col);
    double yCoordinate(int row);
    void zeroToNoData(const char *sourcePath, double noDataValue);

protected:
    int nRows, nCols;
    double transform[6], noData;
    QString m_rasterPath;
    GDALDriver *pDriverTiff;

    void loadDrivers();
};

#endif // RASTER_H
