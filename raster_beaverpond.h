#ifndef RASTER_BEAVERPOND_H
#define RASTER_BEAVERPOND_H

#include "raster.h"


class Raster_BeaverPond : public Raster
{
public:
    Raster_BeaverPond();

    void backwardHAND(GDALDataset *flowDir, GDALDataset *dem, GDALDataset *idOut, GDALDataset *out, int startX, int startY, double startE, float *pondID);
    void groundwaterDepth(const char *startDepth, const char *newDepth, const char *outPath);
    void createHANDInput(const char *pondPath, const char *facPath, const char *outPath);
    void head(const char *demPath, const char *facPath, const char *outPath);
    void head(const char *demPath, const char *facPath, const char *wetPath, const char *outPath);
    void flowDownstream(const char *fdirPath, const char *facPath, const char *demPath, const char *pidPath, const char *gwePath, const char *outPath);
    void flowDownstream(const char *gwPath, const char *fdirPath, const char *facPath, const char *demPath, const char *pidPath, const char *gwePath, const char *outPath);
    void heightAboveNetwork(const char *fdirPath, const char *facPath, const char *outPath, const char *outPondID);
    void heightAboveNetwork(const char *demPath, const char *fdirPath, const char *facPath, const char *outPath, const char *outPondID);
    void heightAboveNetwork_ponds(const char *demPath, const char *fdirPath, const char *facPath
                                  , const char *heightPathLo, const char *heightPathMid, const char *heightPathHi
                                  , const char *outPath, const char *outPondID
                                  , const char *outHeightLo, const char *outHeightMid, const char *outHeightHi);
    void pondDepth_backwardHAND(const char *demPath, const char *fdirPath, const char *idPath, const char *outHtPath, const char *outIdPath);
    void soilRasterCreation(const char *demPath, const char *huc8Path, const char *huc12Path, const char *inSoil, const char *outSoil);
    void subtractHAND(const char *endPath, const char *outPath);
    void subtractHAND(const char *startPath, const char *endPath, const char *outPath);
    //void waterSurfaceElevation(const char *sourcePath, const char *waterDepth, const char *outPath);
};

#endif // RASTER_BEAVERPOND_H
