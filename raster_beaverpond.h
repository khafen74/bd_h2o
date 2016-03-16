#ifndef RASTER_BEAVERPOND_H
#define RASTER_BEAVERPOND_H

#include "raster.h"


class Raster_BeaverPond : public Raster
{
public:
    Raster_BeaverPond();

    void createHANDInput(const char *pondPath, const char *facPath, const char *outPath);
    void heightAboveNetwork(const char *fdirPath, const char *facPath, const char *outPath, const char *outPondID);
    void heightAboveNetwork(const char *demPath, const char *fdirPath, const char *facPath, const char *outPath, const char *outPondID);
    void subtractHAND(const char *endPath, const char *outPath);
    void subtractHAND(const char *startPath, const char *endPath, const char *outPath);
    //void waterSurfaceElevation(const char *sourcePath, const char *waterDepth, const char *outPath);
};

#endif // RASTER_BEAVERPOND_H
