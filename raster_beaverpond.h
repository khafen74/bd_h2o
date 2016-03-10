#ifndef RASTER_BEAVERPOND_H
#define RASTER_BEAVERPOND_H

#include "raster.h"


class Raster_BeaverPond : public Raster
{
public:
    Raster_BeaverPond();

    void createHANDInput(const char *pondPath, const char *facPath, const char *outPath);
    //void waterSurfaceElevation(const char *sourcePath, const char *waterDepth, const char *outPath);
};

#endif // RASTER_BEAVERPOND_H
