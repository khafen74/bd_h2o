#ifndef STORAGEMODEL_H
#define STORAGEMODEL_H

#include "dampolygons.h"
#include "reachlines.h"

class StorageModel
{
public:
    StorageModel(const char *bratPath, const char *outPath, const char *demPath, double capacity);

    void init(const char *bratPath, const char *outPath, const char *demPath, double capacity);

    void cleanOutDir();
    void run();
    void runCompare(const char *damsIn);

private:
    const char *m_bratPath, *m_outPath, *m_demPath;
    double bratCap;
};

#endif // STORAGEMODEL_H
