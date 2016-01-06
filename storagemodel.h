#ifndef STORAGEMODEL_H
#define STORAGEMODEL_H

#include "dampolygons.h"

class StorageModel
{
public:
    StorageModel();

    int MonteCarloRun_Pond();
    int MonteCarloRun_Reach();
    int SingleRun_Pond();
    int SingleRun_Reach();

private:
    int nIterations;
};

#endif // STORAGEMODEL_H
