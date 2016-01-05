#ifndef STORAGEMODEL_H
#define STORAGEMODEL_H


class StorageModel
{
public:
    StorageModel();

    int MonteCarloRun_Pond();
    int MonteCarloRun_Reach();
    int SingleRun_Pond();
    int SingleRun_Reach();
};

#endif // STORAGEMODEL_H
