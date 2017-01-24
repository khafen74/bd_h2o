#ifndef VECTOROPS_H
#define VECTOROPS_H

#include <QtCore>

class VectorOps
{
public:
    VectorOps();

    static double max(QVector<double> data);
    static double sum(QVector<double> data);
};

#endif // VECTOROPS_H
