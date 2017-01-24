#include "vectorops.h"

VectorOps::VectorOps()
{

}

double VectorOps::max(QVector<double> data)
{
    double max;

    max = data[0];

    for (int i=0; i<data.length(); i++)
    {
        if (data[i] > max)
        {
            max = data[i];
        }
    }

    return max;
}

double VectorOps::sum(QVector<double> data)
{
    double sum = 0.0;

    for (int i=0; i<data.length(); i++)
    {
        sum += data[i];
    }

    return sum;
}

