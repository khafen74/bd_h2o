#include "vectorops.h"

VectorOps::VectorOps()
{

}

double VectorOps::max(QVector<double> data)
{
    double max;

    max = data[0];

    for (int i=0 ; i<data.length(); i++)
    {
        if (data[i] > max)
        {
            max = data[i];
        }
    }

    return max;
}

