#include "random.h"

Random::Random(double mean, double stdev)
{
    setMean(mean);
    setStdDev(stdev);
}

double Random::getMean()
{
    return m_mean;
}

double Random::getStdDev()
{
    return m_stdev;
}

void Random::setMean(double mean)
{
    m_mean = mean;
}

void Random::setStdDev(double stdev)
{
    m_stdev = stdev;
}

//generate random number from normal distribution using Box-Muller Transformation
double Random::getRandomNormal(double stdev, double mean)
{
    double x1, x2, w, y1, y2;

    do
    {
        x1 = 2.0 * getUniformRandom() - 1.0;
        x2 = 2.0 * getUniformRandom() - 1.0;
        w = x1 * x1 + x2 * x2;
    }
    while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);

    y1 = x1 * w;
    y2 = x2 * w;

    return y1 * stdev + mean;
}

double Random::getRandomLogNormal(double stdev, double mean)
{
    double rn = getRandomNormal();
    double lnorm = exp(mean+stdev*rn);

    return lnorm;
}

//generate uniform random number between 0 and 1
double Random::getUniformRandom()
{
    return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}

