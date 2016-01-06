#ifndef RANDOM_H
#define RANDOM_H

#include <cmath>
#include <cstdlib>

class Random
{
public:
    Random(double mean = 1.0, double stdev = 1.0);

    double getMean();
    double getStdDev();
    void setMean(double mean);
    void setStdDev(double stdev);

    static double random_normal(double stdev = 1.0, double mean = 0.0);
    static double random_lognormal(double stdev = 1.0, double mean = 0.0);
    static double random_uniform();

private:
    double m_mean;
    double m_stdev;
};

#endif // RANDOM_H
