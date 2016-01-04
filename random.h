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

    static double getRandomNormal(double stdev = 1.0, double mean = 0.0);
    static double getRandomLogNormal(double stdev = 1.0, double mean = 0.0);
    static double getUniformRandom();

private:
    double m_mean;
    double m_stdev;
};

#endif // RANDOM_H
