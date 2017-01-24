#ifndef RANDOM_H
#define RANDOM_H

#include <QtCore>
#include <cmath>
#include <cstdlib>
#include <ctime>

enum DIST_TYPE
{
    RDT_lnorm
    ,RDT_norm
};

class Random
{
public:
    Random(double mean = 1.0, double stdev = 1.0);

    double getMean();
    double getStdDev();
    void setMean(double mean);
    void setStdDev(double stdev);

    static double random_normal(double mean = 0.0, double stdev = 1.0);
    static double random_lognormal(double mean = 0.0, double stdev = 1.0);
    static QVector<double> randomSeries(int count = 1000, DIST_TYPE distr = RDT_norm, double mean = 0.0, double stdev = 1.0);
    static double random_uniform();

private:
    double m_mean;
    double m_stdev;
};

#endif // RANDOM_H
