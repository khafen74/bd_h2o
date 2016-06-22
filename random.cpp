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
double Random::random_normal(double mean, double stdev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stdev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stdev + mean;
    }
}
//{
//    double x1, x2, w, d, y1, y2;

//    do
//    {
//        x1 = 2.0 * random_uniform() - 1.0;
//        x2 = 2.0 * random_uniform() - 1.0;
//        w = x1 * x1 + x2 * x2;
//    }
//    while (w >= 1.0);

//    d = sqrt((-2.0 * log(w)) / w);

//    y1 = x1 * d;
//    y2 = x2 * d;

//    return ((y1 * stdev) + mean);
//}

double Random::random_lognormal(double mean, double stdev)
{
    double rn = random_normal();
    double lnorm = exp(mean+stdev*rn);

    return lnorm;
}

QVector<double> Random::randomSeries(int count, DIST_TYPE distr, double mean, double stdev)
{
    QVector<double> sample;

    switch (distr)
    {
        case RDT_norm:
            for (int i=0; i<count; i++)
            {
                sample.append(random_normal(mean, stdev));
            }
            break;

        case RDT_lnorm:
            for (int i=0; i<count; i++)
            {
                sample.append(random_lognormal(mean, stdev));
            }
            break;
    }

    return sample;
}

//generate uniform random number between 0 and 1
double Random::random_uniform()
{
    return 2.0*rand()/RAND_MAX - 1;
}

