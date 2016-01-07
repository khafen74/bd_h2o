#ifndef STATISTICS_H
#define STATISTICS_H

#include <QtCore>
#include <cmath>

class Statistics
{
public:
    Statistics(QVector<double> sample);
    double calcMeanLognormal();
    double calcMeanNormal();
    double calcStdLognormal();
    double calcStdNormal();
    void setSample(QVector<double> sample);

private:
    double m_mu, m_sigma;
    QVector<double> m_qvSample;
};

#endif // STATISTICS_H
