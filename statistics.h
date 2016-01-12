#ifndef STATISTICS_H
#define STATISTICS_H

#include <QtCore>
#include <cmath>
#include "random.h"
#include "vectorops.h"

enum ConfInt
{
    CI_99
    ,CI_95
    ,CI_90
    ,CI_80
    ,CI_70
};

class Statistics
{
public:
    Statistics(QVector<double> sample, DIST_TYPE distr = RDT_norm);
    void calcConfidenceInterval(ConfInt ci);
    void calcCredibleInterval(ConfInt ci);
    double calcMeanLognormal();
    double calcMeanNormal();
    void calcMu();
    void calcSigma();
    double calcStdLognormal();
    double calcStdNormal();
    QVector<double> getData();
    double getLowerConfidenceLevel();
    double getStdError();
    double getQuantile(double quantile);
    double getUpperConfidenceLevel();
    void setDistributionType(DIST_TYPE distr);
    void setMu(double mu);
    void setSample(QVector<double> sample);
    void setSigma(double sigma);
    void setStdError(double se);
    void setZScore(ConfInt ci);

    static double getZScore(ConfInt ci);
    static double getCIDecimal(ConfInt ci);

private:
    double m_mu, m_sigma, m_se, m_z, m_ciLo, m_ciHi;
    QVector<double> m_qvSample;
    DIST_TYPE m_distr;
};

#endif // STATISTICS_H
