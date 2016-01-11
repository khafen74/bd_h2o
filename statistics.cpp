#include "statistics.h"

Statistics::Statistics(QVector<double> sample, DIST_TYPE distr)
{
    setSample(sample);
    setDistributionType(distr);
    calcMu();
    calcSigma();
}

void Statistics::calcConfidenceInterval(ConfInt ci)
{
    double lower, upper;
    setZScore(ci);

    switch(m_distr)
    {
        case RDT_norm:
            lower = m_mu - (m_z*m_se);
            upper = m_mu + (m_z*m_se);
            m_ciLo = lower;
            m_ciHi = upper;
        case RDT_lnorm:
            lower = m_mu+(m_sigma/2)-(1.96*sqrt((m_sigma/m_qvSample.length())+(pow(m_sigma, 2)/(2*(m_qvSample.length()-1)))));
            upper = m_mu+(m_sigma/2)+(1.96*sqrt((m_sigma/m_qvSample.length())+(pow(m_sigma, 2)/(2*(m_qvSample.length()-1)))));
            m_ciLo = exp(m_mu - (m_z*m_sigma));
            m_ciHi = exp(m_mu + (m_z*m_sigma));
    }
}

double Statistics::calcMeanLognormal()
{
    double sum = 0.0;

    for (int i=0; i<m_qvSample.length(); i++)
    {
        sum += log(m_qvSample[i]);
    }

    return (sum/m_qvSample.length());
}

double Statistics::calcMeanNormal()
{
    double sum = 0.0;

    for (int i=0; i<m_qvSample.length(); i++)
    {
        sum += m_qvSample[i];
    }

    return (sum/m_qvSample.length());
}

void Statistics::calcMu()
{
    switch(m_distr)
    {
        case RDT_norm: setMu(calcMeanNormal());
        case RDT_lnorm: setMu(calcMeanLognormal());
    }
}

void Statistics::calcSigma()
{
    switch(m_distr)
    {
        case RDT_norm: setSigma(calcStdNormal());
        case RDT_lnorm: setSigma(calcStdLognormal());
    }
}

double Statistics::calcStdLognormal()
{
    double logSum = 0.0;
    double logSumSq = 0.0;
    double var;

    for (int i=0; i<m_qvSample.length(); i++)
    {
        logSum += log(m_qvSample[i]);
        logSumSq += pow(log(m_qvSample[i]), 2);
    }

    var = sqrt(((m_qvSample.length()*logSumSq)-(pow(logSum, 2)))/(m_qvSample.length()*(m_qvSample.length()-1)));

    setStdError(sqrt(var/(m_qvSample.length())));

    return var;
}

double Statistics::calcStdNormal()
{

}

double Statistics::getLowerConfidenceLevel()
{
    return m_ciLo;
}

double Statistics::getStdError()
{
    return m_se;
}

double Statistics::getUpperConfidenceLevel()
{
    return m_ciHi;
}

void Statistics::setDistributionType(DIST_TYPE distr)
{
    m_distr = distr;
}

void Statistics::setMu(double mu)
{
    m_mu = mu;
}

void Statistics::setSample(QVector<double> sample)
{
    m_qvSample = sample;
}

void Statistics::setSigma(double sigma)
{
    m_sigma = sigma;
}

void Statistics::setStdError(double se)
{
    m_se = se;
}

void Statistics::setZScore(ConfInt ci)
{
    m_z = getZScore(ci);
}

double Statistics::getZScore(ConfInt ci)
{
    double zval = 0.0;

    switch (ci)
    {
        case CI_99: zval = 2.58;
        case CI_95: zval = 1.96;
        case CI_90: zval = 1.645;
        case CI_80: zval = 1.28;
        case CI_70: zval = 1.04;
    }

    return zval;
}
