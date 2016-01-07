#include "statistics.h"

Statistics::Statistics(QVector<double> sample)
{
    setSample(sample);
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

    return var;
}

void Statistics::setSample(QVector<double> sample)
{
    m_qvSample = sample;
}
