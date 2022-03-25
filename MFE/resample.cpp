#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>
#include "upfirdn.h"

#define m_pi 3.14159265358979323846
using namespace std;
#define OriginSize 4096
#define ResampleSize 1024

extern "C" double *resample(double *inputSignal, int upFactor, int downFactor, const char *path);

// 获得最大公约数
int getgcd(int num1, int num2)
{
    int tmp = 0;
    while (num1 > 0)
    {
        tmp = num1;
        num1 = num2 % num1;
        num2 = tmp;
    }
    return num2;
}

int quotientceil(int num1, int num2)
{
    if (num1 % num2 != 0)
        return num1 / num2 + 1;
    return num1 / num2;
}

double sinc(double x)
{
    if (fabs(x - 0.0) < 0.000001)
        return 1;
    return sin(m_pi * x) / (m_pi * x);
}

void firls(int length, vector<double> freq,
           const vector<double> &amplitude, vector<double> &result)
{

    vector<double> weight;
    int freqsize = freq.size();
    int weightsize = freqsize / 2;

    weight.reserve(weightsize);
    for (int i = 0; i < weightsize; i++)
        weight.push_back(1.0);

    int filterlength = length + 1;

    for (int i = 0; i < freqsize; i++)
        freq[i] /= 2.0;

    vector<double> dfreq;
    for (int i = 1; i < freqsize; i++)
        dfreq.push_back(freq[i] - freq[i - 1]);

    length = (filterlength - 1) / 2;
    int nodd = filterlength % 2;
    double b0 = 0.0;
    vector<double> k;
    if (nodd == 0)
    {
        for (int i = 0; i <= length; i++)
            k.push_back(i + 0.5);
    }
    else
    {
        for (int i = 0; i <= length; i++)
            k.push_back(i);
    }

    vector<double> b;
    int ksize = k.size();
    for (int i = 0; i < ksize; i++)
        b.push_back(0.0);
    for (int i = 0; i < freqsize; i += 2)
    {
        double slope = (amplitude[i + 1] - amplitude[i]) / (freq[i + 1] - freq[i]);
        double b1 = amplitude[i] - slope * freq[i];
        if (nodd == 1)
        {
            b0 += (b1 * (freq[i + 1] - freq[i])) +
                  slope / 2.0 * (freq[i + 1] * freq[i + 1] - freq[i] * freq[i]) *
                      fabs(weight[(i + 1) / 2] * weight[(i + 1) / 2]);
        }
        for (int j = 0; j < ksize; j++)
        {
            b[j] += (slope / (4 * m_pi * m_pi) *
                     (cos(2 * m_pi * k[j] * freq[i + 1]) - cos(2 * m_pi * k[j] * freq[i])) / (k[j] * k[j])) *
                    fabs(weight[(i + 1) / 2] * weight[(i + 1) / 2]);
            b[j] += (freq[i + 1] * (slope * freq[i + 1] + b1) * sinc(2 * k[j] * freq[i + 1]) -
                     freq[i] * (slope * freq[i] + b1) * sinc(2 * k[j] * freq[i])) *
                    fabs(weight[(i + 1) / 2] * weight[(i + 1) / 2]);
        }
    }
    if (nodd == 1)
        b[0] = b0;
    vector<double> a;
    double w0 = weight[0];
    for (int i = 0; i < ksize; i++)
        a.push_back((w0 * w0) * 4 * b[i]);
    if (nodd == 1)
    {
        a[0] /= 2;
        for (int i = length; i >= 1; i--)
            result.push_back(a[i] / 2.0);
        result.push_back(a[0]);
        for (int i = 1; i <= length; i++)
            result.push_back(a[i] / 2.0);
    }
    else
    {
        for (int i = length; i >= 0; i--)
            result.push_back(a[i] / 2.0);
        for (int i = 0; i <= length; i++)
            result.push_back(a[i] / 2.0);
    }
}

int readwindow(vector<double> &window, const char *path)
{
    ifstream windowstream(path);
    if (!windowstream.is_open())
    {
        return -1;
    }
    while (windowstream)
    {
        double data;
        windowstream >> data;
        window.push_back(data);
    }
    windowstream.close();
    return 0;
}

// outframe = resample(xFrame,            1024,          4096,             path);
double *resample(double *inputdata, int upfactor, int downfactor, const char *path)
{
    double *outputdata = (double *)malloc(ResampleSize * sizeof(double));
    vector<double> inputsignal(OriginSize);
    vector<double> outputsignal;
    for (int i = 0; i < OriginSize; i++)
    {
        inputsignal[i] = inputdata[i];
    }
    const int n = 10;
    const double bta = 5.0;
    if (upfactor <= 0 || downfactor <= 0)
        throw std::runtime_error("factors must be positive integer");
    int gcd = getgcd(upfactor, downfactor);
    upfactor /= gcd;
    downfactor /= gcd;

    if (upfactor == downfactor)
    {
        outputdata = inputdata;
        return outputdata;
    }

    int inputsize = inputsignal.size();
    outputsignal.clear();
    int outputsize = quotientceil(inputsize * upfactor, downfactor);
    outputsignal.reserve(outputsize);
    int maxfactor = max(upfactor, downfactor);
    double firlsfreq = 1.0 / 2.0 / static_cast<double>(maxfactor);
    int length = 2 * n * maxfactor + 1;
    double firlsfreqs[] = {0.0, 2.0 * firlsfreq, 2.0 * firlsfreq, 1.0};
    vector<double> firlsfreqsv;
    firlsfreqsv.assign(firlsfreqs, firlsfreqs + 4);
    double firlsamplitude[] = {1.0, 1.0, 0.0, 0.0};
    vector<double> firlsamplitudev;
    firlsamplitudev.assign(firlsamplitude, firlsamplitude + 4);
    vector<double> coefficients;
    firls(length - 1, firlsfreqsv, firlsamplitudev, coefficients);
    vector<double> window;
    // kaiser(length, bta, window);

    if (readwindow(window, path) == -1)
        return outputdata;

    int coefficientssize = coefficients.size();
    for (int i = 0; i < coefficientssize; i++)
        coefficients[i] *= upfactor * window[i];

    int lengthhalf = (length - 1) / 2;
    int nz = downfactor - lengthhalf % downfactor;
    vector<double> h;
    h.reserve(coefficientssize + nz);
    for (int i = 0; i < nz; i++)
        h.push_back(0.0);
    for (int i = 0; i < coefficientssize; i++)
        h.push_back(coefficients[i]);
    int hsize = h.size();
    lengthhalf += nz;
    int delay = lengthhalf / downfactor;
    nz = 0;
    while (quotientceil((inputsize - 1) * upfactor + hsize + nz, downfactor) - delay < outputsize)
        nz++;
    for (int i = 0; i < nz; i++)
        h.push_back(0.0);
    vector<double> y;

    upfirdn(upfactor, downfactor, inputsignal, h, y);
    for (int i = delay; i < outputsize + delay; i++)
    {
        outputsignal.push_back(y[i]);
    }
    for (int i = 0; i < outputsignal.size(); i++)
    {
        outputdata[i] = outputsignal[i];
    }
    return outputdata;
}
