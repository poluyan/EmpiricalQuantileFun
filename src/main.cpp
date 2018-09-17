#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <random>
#include <algorithm>

template<typename T, typename A>
void print2file(std::string fname, std::vector<T, A> u, int step)
{
    std::ofstream fOut;
    fOut.open(fname.c_str());
    if (!fOut.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }
    fOut.precision(10);
    for (size_t i = 0; i < u.size(); i += step)
    {
        fOut << std::scientific << i << '\t' << u[i] << std::endl;
    }
    fOut.close();
    std::cout << fname << std::endl;
}

double empirical_qantile_1d_sorted(std::vector<double> &sorted_sample, double val)
{
    /*size_t i;
    for(i = 0; i != sorted_sample.size();)
    {
        if(sorted_sample[i] < val)
            i++;
        else
            break;
    }
    return i/double(sorted_sample.size());*/

    // if sorted
    auto pos = std::lower_bound(sorted_sample.begin(), sorted_sample.end(), val);
    return std::distance(sorted_sample.begin(), pos)/double(sorted_sample.size());
}

int main()
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::normal_distribution<double> norm(0.0,1.0);
    
    std::vector<double> sample(10);
    for(size_t i = 0; i != sample.size(); i++)
    {
        sample[i] = norm(generator);
    }
    std::vector<double> s_x = sample;
    std::sort(s_x.begin(), s_x.end());

    size_t N = 100;
    std::vector<double> cdf(N);

    double startp = -5.0;
    double endp = 5.0;
    double es = endp - startp;

    for(size_t i = 0; i != N; i++)
    {
        cdf[i] = empirical_qantile_1d_sorted(s_x, startp + es * i / N);
        //std::cout << cdf[i] << std::endl;
    }
    print2file("maps/quantile1d.dat", cdf, 1);
}
