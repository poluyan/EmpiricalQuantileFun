#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <random>
#include <algorithm>

#include "print2file.h"
#include "timer.h"

double empirical_qantile_1d_sorted(std::vector<double> &sorted_sample, double val)
{
    auto pos = std::lower_bound(sorted_sample.begin(), sorted_sample.end(), val);
    return std::distance(sorted_sample.begin(), pos)/double(sorted_sample.size());
}

std::pair<size_t, float> ecdf1d_pair(const std::vector<float> &sample, const std::vector<float> &grid, float val01)
{
    //std::vector<float> sorted_sample = sample;
    //std::sort(sorted_sample.begin(), sorted_sample.end());

    //print2file("maps/sample1d.dat",sorted_sample,1);

    size_t l = 0;
    size_t r = grid.size() - 1;

    size_t m = 0, index1 = 0, index2 = 0;
    float cdf1, cdf2;

    while(l <= r)
    {
        m = l + (r - l) / 2;

        //auto pos1 = std::lower_bound(sorted_sample.begin(), sorted_sample.end(), grid[m]);
        //auto pos2 = std::lower_bound(sorted_sample.begin(), sorted_sample.end(), grid[m + 1]);
        //index1 = std::distance(sorted_sample.begin(), pos1);
        //index2 = std::distance(sorted_sample.begin(), pos2);
        //cdf1 = index1/float(sorted_sample.size());
        //cdf2 = index2/float(sorted_sample.size());

        index1 = 0;
        index2 = 0;
        for(size_t i = 0, n = sample.size(); i != n; ++i)
        {
            if(sample[i] < grid[m])
            {
                ++index1;
            }
            if(sample[i] < grid[m + 1])
            {
                ++index2;
            }
        }
        cdf1 = index1/float(sample.size());
        cdf2 = index2/float(sample.size());

        //std::cout << grid[m] << '\t' << grid[m + 1] << '\t' << cdf1 << '\t' << cdf2 << std::endl;

        if((val01 < cdf2) && (val01 > cdf1))
            break;

        if(val01 > cdf1)
            l = m + 1;
        else
            r = m - 1;
    }

    float x0 = grid[m], y0 = cdf1, x1 = grid[m + 1], y1 = cdf2;

    //if(index1 == index2 || m == grid.size() - 1 || m == grid.size() || m < 0)
    //{
    //    std::cout << "index1 == index2\t" << index1 << '\t' << index2 << std::endl;
    //    std::cin.get();
    //}
    return std::make_pair(m, x0 + (val01 - y0) * (x1 - x0) / (y1 - y0));
}

void ecdfNd_one_MultipleGrids(const std::vector<std::vector<float> > &sample,
                              const std::vector<std::vector<float> > &grids,
                              const std::vector<float> &val01,
                              std::vector<float> &rez)
{
    std::vector<size_t> m;
    //float delta = std::abs(grid[1] - grid[0])/2.0;
    //if((sample[j][k] > rez[i - 1] - delta)  && sample[j][k] < (rez[i - 1] + delta))

    for(size_t i = 0, g = val01.size(); i != g; i++)
    {
        std::vector<float> row(sample.size());
        size_t index = 0;
        for(size_t j = 0, n = sample.size(); j != n; j++)
        {
            bool flag = true;
            for(size_t k = 0, t = m.size(); k != t; k++)
            {
                //if((sample[j][k] > grid[m[k]] - i*delta) && (sample[j][k] < grid[m[k] + 1] + i*delta))

                //std::cout << sample[j][k] << '\t' << grids[k][m[k]] << '\t' << sample[j][k] << '\t' << grids[k][m[k] + 1] << std::endl;
                //std::cout << k + m[k] << '\t' << k + m[k] + 1 << std::endl;
                //std::cin.get();
                if(sample[j][k] > grids[k][m[k]] && sample[j][k] < grids[k][m[k] + 1])
                {
                    flag = true;
                }
                else
                {
                    flag = false;
                    break;
                }
            }
            if(flag)
            {
                row[index] = sample[j][i];
                ++index;
                //row.push_back(sample[j][i]);
            }
        }
        row.resize(index);
        //std::cout << row.size() << std::endl;
        //std::cin.get();
        //if(row.empty())
        //{
        //    std::cout << "row empty!" << row.size() << std::endl;

        //row.push_back(grid.front());
        //row.push_back(grid.back());
        //}

        //for(auto u : row)
        //    std::cout << u << std::endl;

        //std::cin.get();

        //std::cout << "here  v " << row.size() << '\t' << grids[i].size() << '\t' << val01[i] << std::endl;
        auto rez2 = ecdf1d_pair(row,grids[i],val01[i]);//ecdf1d_pair(row,grids[i],val01[i]);
        rez[i] = rez2.second;
        m.push_back(rez2.first);
    }
    //return rez;
}

void explicit_quantile(std::vector<std::vector<float> > &sample, std::vector<std::vector<float> > &grids)
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);
    
    timer::Timer time_cpp11;
    time_cpp11.reset();
    std::vector<std::vector<float> > sampled;
    long long nrolls = 2e+3;  // number of experiments

    std::vector<std::vector<float>> u01zvectors;

    std::vector<float> temp1(sample[0].size());
    std::vector<float> temp2(temp1.size());
    for(long long i = 0; i != nrolls; ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        u01zvectors.push_back(temp1);
        ecdfNd_one_MultipleGrids(sample,grids,temp1,temp2);
        sampled.push_back(temp2);
        //std::cout << std::endl;
    }
    std::cout << "total time: " << time_cpp11.elapsed_seconds() << std::endl;
    std::cout << "time per transform: " << time_cpp11.elapsed_seconds()/double(nrolls) << std::endl;
    print2file2d("maps/sampled_explicit.dat",sampled);
    print2file2d("maps/z.dat",u01zvectors);
}

void implicit_quantile(std::vector<std::vector<int> > &sample, std::vector<std::vector<float> > &grids)
{
    
}

int main()
{
//    std::mt19937_64 generator;
//    generator.seed(1);
//    std::normal_distribution<double> norm(0.0,1.0);
//
//    std::vector<double> sample(10);
//    for(size_t i = 0; i != sample.size(); i++)
//    {
//        sample[i] = norm(generator);
//    }
//    std::vector<double> s_x = sample;
//    std::sort(s_x.begin(), s_x.end());
//
//    size_t N = 100;
//    std::vector<double> cdf(N);
//
//    double startp = -5.0;
//    double endp = 5.0;
//    double es = endp - startp;
//
//    for(size_t i = 0; i != N; i++)
//    {
//        cdf[i] = empirical_qantile_1d_sorted(s_x, startp + es * i / N);
//        //std::cout << cdf[i] << std::endl;
//    }
//    print2file("maps/quantile1d.dat", cdf, 1);


    //

    std::vector<size_t> grid_number = {9, 10};

    std::vector<std::vector<float>> grids(grid_number.size());
    std::vector<float> dx(grid_number.size());

    for(size_t i = 0; i != grids.size(); i++)
    {
        std::vector<float> grid(grid_number[i] + 1);
        float startp = -3;
        float endp = 3;
        float es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/float(grid_number[i]);
        }
        grids[i] = grid;
        dx[i] = es/(float(grid_number[i])*2);
    }

    std::vector<std::vector<int>> sample_implicit;
    sample_implicit.push_back(std::vector{2,6});

    sample_implicit.push_back(std::vector{3,2});
    sample_implicit.push_back(std::vector{3,3});
    sample_implicit.push_back(std::vector{3,5});
    sample_implicit.push_back(std::vector{3,6});
    sample_implicit.push_back(std::vector{3,7});

    sample_implicit.push_back(std::vector{4,5});
    sample_implicit.push_back(std::vector{4,6});
    sample_implicit.push_back(std::vector{4,7});

    sample_implicit.push_back(std::vector{5,3});
    sample_implicit.push_back(std::vector{5,4});
    sample_implicit.push_back(std::vector{5,5});
    sample_implicit.push_back(std::vector{5,6});
    sample_implicit.push_back(std::vector{5,7});

    sample_implicit.push_back(std::vector{6,3});
    sample_implicit.push_back(std::vector{6,4});

    std::vector<std::vector<float>> sample_explicit;
    for(size_t i = 0; i != sample_implicit.size(); ++i)
    {
        std::vector<float> temp;
        for(size_t j = 0; j != sample_implicit[i].size(); ++j)
        {
            temp.push_back(grids[j][sample_implicit[i][j]] + dx[j]);
        }
        sample_explicit.push_back(temp);
    }

//    for(size_t i = 0; i != sample_explicit.size(); ++i)
//    {
//        for(size_t j = 0; j != sample_explicit[i].size(); ++j)
//        {
//            std::cout << sample_explicit[i][j] << '\t';
//        }
//        std::cout << std::endl;
//    }


    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(sample_explicit, grids);
    implicit_quantile(sample_implicit, grids);
}
