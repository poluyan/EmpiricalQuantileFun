/**************************************************************************

   Copyright © 2018 Sergey Poluyan <svpoluyan@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

**************************************************************************/

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <random>
#include <algorithm>


#include "timer.h"
#include "trie_based.h"
#include "quantile.h"
#include "test.h"

float empirical_cdf(std::vector<float> &sorted_sample, float val)
{
    auto pos = std::lower_bound(sorted_sample.begin(), sorted_sample.end(), val);
    return std::distance(sorted_sample.begin(), pos)/float(sorted_sample.size());
}

int empirical_cdf_int(std::vector<float> &sorted_sample, float val)
{
    auto pos = std::lower_bound(sorted_sample.begin(), sorted_sample.end(), val);
    return std::distance(sorted_sample.begin(), pos);
}

float empirical_qantile_from_real_cdf(std::vector<std::pair<float,float>> &cdf, float val01)
{
    auto pos = std::lower_bound(cdf.begin(), cdf.end(), val01,
                                [](const std::pair<float,float> &l, const float &r)
    {
        return l.second < r;
    });

    size_t stop = std::distance(cdf.begin(), pos);
    size_t start = stop - 1;
    /*if(std::distance(cdf.begin(), pos) == 0)
    {
        start = 0;
        stop = 1;
    }*/
    float x0 = cdf[start].first, x1 = cdf[stop].first, y0 = cdf[start].second, y1 = cdf[stop].second;
    return x0 + (val01 - y0)*(x1-x0)/(y1-y0);
}

float empirical_qantile_1d(std::vector<float> &cdf, float val01)
{
    auto pos = std::lower_bound(cdf.begin(), cdf.end(), val01);
    //return cdf[std::distance(cdf.begin(), pos)];

    return cdf[std::distance(cdf.begin(), pos)];
}


std::pair<size_t, float> ecdf1d_pair_fromgrid_trie(const std::vector<std::pair<int,int>> &sample, size_t sample_size, const std::vector<float> &grid, float val01)
{
    size_t l = 0, r = grid.size() - 1;

    size_t m = 0, index1 = 0, index2 = 0;
    float cdf1, cdf2;

    while(l <= r)
    {
        m = l + (r - l) / 2;

        index1 = 0;
        index2 = 0;
        for(size_t i = 0, n = sample.size(); i != n; ++i)
        {
            if(static_cast<size_t>(sample[i].first) < m)
            {
                index1 += sample[i].second;
            }
            if(static_cast<size_t>(sample[i].first) < m + 1)
            {
                index2 += sample[i].second;
            }
        }
        cdf1 = index1/float(sample_size);
        cdf2 = index2/float(sample_size);

        if((val01 > cdf1) && (val01 < cdf2))
            break;

        if(val01 > cdf1)
            l = m + 1;
        else
            r = m - 1;
    }

    float x0 = grid[m], y0 = cdf1, x1 = grid[m + 1], y1 = cdf2;
    return std::make_pair(m, x0 + (val01 - y0) * (x1 - x0) / (y1 - y0));
}

void ecdfNd_one_MultipleGrids_fromgrid_Trie(trie_based::TrieBased<trie_based::NodeCount<int>,int> &sample,
        const std::vector<std::vector<float> > &grids,
        const std::vector<float> &val01,
        std::vector<float> &rez)
{
    auto p = sample.root.get();
    for(size_t i = 0; i != val01.size(); i++)
    {
        std::vector<std::pair<int,int>> row;
        int cc = 0;
        for(size_t j = 0; j != p->children.size(); j++)
        {
            row.push_back(std::make_pair(p->children[j]->index,p->children[j]->count));
            cc += p->children[j]->count;
        }

        auto rez2 = ecdf1d_pair_fromgrid_trie(row,cc,grids[i],val01[i]);
        rez[i] = rez2.second;

        int index = 0;
        for(size_t j = 1; j < p->children.size(); j++)
        {
            if(static_cast<size_t>(p->children[j]->index) == rez2.first)
                index = j;
        }
        p = p->children[index].get();
    }
}

void implicit_quantile(std::vector<std::vector<int> > &sample, std::vector<std::vector<float> > &grids)
{
    trie_based::TrieBased<trie_based::NodeCount<int>,int> sample_trie(grids.size());
    for(auto i : sample)
        sample_trie.insert(i);
    sample_trie.fill_tree_count();

    std::cout << "total samples " << sample_trie.get_total_count() << std::endl;

    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    timer::Timer time_cpp11;
    time_cpp11.reset();
    std::vector<std::vector<float> > sampled;
    long long nrolls = 1e+2;  // number of experiments

    std::vector<float> temp1(grids.size());
    std::vector<float> temp2(temp1.size());
    for(long long i = 0; i != nrolls; ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        ecdfNd_one_MultipleGrids_fromgrid_Trie(sample_trie,grids,temp1,temp2);
        sampled.push_back(temp2);
    }
    std::cout << "total time: " << time_cpp11.elapsed_seconds() << std::endl;
    std::cout << "time per transform: " << time_cpp11.elapsed_seconds()/double(nrolls) << std::endl;
    data_io::write_default2d("maps/sampled_implicit.dat",sampled,5);

    sample_trie.remove_tree();
}

float getval(std::vector<float> &sample, std::vector<float> &grid, float val01)
{
    size_t count = grid.size() - 1, step, c1 = 0, c2 = 0;
    float f1, f2, n = sample.size();
    std::vector<float>::iterator first = grid.begin(), it;
    while(count > 0)
    {
        it = first;
        step = count / 2;
        std::advance(it, step);
        c1 = std::count_if(sample.begin(), sample.end(), [&it](const float &i)
        {
            return i < *it;
        });
        c2 = std::count_if(sample.begin(), sample.end(), [&it](const float &i)
        {
            return i < *(it + 1);
        });
        f1 = c1/n;
        f2 = c2/n;
        if(val01 > f1 && val01 < f2)
            break;
        if(f1 < val01)
        {
            first = ++it;
            count -= step + 1;
        }
        else
            count = step;
    }
    if(c1 == c2)
        return it == grid.begin() ? grid.front() : grid.back();
    return *it + (val01 - f1) * (*(it + 1) - *it) / (f2 - f1);
}

std::pair<size_t, float> ecdf1d_pair_fromgrid_trie2(const std::vector<std::pair<int,int>> &sample, size_t sample_size, const std::vector<float> &grid, float val01)
{
    size_t m = 0, count = grid.size() - 1, step, index1 = 0, index2 = 0;
    float cdf1, cdf2;
    std::vector<float>::const_iterator first = grid.begin(), it;

    while(count > 0)
    {
        it = first;
        step = count / 2;
        std::advance(it, step);
        m = std::distance(grid.begin(), it);

        index1 = 0;
        index2 = 0;
        for(const auto & i : sample)
        {
            if(static_cast<size_t>(i.first) < m)
            {
                index1 += i.second;
            }
            if(static_cast<size_t>(i.first) < m + 1)
            {
                index2 += i.second;
            }
        }
        cdf1 = index1/float(sample_size);
        cdf2 = index2/float(sample_size);

        if(val01 > cdf1 && val01 < cdf2)
            break;

        if(cdf1 < val01)
        {
            first = ++it;
            count -= step + 1;
        }
        else
            count = step;
    }
    if(index1 == index2)
        return it == grid.begin() ? std::make_pair(0, grid.front()) : std::make_pair(int(grid.size() - 1), grid.back());
    float x0 = grid[m], y0 = cdf1, x1 = grid[m + 1], y1 = cdf2;
    return std::make_pair(m, x0 + (val01 - y0) * (x1 - x0) / (y1 - y0));
}


void one_dim_check()
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::normal_distribution<float> norm(0.0,1.0);
    std::uniform_real_distribution<float> ureal(0.0,1.0);

    std::vector<float> grid = {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<float> sample = {-1.5, 1.5, 2.5};
    std::vector<std::pair<int,int>> ssample = {std::pair<int,int>{1,1}, std::pair<int,int>{4,1},std::pair<int,int>{5,1}};

//    std::vector<float> sample = {-2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5};

    std::vector<float> sampled1,sampled2,sampled3,sampled4;
    for(size_t i = 0; i != 100; i++)
    {
        float u = ureal(generator);
        sampled1.push_back(getval(sample, grid, u));
        //sampled2.push_back(getval3(sample, grid, u));
//        sampled3.push_back(ecdf1d_pair_fromgrid_trie(ssample, ssample.size(), grid, u).second);
        //sampled4.push_back(ecdf1d_pair_fromgrid_trie2(ssample, ssample.size(), grid, u).second);
    }
    //std::cout << getval(sample, grid, 0.999999) << std::endl;
    //std::cout << getval(sample, grid, 1.0) << std::endl;
    /*std::cout << getval2(sample, grid, 0.0) << std::endl;
    std::cout << getval2(sample, grid, 0.001) << std::endl;
    std::cout << getval(sample, grid, 0.0) << std::endl;
    std::cout << std::endl;
    std::cout << getval2(sample, grid, 0.9999) << std::endl;
    //const std::vector<std::pair<int,int>> &sample, size_t sample_size, const std::vector<float> &grid, float val01)
    std::cout << ecdf1d_pair_fromgrid_trie(ssample, ssample.size(), grid, 0.9999).second << std::endl;
    std::cout << getval(sample, grid, 0.9999) << std::endl;*/

//    std::cout.precision(15);
    //std::cout << std::fixed << getval(sample,grid, 0.999) << std::endl;
//    std::cout << '\n' << std::endl;
//    std::cout << std::fixed << getval3(sample,grid, /*0.350898*/0.9999) << std::endl;
//    std::cout << ecdf1d_pair_fromgrid_trie2(ssample, ssample.size(), grid, 0.9999).second << std::endl;
//    std::cout << std::fixed << getval3(sample,grid, /*0.350898*/0.0) << std::endl;
//    std::cout << ecdf1d_pair_fromgrid_trie2(ssample, ssample.size(), grid, 0.0).second << std::endl;
//    std::cout << std::fixed << getval3(sample,grid, /*0.350898*/1.0) << std::endl;
//    std::cout << ecdf1d_pair_fromgrid_trie2(ssample, ssample.size(), grid, 1.0).second << std::endl;
    //std::cout << std::fixed << getval2(sample,grid,0.350898) << std::endl;
    /*std::cout << std::fixed << getval2(sample,grid,0.0) << std::endl;
    std::cout << std::fixed << getval(sample,grid,0.0) << std::endl;*/

//    print2file("maps/1d1.dat",sampled1,1);
//    print2file("maps/1d2.dat",sampled2,1);
//    print2file("maps/1d3.dat",sampled3,1);
//    print2file("maps/1d4.dat",sampled4,1);

}


void simple_empirical_1d()
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::normal_distribution<float> norm(0.0,1.0);
    std::vector<float> sample(50);
    for(size_t i = 0; i != sample.size(); i++)
    {
        sample[i] = norm(generator);
    }
    std::vector<float> s_x = sample;
    std::sort(s_x.begin(), s_x.end());

    size_t N = 100;
    std::vector<float> cdf(N);

    float startp = -5.0;
    float endp = 5.0;
    float es = endp - startp;

    for(size_t i = 0; i != N; i++)
    {
        cdf[i] = empirical_cdf(s_x, startp + es * i / N);
        //std::cout << cdf[i] << std::endl;
    }
    data_io::write_default1d("maps/1d/quantile1d.dat", cdf, 1, 5);
}

float objective_function(float x)
{
    //return -1000.0+std::pow(x, 2.0);
    float A = 2000.0;
    return 100000 + std::pow(x-5.0,2.0) + A*(1.0-cos(acos(-1.0)*(x-5.0)/5.0));
}

void simple1d_example()
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> urand(-50.0,50.0);
    std::vector<std::vector<float>> sample;
    for(size_t i = 0; i != 10; i++)
    {
        std::vector<float> temp(2);
        temp[0] = urand(generator);
        temp[1] = objective_function(temp[0]);
        sample.push_back(temp);
    }
    data_io::write_default2d("maps/1d/pdf.dat", sample, 4);

    auto max_pdf = *std::max_element(sample.begin(), sample.end(),
                                     [](const std::vector<float> &a,const std::vector<float> &b)
    {
        return a[1] < b[1];
    });
    auto min_pdf = *std::min_element(sample.begin(), sample.end(),
                                     [](const std::vector<float> &a,const std::vector<float> &b)
    {
        return a[1] < b[1];
    });
    std::cout << max_pdf.back() << '\t' << min_pdf.back() << std::endl;

    if(max_pdf.back() > 0 && min_pdf.back() < 0)
    {
        for(size_t i = 0; i != sample.size(); i++)
        {
            sample[i][1] -= max_pdf.back();
        }
    }
    else if(max_pdf.back() > 0 && min_pdf.back() > 0)
    {
        for(size_t i = 0; i != sample.size(); i++)
        {
            sample[i][1] -= max_pdf.back();
            //sample[i][1] = sample[i][1] - min_pdf.back();
            //sample[i][1] = sample[i][1] - max_pdf.back() + min_pdf.back();
        }
    }
    else if(max_pdf.back() < 0 && min_pdf.back() < 0)
    {
        for(size_t i = 0; i != sample.size(); i++)
        {
            sample[i][1] += std::abs(max_pdf.back());
        }
    }

    for(size_t i = 0; i != sample.size(); i++)
    {
        sample[i][1] = -sample[i][1];
    }

    std::sort(sample.begin(), sample.end(),
              [](const std::vector<float>& x, const std::vector<float>& y)
    {
        return x.front() < y.front();
    });

    auto new_max_pdf = *std::max_element(sample.begin(), sample.end(),
                                         [](const std::vector<float> &a,const std::vector<float> &b)
    {
        return a[1] < b[1];
    });
    auto new_min_pdf = *std::min_element(sample.begin(), sample.end(),
                                         [](const std::vector<float> &a,const std::vector<float> &b)
    {
        return a[1] < b[1];
    });
    std::cout << new_max_pdf.back() << '\t' << new_min_pdf.back() << std::endl;

    for(size_t i = 0; i != sample.size(); i++)
    {
        sample[i][1] /= new_max_pdf.back();
    }

    data_io::write_default2d("maps/1d/pdf.dat", sample, 4);

    std::vector<float> s_x;
    for(size_t i = 0; i != sample.size(); i++)
    {
        size_t n = static_cast<size_t>(sample[i][1]*100);
        std::cout << n << std::endl;
        for(size_t j = 0; j != n; j++)
            s_x.push_back(sample[i][0]);
    }
    std::sort(s_x.begin(), s_x.end());

    data_io::write_default1d("maps/1d/sorted.dat", s_x, 1, 4);

    size_t N = 100000;
    std::vector<float> cdf(N);
    std::vector<std::pair<int,int>> cdf_int;
    std::vector<std::pair<float,float>> cdf_float;

    float startp = -50.0;
    float endp = 50.0;
    float es = endp - startp;

    for(size_t i = 0; i != N; i++)
    {
        cdf[i] = empirical_cdf(s_x, startp + es * i / N);
        //std::cout << cdf[i] << std::endl;
    }
    data_io::write_default1d("maps/1d/cdf.dat", cdf, 1, 4);


    for(size_t i = 0; i != N; i++)
    {
        int t = empirical_cdf_int(s_x, startp + es * i / N);
        if(std::find_if(cdf_int.begin(),cdf_int.end(),[&t](const std::pair<int,int> &l)
    {
        return t == l.second;
    }) == cdf_int.end())
        {
            cdf_int.push_back(std::make_pair(i, t));
            cdf_float.push_back(std::make_pair(startp + es * i / N, empirical_cdf(s_x, startp + es * i / N)));
        }
    }
    std::cout << cdf_int.size() << std::endl;
    for(const auto & i : cdf_int)
        std::cout << i.first << '\t' << i.second << std::endl;

    std::vector<std::vector<float>> cdf_simple;
    for(const auto & i : cdf_float)
    {
        std::cout << i.first << '\t' << i.second << std::endl;
        cdf_simple.push_back(std::vector<float> {i.first, i.second});
    }
    data_io::write_default2d("maps/1d/cdf_simple.dat", cdf_simple, 5);


    std::uniform_real_distribution<float> urand01(0.0,1);
    std::vector<float> sampled;
    for(size_t i = 0; i != 1e4; i++)
    {
        sampled.push_back(empirical_qantile_from_real_cdf(cdf_float, urand01(generator)));
    }
    data_io::write_default1d("maps/1d/sampled.dat", sampled, 1, 5);

    N = 1000;
    std::vector<float> quant(N);

    startp = 0.0;
    endp = 1.0;
    es = endp - startp;
    for(size_t i = 0; i != N; i++)
    {
        quant[i] = empirical_qantile_1d(cdf, startp + es * i / N);
    }
    data_io::write_default1d("maps/1d/quant.dat", cdf, 1, 4);
}

int main()
{ 
//
//    std::cout << tt.front() == tt.back() << std::endl;

//    simple_empirical_1d();
//    simple1d_example();

//    test_1d1();
//    test_1d2();
//    test_1d3();
//    test_1d4();


    test_2d1();
//    test_2d2();

//    test_grid_10d();
}
