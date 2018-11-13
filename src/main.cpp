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

#include "print2file.h"
#include "timer.h"
#include "trie_based.h"

#include "quantile.h"

double empirical_qantile_1d_sorted(std::vector<double> &sorted_sample, double val)
{
    auto pos = std::lower_bound(sorted_sample.begin(), sorted_sample.end(), val);
    return std::distance(sorted_sample.begin(), pos)/double(sorted_sample.size());
}

std::pair<size_t, float> ecdf1d_pair(const std::vector<float> &sample, const std::vector<float> &grid, float val01)
{
    size_t count = grid.size() - 1, step, c1 = 0, c2 = 0, m = 0;
    float f1, f2, n = sample.size();
    std::vector<float>::const_iterator first = grid.begin(), it;
    while(count > 0)
    {
        it = first;
        step = count / 2;
        std::advance(it, step);
        m = std::distance(grid.begin(), it);
        c1 = std::count_if(sample.begin(), sample.end(),
                           [&it](const float &v)
        {
            return v < *it;
        });
        c2 = std::count_if(sample.begin(), sample.end(),
                           [&it](const float &v)
        {
            return v < *(it + 1);
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
        return it == grid.begin() ? std::make_pair(size_t(0), grid.front()) : std::make_pair(grid.size() - 1,grid.back());
    return std::make_pair(m, *it + (val01 - f1) * (*(it + 1) - *it) / (f2 - f1));
}

void ecdfNd_one_MultipleGrids(const std::vector<std::vector<float> > &sample,
                              const std::vector<std::vector<float> > &grids,
                              const std::vector<float> &val01,
                              std::vector<float> &rez)
{
    std::vector<size_t> m;
    for(size_t i = 0, g = val01.size(); i != g; i++)
    {
        std::vector<float> row(sample.size());
        size_t index = 0;
        for(size_t j = 0, n = sample.size(); j != n; j++)
        {
            bool flag = true;
            for(size_t k = 0, t = m.size(); k != t; k++)
            {
                if(!(sample[j][k] > grids[k][m[k]] && sample[j][k] < grids[k][m[k] + 1]))
                {
                    flag = false;
                    break;
                }
            }
            if(flag)
            {
                row[index] = sample[j][i];
                ++index;
            }
        }
        row.resize(index);
        auto rez2 = ecdf1d_pair(row,grids[i],val01[i]);
        rez[i] = rez2.second;
        m.push_back(rez2.first);
    }
}

void explicit_quantile(std::vector<std::vector<float> > &sample, std::vector<std::vector<float> > &grids)
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    timer::Timer time_cpp11;
    time_cpp11.reset();
    std::vector<std::vector<float> > sampled;
    long long nrolls = 1e+2;  // number of experiments

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
    print2file2d("maps/sampled_explicit.dat",sampled, 5);
    print2file2d("maps/z.dat",u01zvectors, 5);
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
    print2file2d("maps/sampled_implicit.dat",sampled,5);

    sample_trie.remove_tree();
}

void implicit_quantile_class(float lb,
                             float ub,
                             std::vector<size_t> gridn,
                             std::vector<std::vector<int> > &sample)
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    timer::Timer time_cpp11;
    time_cpp11.reset();
    std::vector<std::vector<float> > sampled;
    long long nrolls = 1e+4;  // number of experiments

    empirical_quantile::ImplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn, sample);

    std::vector<float> temp1(gridn.size());
    std::vector<float> temp2(temp1.size());
    for(long long i = 0; i != nrolls; ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        quant.transform(temp1,temp2);
        sampled.push_back(temp2);
    }
    std::cout << "total time: " << time_cpp11.elapsed_seconds() << std::endl;
    std::cout << "time per transform: " << time_cpp11.elapsed_seconds()/double(nrolls) << std::endl;
    print2file2d("maps/sampled_implicit_class.dat",sampled,3);
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

void example_3d1()
{
    //    std::mt19937_64 generator;
//    generator.seed(1);
//    std::normal_distribution<float> norm(0.0,1.0);
//    std::vector<double> sample(50);
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

    //std::vector<size_t> grid_number = {9, 10};
    //std::vector<size_t> grid_number = {9, 12, 13, 8, 19, 44, 8, 4, 6, 7};
    std::vector<size_t> grid_number = {5,5,5};

    std::vector<std::vector<float>> grids(grid_number.size());
    std::vector<float> dx(grid_number.size());

    for(size_t i = 0; i != grids.size(); i++)
    {
        std::vector<float> grid(grid_number[i] + 1);
        float startp = 0;//-3;
        float endp = 5;//3
        float es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/float(grid_number[i]);
        }
        grids[i] = grid;
        dx[i] = es/(float(grid_number[i])*2);
    }

    std::vector<std::vector<int>> sample_implicit;
//    sample_implicit.push_back(std::vector{2,6});
//
//    sample_implicit.push_back(std::vector{3,2});
//    sample_implicit.push_back(std::vector{3,3});
//    sample_implicit.push_back(std::vector{3,5});
//    sample_implicit.push_back(std::vector{3,6});
//    sample_implicit.push_back(std::vector{3,7});
//
//    sample_implicit.push_back(std::vector{4,5});
//    sample_implicit.push_back(std::vector{4,6});
//    sample_implicit.push_back(std::vector{4,7});
//
//    sample_implicit.push_back(std::vector{5,3});
//    sample_implicit.push_back(std::vector{5,4});
//    sample_implicit.push_back(std::vector{5,5});
//    sample_implicit.push_back(std::vector{5,6});
//    sample_implicit.push_back(std::vector{5,7});
//
//    sample_implicit.push_back(std::vector{6,3});
//    sample_implicit.push_back(std::vector{6,4});


//    ///sample_implicit.push_back(std::vector{9, 12, 13, 8, 19, 44, 8, 4, 6, 7});
//    sample_implicit.push_back(std::vector{4, 2, 6, 3, 2, 3, 2, 1, 2, 4});
//    sample_implicit.push_back(std::vector{4, 3, 4, 3, 2, 3, 2, 1, 2, 5});
//    sample_implicit.push_back(std::vector{4, 4, 1, 3, 2, 3, 2, 1, 2, 1});
//    sample_implicit.push_back(std::vector{4, 5, 6, 3, 2, 3, 2, 1, 2, 2});
//    sample_implicit.push_back(std::vector{4, 5, 0, 3, 2, 3, 2, 1, 2, 3});
//    sample_implicit.push_back(std::vector{1, 2, 3, 4, 5, 6, 2, 0, 0, 2});

    // 0 a,       1 b,         2 c,           3 d,             4 e
    sample_implicit.push_back(std::vector{1,0,0}); // baa
    sample_implicit.push_back(std::vector{2,0,0}); // caa
    sample_implicit.push_back(std::vector{4,0,0}); // eaa
    sample_implicit.push_back(std::vector{0,2,0}); // aca
    sample_implicit.push_back(std::vector{4,4,0}); // eea
    sample_implicit.push_back(std::vector{4,3,0}); // eda
    sample_implicit.push_back(std::vector{3,3,0}); // dda
    
    sample_implicit.push_back(std::vector{0,0,1}); // aab
    
    sample_implicit.push_back(std::vector{3,0,2}); // dac
    sample_implicit.push_back(std::vector{0,3,2}); // adc
    
    sample_implicit.push_back(std::vector{0,3,3}); // add
    
    sample_implicit.push_back(std::vector{2,0,4}); // cae
    sample_implicit.push_back(std::vector{2,1,4}); // cbe
    sample_implicit.push_back(std::vector{2,2,4}); // cce


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
//    explicit_quantile(sample_explicit, grids);
//    implicit_quantile(sample_implicit, grids);

    //implicit_quantile_class(-3, 3, grid_number, sample_implicit);
    implicit_quantile_class(0, 5, grid_number, sample_implicit);
}


int main()
{
//    std::mt19937_64 generator;
//    generator.seed(1);
//    std::normal_distribution<float> norm(0.0,1.0);
//    std::vector<double> sample(50);
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

    //std::vector<size_t> grid_number = {9, 10};
    //std::vector<size_t> grid_number = {9, 12, 13, 8, 19, 44, 8, 4, 6, 7};
    std::vector<size_t> grid_number = {3,3,3};

    std::vector<std::vector<float>> grids(grid_number.size());
    std::vector<float> dx(grid_number.size());

    for(size_t i = 0; i != grids.size(); i++)
    {
        std::vector<float> grid(grid_number[i] + 1);
        float startp = 0;//-3;
        float endp = 3;//3
        float es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/float(grid_number[i]);
        }
        grids[i] = grid;
        dx[i] = es/(float(grid_number[i])*2);
    }

    std::vector<std::vector<int>> sample_implicit;
//    sample_implicit.push_back(std::vector{2,6});
//
//    sample_implicit.push_back(std::vector{3,2});
//    sample_implicit.push_back(std::vector{3,3});
//    sample_implicit.push_back(std::vector{3,5});
//    sample_implicit.push_back(std::vector{3,6});
//    sample_implicit.push_back(std::vector{3,7});
//
//    sample_implicit.push_back(std::vector{4,5});
//    sample_implicit.push_back(std::vector{4,6});
//    sample_implicit.push_back(std::vector{4,7});
//
//    sample_implicit.push_back(std::vector{5,3});
//    sample_implicit.push_back(std::vector{5,4});
//    sample_implicit.push_back(std::vector{5,5});
//    sample_implicit.push_back(std::vector{5,6});
//    sample_implicit.push_back(std::vector{5,7});
//
//    sample_implicit.push_back(std::vector{6,3});
//    sample_implicit.push_back(std::vector{6,4});


//    ///sample_implicit.push_back(std::vector{9, 12, 13, 8, 19, 44, 8, 4, 6, 7});
//    sample_implicit.push_back(std::vector{4, 2, 6, 3, 2, 3, 2, 1, 2, 4});
//    sample_implicit.push_back(std::vector{4, 3, 4, 3, 2, 3, 2, 1, 2, 5});
//    sample_implicit.push_back(std::vector{4, 4, 1, 3, 2, 3, 2, 1, 2, 1});
//    sample_implicit.push_back(std::vector{4, 5, 6, 3, 2, 3, 2, 1, 2, 2});
//    sample_implicit.push_back(std::vector{4, 5, 0, 3, 2, 3, 2, 1, 2, 3});
//    sample_implicit.push_back(std::vector{1, 2, 3, 4, 5, 6, 2, 0, 0, 2});

    // 0 a,       1 b,         2 c,           3 d,             4 e
    sample_implicit.push_back(std::vector{0,0,0}); // aaa
    sample_implicit.push_back(std::vector{0,0,0}); // aab
    sample_implicit.push_back(std::vector{0,0,0}); // aac
    
    sample_implicit.push_back(std::vector{0,2,0}); // aba
    sample_implicit.push_back(std::vector{4,4,0}); // abb
    sample_implicit.push_back(std::vector{4,3,0}); // abc
    
    sample_implicit.push_back(std::vector{3,3,0}); // aca
    sample_implicit.push_back(std::vector{3,3,0}); // acb
    sample_implicit.push_back(std::vector{3,3,0}); // acc
    
    sample_implicit.push_back(std::vector{3,3,0}); // baa
    sample_implicit.push_back(std::vector{3,3,0}); // bab
    sample_implicit.push_back(std::vector{3,3,0}); // bac
    
    sample_implicit.push_back(std::vector{3,3,0}); // bba
    sample_implicit.push_back(std::vector{3,3,0}); // bbb
    sample_implicit.push_back(std::vector{3,3,0}); // bbc
    
    sample_implicit.push_back(std::vector{3,3,0}); // bca
    sample_implicit.push_back(std::vector{3,3,0}); // bcb
    sample_implicit.push_back(std::vector{3,3,0}); // bcc
    
    sample_implicit.push_back(std::vector{3,3,0}); // caa
    sample_implicit.push_back(std::vector{3,3,0}); // cab
    sample_implicit.push_back(std::vector{0,0,1}); // cac
    
    sample_implicit.push_back(std::vector{0,0,1}); // cba
    sample_implicit.push_back(std::vector{0,0,1}); // cbb
    sample_implicit.push_back(std::vector{0,0,1}); // cbc
    
    sample_implicit.push_back(std::vector{0,0,1}); // cca
    sample_implicit.push_back(std::vector{0,0,1}); // ccb
    sample_implicit.push_back(std::vector{0,0,1}); // ccc
    
    /*
     digraph G {
  "root" -> "0"
  "root" -> "1"
  "root" -> "2"

  "0" -> "01"
  "0" -> "02"
  "0" -> "03"

  "1" -> "11"
  "1" -> "12"
  "1" -> "13"

  "2" -> "21"
  "2" -> "22"
  "2" -> "23"

  "01" -> "000"
  "01" -> "111"
  "01" -> "222"

  "02" -> "000"
  "02" -> "111"
  "02" -> "222"

  "03" -> "000"
  "03" -> "111"
  "03" -> "222"

  "11" -> "000"
  "11" -> "111"
  "11" -> "222"

  "12" -> "000"
  "12" -> "111"
  "12" -> "222"

  "13" -> "000"
  "13" -> "111"
  "13" -> "222"

  "13" -> "000"
  "13" -> "111"
  "13" -> "222"

  "21" -> "000"
  "21" -> "111"
  "21" -> "222"

  "22" -> "000"
  "22" -> "111"
  "22" -> "222"

  "23" -> "000"
  "23" -> "111"
  "23" -> "222"
}

 * /

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
//    explicit_quantile(sample_explicit, grids);
//    implicit_quantile(sample_implicit, grids);

    //implicit_quantile_class(-3, 3, grid_number, sample_implicit);
    implicit_quantile_class(0, 5, grid_number, sample_implicit);
}
