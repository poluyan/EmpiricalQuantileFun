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

#include "test.h"
#include "data_io.h"
#include <random>

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

        if(f1 < val01)
        {
            if(val01 < f2)
                break;

            first = ++it;
            count -= step + 1;
        }
        else
            count = step;
    }
    if(c1 == c2)
        return it == grid.begin() ? std::make_pair(size_t(0), grid.front()) : std::make_pair(grid.size() - 1, grid.back());
    //std::make_pair(m, *it + (val01 - f1) * (*(it + 1) - *it) / (f2 - f1));
    return std::make_pair(m, grid[m] + (val01 - f1) * (grid[m + 1] - grid[m]) / (f2 - f1));
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
void explicit_quantile(std::vector<std::vector<float> > &sample, std::vector<std::vector<float> > &grids, size_t nrolls)
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    timer::Timer time_cpp11;
    time_cpp11.reset();
    std::vector<std::vector<float> > sampled;

    std::vector<std::vector<float>> u01zvectors;

    std::vector<float> temp1(sample[0].size());
    std::vector<float> temp2(temp1.size());
    for(size_t i = 0; i != nrolls; ++i)
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
    data_io::write_default2d("maps/sampled_explicit.dat", sampled, 15);
    data_io::write_default2d("maps/z.dat", u01zvectors, 15);
}


void implicit_quantile_class(float lb, float ub, std::vector<size_t> gridn,std::vector<std::vector<int> > &sample, size_t nrolls)
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    std::vector<std::vector<float> > sampled;
    std::vector<std::vector<float> > values01;

    empirical_quantile::ImplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn, sample);

    timer::Timer time_cpp11;
    time_cpp11.reset();

    std::vector<float> temp1(gridn.size());
    std::vector<float> temp2(temp1.size());
    for(size_t i = 0; i != nrolls; ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }
    std::cout << "total time: " << time_cpp11.elapsed_seconds() << std::endl;
    std::cout << "time per transform: " << time_cpp11.elapsed_seconds()/double(nrolls) << std::endl;
    data_io::write_default2d("maps/values01.dat", values01, 15);
    data_io::write_default2d("maps/sampled_implicit_class.dat", sampled, 15);
}

void implicit_quantile_class_sorted(float lb, float ub, std::vector<size_t> gridn,std::vector<std::vector<int> > &sample, size_t nrolls)
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    std::vector<std::vector<float> > sampled;
    std::vector<std::vector<float> > values01;

    empirical_quantile::ImplicitQuantileSorted<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn, sample);

    timer::Timer time_cpp11;
    time_cpp11.reset();

    std::vector<float> temp1(gridn.size());
    std::vector<float> temp2(temp1.size());
    for(size_t i = 0; i != nrolls; ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }
    std::cout << "total time: " << time_cpp11.elapsed_seconds() << std::endl;
    std::cout << "time per transform: " << time_cpp11.elapsed_seconds()/double(nrolls) << std::endl;
    data_io::write_default2d("maps/values01.dat", values01, 15);
    data_io::write_default2d("maps/sampled_implicit_class_sorted.dat", sampled, 15);
}

void test_1d1()
{
    std::vector<size_t> grid_number = {6};

    std::vector<std::vector<float>> grids(grid_number.size());
    std::vector<float> dx(grid_number.size());

    for(size_t i = 0; i != grids.size(); i++)
    {
        std::vector<float> grid(grid_number[i] + 1);
        float startp = -2;
        float endp = 4;
        float es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/float(grid_number[i]);
        }
        grids[i] = grid;
        dx[i] = es/(float(grid_number[i])*2);
    }

    std::vector<std::vector<int>> sample_implicit;
    sample_implicit.push_back(std::vector{0});
    sample_implicit.push_back(std::vector{1});
    sample_implicit.push_back(std::vector{2});
    sample_implicit.push_back(std::vector{3});
    sample_implicit.push_back(std::vector{4});
    sample_implicit.push_back(std::vector{5});

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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(sample_explicit, grids, 1e+3);
    implicit_quantile_class(-2, 4, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class_sorted(-2, 4, grid_number, sample_implicit, 1e+3);
}

void test_1d2()
{
    std::vector<size_t> grid_number = {6};

    std::vector<std::vector<float>> grids(grid_number.size());
    std::vector<float> dx(grid_number.size());

    for(size_t i = 0; i != grids.size(); i++)
    {
        std::vector<float> grid(grid_number[i] + 1);
        float startp = -2;
        float endp = 4;
        float es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/float(grid_number[i]);
        }
        grids[i] = grid;
        dx[i] = es/(float(grid_number[i])*2);
    }

    std::vector<std::vector<int>> sample_implicit;
    sample_implicit.push_back(std::vector{4});
    sample_implicit.push_back(std::vector{1});
    sample_implicit.push_back(std::vector{3});
    //sample_implicit.push_back(std::vector{8});

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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(sample_explicit, grids, 1e+3);
    implicit_quantile_class(-2, 4, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class_sorted(-2, 4, grid_number, sample_implicit, 1e+3);
}


void test_2d1()
{
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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(sample_explicit, grids, 2e+3);
//    implicit_quantile(sample_implicit, grids);
//
    implicit_quantile_class(-3, 3, grid_number, sample_implicit, 2e+3);
    implicit_quantile_class_sorted(-3, 3, grid_number, sample_implicit, 2e+3);
}

void test_2d2()
{
    std::vector<size_t> grid_number;
    std::vector<std::vector<int>> sample_implicit;
    data_io::load_grid_and_sample("input/2d/grid.dat", "input/2d/points.dat", grid_number, sample_implicit);

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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    timer::Timer time_cpp11;
    time_cpp11.reset();
    explicit_quantile(sample_explicit, grids, 2e+1);
    std::cout << "--------->   total time: " << time_cpp11.elapsed_seconds() << std::endl;
    time_cpp11.reset();
    implicit_quantile_class(-3, 3, grid_number, sample_implicit, 2e+1);
    std::cout << "--------->   total time: " << time_cpp11.elapsed_seconds() << std::endl;
    time_cpp11.reset();
    implicit_quantile_class_sorted(-3, 3, grid_number, sample_implicit, 2e+1);
    std::cout << "--------->   total time: " << time_cpp11.elapsed_seconds() << std::endl;
}

void test_3d1()
{
    std::vector<size_t> grid_number = {5,5,5};

    std::vector<std::vector<float>> grids(grid_number.size());
    std::vector<float> dx(grid_number.size());

    for(size_t i = 0; i != grids.size(); i++)
    {
        std::vector<float> grid(grid_number[i] + 1);
        float startp = 0;
        float endp = 5;
        float es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/float(grid_number[i]);
        }
        grids[i] = grid;
        dx[i] = es/(float(grid_number[i])*2);
    }

    std::vector<std::vector<int>> sample_implicit;


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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
//    explicit_quantile(sample_explicit, grids);
//    implicit_quantile(sample_implicit, grids);

    implicit_quantile_class(0, 5, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class_sorted(0, 5, grid_number, sample_implicit, 1e+3);
}

void test_3d2()
{
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

    */

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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
//    explicit_quantile(sample_explicit, grids);
//    implicit_quantile(sample_implicit, grids);

    //implicit_quantile_class(-3, 3, grid_number, sample_implicit);
    implicit_quantile_class(0, 3, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class_sorted(0, 3, grid_number, sample_implicit, 1e+3);
}


void test_grid_10d()
{
    std::vector<size_t> grid_number;
    std::vector<std::vector<int>> sample_implicit;
    //400 temp1 = {0.99935, 0.546268, 0.140131, 0.692333, 0.441771, 0.890283, 0.0597646, 0.607688, 0.566813, 0.61283};
    data_io::load_grid_and_sample("input/grid_test/10000/grid.dat", "input/grid_test/10000/sample.dat", grid_number, sample_implicit);

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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    timer::Timer time_cpp11;
    time_cpp11.reset();
//    explicit_quantile(sample_explicit, grids, 1e+1);
    std::cout << "--------->   total time: " << time_cpp11.elapsed_seconds() << std::endl;
    time_cpp11.reset();
    implicit_quantile_class(-3, 3, grid_number, sample_implicit, 1e+5);
    std::cout << "--------->   total time: " << time_cpp11.elapsed_seconds() << std::endl;
    time_cpp11.reset();
    implicit_quantile_class_sorted(-3, 3, grid_number, sample_implicit, 1e+5);
    std::cout << "--------->   total time: " << time_cpp11.elapsed_seconds() << std::endl;
}

void test_1d3()
{
    std::vector<size_t> grid_number = {400};

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

    sample_implicit.push_back(std::vector{0});
    sample_implicit.push_back(std::vector{2});
    sample_implicit.push_back(std::vector{4});
    sample_implicit.push_back(std::vector{11});
    sample_implicit.push_back(std::vector{12});
    sample_implicit.push_back(std::vector{16});
    sample_implicit.push_back(std::vector{19});
    sample_implicit.push_back(std::vector{20});
    sample_implicit.push_back(std::vector{24});
    sample_implicit.push_back(std::vector{25});
    sample_implicit.push_back(std::vector{26});
    sample_implicit.push_back(std::vector{27});
    sample_implicit.push_back(std::vector{28});
    sample_implicit.push_back(std::vector{30});
    sample_implicit.push_back(std::vector{31});
    sample_implicit.push_back(std::vector{33});
    sample_implicit.push_back(std::vector{34});
    sample_implicit.push_back(std::vector{36});
    sample_implicit.push_back(std::vector{39});
    sample_implicit.push_back(std::vector{41});
    sample_implicit.push_back(std::vector{44});
    sample_implicit.push_back(std::vector{46});
    sample_implicit.push_back(std::vector{47});
    sample_implicit.push_back(std::vector{48});
    sample_implicit.push_back(std::vector{50});
    sample_implicit.push_back(std::vector{51});
    sample_implicit.push_back(std::vector{54});
    sample_implicit.push_back(std::vector{56});
    sample_implicit.push_back(std::vector{61});
    sample_implicit.push_back(std::vector{62});
    sample_implicit.push_back(std::vector{63});
    sample_implicit.push_back(std::vector{64});
    sample_implicit.push_back(std::vector{66});
    sample_implicit.push_back(std::vector{68});
    sample_implicit.push_back(std::vector{71});
    sample_implicit.push_back(std::vector{77});
    sample_implicit.push_back(std::vector{78});
    sample_implicit.push_back(std::vector{79});
    sample_implicit.push_back(std::vector{80});
    sample_implicit.push_back(std::vector{81});
    sample_implicit.push_back(std::vector{82});
    sample_implicit.push_back(std::vector{84});
    sample_implicit.push_back(std::vector{85});
    sample_implicit.push_back(std::vector{87});
    sample_implicit.push_back(std::vector{88});
    sample_implicit.push_back(std::vector{93});
    sample_implicit.push_back(std::vector{94});
    sample_implicit.push_back(std::vector{95});
    sample_implicit.push_back(std::vector{96});
    sample_implicit.push_back(std::vector{99});
    sample_implicit.push_back(std::vector{101});
    sample_implicit.push_back(std::vector{103});
    sample_implicit.push_back(std::vector{104});
    sample_implicit.push_back(std::vector{106});
    sample_implicit.push_back(std::vector{107});
    sample_implicit.push_back(std::vector{110});
    sample_implicit.push_back(std::vector{112});
    sample_implicit.push_back(std::vector{114});
    sample_implicit.push_back(std::vector{118});
    sample_implicit.push_back(std::vector{121});
    sample_implicit.push_back(std::vector{123});
    sample_implicit.push_back(std::vector{126});
    sample_implicit.push_back(std::vector{129});
    sample_implicit.push_back(std::vector{130});
    sample_implicit.push_back(std::vector{133});
    sample_implicit.push_back(std::vector{137});
    sample_implicit.push_back(std::vector{138});
    sample_implicit.push_back(std::vector{140});
    sample_implicit.push_back(std::vector{143});
    sample_implicit.push_back(std::vector{146});
    sample_implicit.push_back(std::vector{150});
    sample_implicit.push_back(std::vector{151});
    sample_implicit.push_back(std::vector{153});
    sample_implicit.push_back(std::vector{154});
    sample_implicit.push_back(std::vector{157});
    sample_implicit.push_back(std::vector{162});
    sample_implicit.push_back(std::vector{163});
    sample_implicit.push_back(std::vector{164});
    sample_implicit.push_back(std::vector{166});


    //sample_implicit.push_back(std::vector{8});

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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(sample_explicit, grids, 1e+3);
    implicit_quantile_class(-3, 3, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class_sorted(-3, 3, grid_number, sample_implicit, 1e+3);
}
