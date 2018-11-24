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
#include <random>

void implicit_quantile_class(float lb, float ub, std::vector<size_t> gridn,std::vector<std::vector<int> > &sample, size_t nrolls)
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    timer::Timer time_cpp11;
    time_cpp11.reset();
    std::vector<std::vector<float> > sampled;
    std::vector<std::vector<float> > values01;

    empirical_quantile::ImplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn, sample);

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
    write_default2d("maps/values01.dat", values01, 4);
    write_default2d("maps/sampled_implicit_class.dat", sampled, 4);
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
//    explicit_quantile(sample_explicit, grids);
//    implicit_quantile(sample_implicit, grids);

    implicit_quantile_class(-3, 3, grid_number, sample_implicit, 2e+3);
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
}
