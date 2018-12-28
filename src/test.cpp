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

void explicit_quantile(float lb, float ub, std::vector<size_t> gridn, std::vector<std::vector<int> > &sample, size_t nrolls)
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    std::vector<std::vector<float> > sampled;
    std::vector<std::vector<float> > values01;

    empirical_quantile::ExplicitQuantile<int, float> quant(
        std::vector<float>(gridn.size(), lb),
        std::vector<float>(gridn.size(), ub),
        gridn);
    quant.set_sample(sample);

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

    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 0.0;
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }
    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 1.0;
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }

    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = 0.0;
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }
    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = 1.0;
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }

    for(size_t i = 0; i != temp1.size() - 1; ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = 0.0;
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }
    for(size_t i = 0; i != temp1.size() - 1; ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = 1.0;
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }

//    std::vector<float> a = {7.546996e-01,
//                            2.108486e-01,
//                            9.091859e-01,
//                            9.425417e-01,
//                            6.641316e-01,
//                            8.285044e-01,
//                            6.402776e-01,
//                            1.946763e-01,
//                            2.501838e-01,
//                            8.654959e-01
//                           };
//    std::vector<float> b(a.size(), 0);
//    quant.transform(a, b);
//    sampled.push_back(b);

//    a = {0.68};
//    b = {0.0};
//    quant.transform(a, b);
//    sampled.push_back(b);

    std::cout << "total time: " << time_cpp11.elapsed_seconds() << std::endl;
    std::cout << "time per transform: " << std::scientific << time_cpp11.elapsed_seconds()/double(sampled.size()) << std::endl;
    data_io::write_default2d("maps/sampled_explicit.dat", sampled, 5);
    //data_io::write_default2d("maps/z.dat", values01, 15);
}


void implicit_quantile_class(float lb, float ub, std::vector<size_t> gridn,std::vector<std::vector<int> > &sample, size_t nrolls)
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    std::vector<std::vector<float> > sampled;
    std::vector<std::vector<float> > values01;

    empirical_quantile::ImplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
    quant.set_sample(sample);

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

    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 0.0;
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }
    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 1.0;
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }

    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = 0.0;
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }
    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = 1.0;
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }

    for(size_t i = 0; i != temp1.size() - 1; ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = 0.0;
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }
    for(size_t i = 0; i != temp1.size() - 1; ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = 1.0;
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }

//    std::vector<float> a = {7.546996e-01,
//                            2.108486e-01,
//                            9.091859e-01,
//                            9.425417e-01,
//                            6.641316e-01,
//                            8.285044e-01,
//                            6.402776e-01,
//                            1.946763e-01,
//                            2.501838e-01,
//                            8.654959e-01
//                           };
//    std::vector<float> b(a.size(), 0);
//    quant.transform(a, b);
//    sampled.push_back(b);
//
//    a = {0.68};
//    b = {0.0};
//    quant.transform(a, b);
//    sampled.push_back(b);

    std::cout << "total time: " << time_cpp11.elapsed_seconds() << std::endl;
    std::cout << "time per transform: " << std::scientific << time_cpp11.elapsed_seconds()/double(sampled.size()) << std::endl;
//    data_io::write_default2d("maps/values01.dat", values01, 15);
    data_io::write_default2d("maps/sampled_implicit.dat", sampled, 5);
}

void implicit_quantile_class_sorted(float lb, float ub, std::vector<size_t> gridn, std::vector<std::vector<int> > &sample, size_t nrolls)
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    std::vector<std::vector<float> > sampled;
    std::vector<std::vector<float> > values01;

    empirical_quantile::ImplicitQuantileSorted<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
    quant.set_sample(sample);

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

    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 0.0;
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }
    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 1.0;
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }

    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = 0.0;
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }
    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = 1.0;
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }

    for(size_t i = 0; i != temp1.size() - 1; ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = 0.0;
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }
    for(size_t i = 0; i != temp1.size() - 1; ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = 1.0;
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }

//    std::vector<float> a = {7.546996e-01,
//                            2.108486e-01,
//                            9.091859e-01,
//                            9.425417e-01,
//                            6.641316e-01,
//                            8.285044e-01,
//                            6.402776e-01,
//                            1.946763e-01,
//                            2.501838e-01,
//                            8.654959e-01
//                           };
//    std::vector<float> b(a.size(), 0);
//    quant.transform(a, b);
//    sampled.push_back(b);

//    a = {0.68};
//    b = {0.0};
//    quant.transform(a, b);
//    sampled.push_back(b);

    std::cout << "total time: " << time_cpp11.elapsed_seconds() << std::endl;
    std::cout << "time per transform: " << std::scientific << time_cpp11.elapsed_seconds()/double(sampled.size()) << std::endl;
//    data_io::write_default2d("maps/values01.dat", values01, 15);
    data_io::write_default2d("maps/sampled_implicit_sorted.dat", sampled, 5);
}

void test_1d1()
{
    std::vector<size_t> grid_number = {6};

    std::vector<std::vector<int>> sample_implicit;
    sample_implicit.push_back(std::vector{0});
    sample_implicit.push_back(std::vector{1});
    sample_implicit.push_back(std::vector{2});
    sample_implicit.push_back(std::vector{3});
    sample_implicit.push_back(std::vector{4});
    sample_implicit.push_back(std::vector{5});

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(-2, 4, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class(-2, 4, grid_number, sample_implicit, 1e3);
    implicit_quantile_class_sorted(-2, 4, grid_number, sample_implicit, 1e3);
}

void test_1d5()
{
    std::vector<size_t> grid_number = {6};

    std::vector<std::vector<int>> sample_implicit;
    sample_implicit.push_back(std::vector{3});

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(-2, 4, grid_number, sample_implicit, 1e3);
    implicit_quantile_class(-2, 4, grid_number, sample_implicit, 1e3);
    implicit_quantile_class_sorted(-2, 4, grid_number, sample_implicit, 1e3);
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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(-2, 4, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class(-2, 4, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class_sorted(-2, 4, grid_number, sample_implicit, 1e+3);
}


void test_2d1()
{
    std::vector<size_t> grid_number = {9, 10};

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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n nrolls = 2e3
    explicit_quantile(-3, 3, grid_number, sample_implicit, 1);
    implicit_quantile_class(-3, 3, grid_number, sample_implicit, 1);
    implicit_quantile_class_sorted(-3, 3, grid_number, sample_implicit, 1);
}

void test_2d2()
{
    std::vector<size_t> grid_number;
    std::vector<std::vector<int>> sample_implicit;
    data_io::load_grid_and_sample("input/2d/grid.dat", "input/2d/points.dat", grid_number, sample_implicit);

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    timer::Timer time_cpp11;
    time_cpp11.reset();
    explicit_quantile(-3, 3, grid_number, sample_implicit, 2e+1);
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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(0, 5, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class(0, 5, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class_sorted(0, 5, grid_number, sample_implicit, 1e+3);
}

void test_3d2()
{
    std::vector<size_t> grid_number = {3,3,3};

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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(0, 3, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class(0, 3, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class_sorted(0, 3, grid_number, sample_implicit, 1e+3);
}


void test_grid_10d()
{
    std::vector<size_t> grid_number;
    std::vector<std::vector<int>> sample_implicit;
    //400 temp1 = {0.99935, 0.546268, 0.140131, 0.692333, 0.441771, 0.890283, 0.0597646, 0.607688, 0.566813, 0.61283};
    data_io::load_grid_and_sample("input/grid_test/10000/grid.dat", "input/grid_test/10000/sample.dat", grid_number, sample_implicit);

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    timer::Timer time_cpp11;
    time_cpp11.reset();
    explicit_quantile(-3, 3, grid_number, sample_implicit, 1e+2);
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

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(-3, 3, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class(-3, 3, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class_sorted(-3, 3, grid_number, sample_implicit, 1e+3);
}

void test_1d4()
{
    std::vector<size_t> grid_number = {400};

    std::vector<std::vector<float>> grids(grid_number.size());
    std::vector<float> dx(grid_number.size());

    /// must find 0.68

    std::vector<std::vector<int>> sample_implicit;

    sample_implicit.push_back(std::vector{0});
    sample_implicit.push_back(std::vector{2});
    sample_implicit.push_back(std::vector{3});
    sample_implicit.push_back(std::vector{4});
    sample_implicit.push_back(std::vector{5});
    sample_implicit.push_back(std::vector{6});
    sample_implicit.push_back(std::vector{7});
    sample_implicit.push_back(std::vector{8});
    sample_implicit.push_back(std::vector{9});
    sample_implicit.push_back(std::vector{10});
    sample_implicit.push_back(std::vector{11});
    sample_implicit.push_back(std::vector{12});
    sample_implicit.push_back(std::vector{15});
    sample_implicit.push_back(std::vector{16});
    sample_implicit.push_back(std::vector{22});
    sample_implicit.push_back(std::vector{24});
    sample_implicit.push_back(std::vector{25});
    sample_implicit.push_back(std::vector{29});
    sample_implicit.push_back(std::vector{32});
    sample_implicit.push_back(std::vector{33});
    sample_implicit.push_back(std::vector{35});
    sample_implicit.push_back(std::vector{38});
    sample_implicit.push_back(std::vector{39});
    sample_implicit.push_back(std::vector{41});
    sample_implicit.push_back(std::vector{43});
    sample_implicit.push_back(std::vector{48});
    sample_implicit.push_back(std::vector{51});
    sample_implicit.push_back(std::vector{54});
    sample_implicit.push_back(std::vector{55});
    sample_implicit.push_back(std::vector{56});
    sample_implicit.push_back(std::vector{58});
    sample_implicit.push_back(std::vector{59});
    sample_implicit.push_back(std::vector{62});
    sample_implicit.push_back(std::vector{65});
    sample_implicit.push_back(std::vector{69});
    sample_implicit.push_back(std::vector{73});
    sample_implicit.push_back(std::vector{76});
    sample_implicit.push_back(std::vector{77});
    sample_implicit.push_back(std::vector{78});
    sample_implicit.push_back(std::vector{80});
    sample_implicit.push_back(std::vector{81});
    sample_implicit.push_back(std::vector{82});
    sample_implicit.push_back(std::vector{83});
    sample_implicit.push_back(std::vector{84});
    sample_implicit.push_back(std::vector{85});
    sample_implicit.push_back(std::vector{89});
    sample_implicit.push_back(std::vector{93});
    sample_implicit.push_back(std::vector{96});
    sample_implicit.push_back(std::vector{98});
    sample_implicit.push_back(std::vector{99});
    sample_implicit.push_back(std::vector{101});
    sample_implicit.push_back(std::vector{103});
    sample_implicit.push_back(std::vector{104});
    sample_implicit.push_back(std::vector{105});
    sample_implicit.push_back(std::vector{107});
    sample_implicit.push_back(std::vector{109});
    sample_implicit.push_back(std::vector{112});
    sample_implicit.push_back(std::vector{115});
    sample_implicit.push_back(std::vector{117});
    sample_implicit.push_back(std::vector{119});
    sample_implicit.push_back(std::vector{121});
    sample_implicit.push_back(std::vector{124});
    sample_implicit.push_back(std::vector{126});
    sample_implicit.push_back(std::vector{127});
    sample_implicit.push_back(std::vector{129});
    sample_implicit.push_back(std::vector{137});
    sample_implicit.push_back(std::vector{138});
    sample_implicit.push_back(std::vector{139});
    sample_implicit.push_back(std::vector{142});
    sample_implicit.push_back(std::vector{143});
    sample_implicit.push_back(std::vector{144});
    sample_implicit.push_back(std::vector{146});
    sample_implicit.push_back(std::vector{147});
    sample_implicit.push_back(std::vector{149});
    sample_implicit.push_back(std::vector{153});
    sample_implicit.push_back(std::vector{154});
    sample_implicit.push_back(std::vector{155});
    sample_implicit.push_back(std::vector{156});
    sample_implicit.push_back(std::vector{158});
    sample_implicit.push_back(std::vector{162});
    sample_implicit.push_back(std::vector{163});
    sample_implicit.push_back(std::vector{164});
    sample_implicit.push_back(std::vector{165});
    sample_implicit.push_back(std::vector{166});
    sample_implicit.push_back(std::vector{171});
    sample_implicit.push_back(std::vector{172});
    sample_implicit.push_back(std::vector{174});
    sample_implicit.push_back(std::vector{176});
    sample_implicit.push_back(std::vector{177});
    sample_implicit.push_back(std::vector{179});
    sample_implicit.push_back(std::vector{184});
    sample_implicit.push_back(std::vector{185});
    sample_implicit.push_back(std::vector{189});
    sample_implicit.push_back(std::vector{190});
    sample_implicit.push_back(std::vector{193});
    sample_implicit.push_back(std::vector{194});
    sample_implicit.push_back(std::vector{195});
    sample_implicit.push_back(std::vector{197});
    sample_implicit.push_back(std::vector{200});
    sample_implicit.push_back(std::vector{201});
    sample_implicit.push_back(std::vector{202});
    sample_implicit.push_back(std::vector{204});
    sample_implicit.push_back(std::vector{208});
    sample_implicit.push_back(std::vector{209});
    sample_implicit.push_back(std::vector{210});
    sample_implicit.push_back(std::vector{211});
    sample_implicit.push_back(std::vector{213});
    sample_implicit.push_back(std::vector{215});
    sample_implicit.push_back(std::vector{216});
    sample_implicit.push_back(std::vector{218});
    sample_implicit.push_back(std::vector{222});
    sample_implicit.push_back(std::vector{223});
    sample_implicit.push_back(std::vector{224});
    sample_implicit.push_back(std::vector{225});
    sample_implicit.push_back(std::vector{226});
    sample_implicit.push_back(std::vector{228});
    sample_implicit.push_back(std::vector{231});
    sample_implicit.push_back(std::vector{234});
    sample_implicit.push_back(std::vector{239});
    sample_implicit.push_back(std::vector{240});
    sample_implicit.push_back(std::vector{246});
    sample_implicit.push_back(std::vector{247});
    sample_implicit.push_back(std::vector{249});
    sample_implicit.push_back(std::vector{251});
    sample_implicit.push_back(std::vector{252});
    sample_implicit.push_back(std::vector{254});
    sample_implicit.push_back(std::vector{255});
    sample_implicit.push_back(std::vector{257});
    sample_implicit.push_back(std::vector{258});
    sample_implicit.push_back(std::vector{259});
    sample_implicit.push_back(std::vector{262});
    sample_implicit.push_back(std::vector{268});
    sample_implicit.push_back(std::vector{269});
    sample_implicit.push_back(std::vector{271});
    sample_implicit.push_back(std::vector{275});
    sample_implicit.push_back(std::vector{277});
    sample_implicit.push_back(std::vector{280});
    sample_implicit.push_back(std::vector{281});
    sample_implicit.push_back(std::vector{283});
    sample_implicit.push_back(std::vector{284});
    sample_implicit.push_back(std::vector{286});
    sample_implicit.push_back(std::vector{287});
    sample_implicit.push_back(std::vector{289});
    sample_implicit.push_back(std::vector{290});
    sample_implicit.push_back(std::vector{292});
    sample_implicit.push_back(std::vector{293});
    sample_implicit.push_back(std::vector{298});
    sample_implicit.push_back(std::vector{299});
    sample_implicit.push_back(std::vector{300});
    sample_implicit.push_back(std::vector{304});
    sample_implicit.push_back(std::vector{306});
    sample_implicit.push_back(std::vector{311});
    sample_implicit.push_back(std::vector{315});
    sample_implicit.push_back(std::vector{320});
    sample_implicit.push_back(std::vector{324});
    sample_implicit.push_back(std::vector{326});
    sample_implicit.push_back(std::vector{327});
    sample_implicit.push_back(std::vector{331});
    sample_implicit.push_back(std::vector{332});
    sample_implicit.push_back(std::vector{334});
    sample_implicit.push_back(std::vector{337});
    sample_implicit.push_back(std::vector{338});
    sample_implicit.push_back(std::vector{341});
    sample_implicit.push_back(std::vector{342});
    sample_implicit.push_back(std::vector{344});
    sample_implicit.push_back(std::vector{346});
    sample_implicit.push_back(std::vector{347});
    sample_implicit.push_back(std::vector{348});
    sample_implicit.push_back(std::vector{349});
    sample_implicit.push_back(std::vector{352});
    sample_implicit.push_back(std::vector{353});
    sample_implicit.push_back(std::vector{354});
    sample_implicit.push_back(std::vector{356});
    sample_implicit.push_back(std::vector{357});
    sample_implicit.push_back(std::vector{358});
    sample_implicit.push_back(std::vector{359});
    sample_implicit.push_back(std::vector{360});
    sample_implicit.push_back(std::vector{361});
    sample_implicit.push_back(std::vector{364});
    sample_implicit.push_back(std::vector{366});
    sample_implicit.push_back(std::vector{367});
    sample_implicit.push_back(std::vector{369});
    sample_implicit.push_back(std::vector{370});
    sample_implicit.push_back(std::vector{373});
    sample_implicit.push_back(std::vector{376});
    sample_implicit.push_back(std::vector{378});
    sample_implicit.push_back(std::vector{379});
    sample_implicit.push_back(std::vector{380});
    sample_implicit.push_back(std::vector{381});
    sample_implicit.push_back(std::vector{387});
    sample_implicit.push_back(std::vector{388});
    sample_implicit.push_back(std::vector{389});
    sample_implicit.push_back(std::vector{391});
    sample_implicit.push_back(std::vector{393});
    sample_implicit.push_back(std::vector{394});
    sample_implicit.push_back(std::vector{395});
    sample_implicit.push_back(std::vector{396});
    sample_implicit.push_back(std::vector{397});
    sample_implicit.push_back(std::vector{398});
    sample_implicit.push_back(std::vector{399});

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(-3, 3, grid_number, sample_implicit, 1e3);
    implicit_quantile_class(-3, 3, grid_number, sample_implicit, 1e3);
    implicit_quantile_class_sorted(-3, 3, grid_number, sample_implicit, 1e3);
}


void test_Nd(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nsamples, size_t Nrolls)
{
    size_t max_sample_size = gridN.front();
    for(size_t i = 1; i != gridN.size(); i++)
    {
        max_sample_size *= gridN[i];
        if(max_sample_size > Nsamples)
            break;
    }
    if(Nsamples > max_sample_size)
    {
        std::cout << "Nsamples > max_sample_size" << std::endl;
        std::cout << Nsamples << " > " << max_sample_size << std::endl;
        return;
    }
    for(size_t i = 0; i != lb.size(); i++)
    {
        if(lb[i] > ub[i])
            return;
    }

    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    typedef trie_based::TrieBased<trie_based::NodeCount<int>,int> sample_type;
    std::shared_ptr<sample_type> sample = std::make_shared<sample_type>();

    std::vector<std::vector<int> > sample_implicit;

    std::vector<int> temp(gridN.size());
    for(size_t i = 0; i != Nsamples;)
    {
        for(size_t j = 0; j != gridN.size(); j++)
        {
            temp[j] = static_cast<int>(std::round(ureal01(generator)*(gridN[j] - 1.0)));
        }
        if(!sample->search(temp))
        {
            sample->insert(temp);
            sample_implicit.push_back(temp);
            i++;
        }
    }

    std::vector<std::vector<float> > values01;

    std::vector<float> temp1(gridN.size());
    std::vector<float> temp2(temp1.size());

    for(size_t i = 0; i != Nrolls; ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        values01.push_back(temp1);
    }

    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 0.0;
        values01.push_back(temp1);
    }
    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 1.0;
        values01.push_back(temp1);
    }

    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = 0.0;
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        values01.push_back(temp1);
    }
    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = 1.0;
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        values01.push_back(temp1);
    }

    for(size_t i = 0; i != temp1.size() - 1; ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = 0.0;
        }
        values01.push_back(temp1);
    }
    for(size_t i = 0; i != temp1.size() - 1; ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = 1.0;
        }
        values01.push_back(temp1);
    }
    data_io::write_default2d("maps/values01.dat", values01, 15);


//    values01.clear();

//    values01.push_back(std::vector<float>{7.128837108612061e-01	,
//    8.562362790107727e-01,	1.416520029306412e-01,	1.646486222743988e-01,	3.487232625484467e-01,
//    2.512360513210297e-01,	9.061432480812073e-01,	7.048301100730896e-01,	3.441495597362518e-01,
//	0.000000000000000e+00});
//
//    values01.push_back(
//    std::vector<float>{
//    5.726840496063232e-01,	3.626562356948853e-01,	8.979787230491638e-01,	3.496656417846680e-01,
//    1.190820038318634e-01,	4.006883800029755e-01,	7.749708294868469e-01,	4.034065902233124e-0,
//    7.621356248855591e-01,	2.899370491504669e-01});

//    values01.push_back(std::vector<float>{
//        3.842369019985199e-01,	3.113699853420258e-01,	5.914384126663208e-01,	6.040902137756348e-01,
//        7.430433034896851e-01,	8.774918913841248e-01,	1.933544427156448e-01,	3.517977595329285e-01,
//        9.651318788528442e-01,	3.330467343330383e-01,	6.593098044395447e-01,	7.970626354217529e-01,
//        1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,
//	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,
//    1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,
//    1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,
//    1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,
//    1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,
//    1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,
//    1.000000000000000e+00,	1.000000000000000e+00	,1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,
//    1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,
//    1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00,	1.000000000000000e+00});

    std::vector<std::vector<float> > sampled(values01.size(), std::vector<float>(values01.front().size()));

    timer::Timer time_cpp11;

    empirical_quantile::ExplicitQuantile<int, float> quant_expl(lb, ub, gridN);
    quant_expl.set_sample(sample_implicit);
    time_cpp11.reset();
    for(size_t i = 0; i != values01.size(); i++)
        quant_expl.transform(values01[i], sampled[i]);
    std::cout << "\ntotal time explicit       : " << std::scientific << time_cpp11.elapsed_seconds() << std::endl;
    std::cout << "time per transform: " << std::scientific << time_cpp11.elapsed_seconds()/double(sampled.size()) << std::endl;

    for(size_t i = 0; i != sampled.size(); i++)
    {
        for(size_t j = 0; j != sampled[i].size(); j++)
        {
            if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
            {
                std::cout << "beyond bounds" << std::endl;
                std::cout << sampled[i][j] << std::endl;
            }
        }
    }
    data_io::write_default2d("maps/sampled_explicit.dat", sampled, 5);

    empirical_quantile::ImplicitQuantile<int, float> quant_impl(lb, ub, gridN);
    quant_impl.set_sample_shared(sample);
    time_cpp11.reset();
    for(size_t i = 0; i != values01.size(); i++)
        quant_impl.transform(values01[i], sampled[i]);
    std::cout << "\ntotal time implicit       : " << std::scientific << time_cpp11.elapsed_seconds() << std::endl;
    std::cout << "time per transform: " << std::scientific << time_cpp11.elapsed_seconds()/double(sampled.size()) << std::endl;

    for(size_t i = 0; i != sampled.size(); i++)
    {
        for(size_t j = 0; j != sampled[i].size(); j++)
        {
            if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
            {
                std::cout << "beyond bounds" << std::endl;
                std::cout << sampled[i][j] << std::endl;
            }
        }
    }
    data_io::write_default2d("maps/sampled_implicit.dat", sampled, 5);

    empirical_quantile::ImplicitQuantileSorted<int, float> quant_impls(lb, ub, gridN);
    quant_impls.set_sample_shared(sample);
    time_cpp11.reset();
    for(size_t i = 0; i != values01.size(); i++)
        quant_impls.transform(values01[i], sampled[i]);
    std::cout << "\ntotal time implicit sorted: " << std::scientific << time_cpp11.elapsed_seconds() << std::endl;
    std::cout << "time per transform: " << std::scientific << time_cpp11.elapsed_seconds()/double(sampled.size()) << std::endl;

    for(size_t i = 0; i != sampled.size(); i++)
    {
        for(size_t j = 0; j != sampled[i].size(); j++)
        {
            if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
            {
                std::cout << "beyond bounds" << std::endl;
                std::cout << sampled[i][j] << std::endl;
            }
        }
    }
    data_io::write_default2d("maps/sampled_implicit_sorted.dat", sampled, 5);
}



std::pair<double, double> test_Nd_time(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nsamples, size_t Nrolls)
{
    std::mt19937_64 generator;
    generator.seed(1 + gridN.size() + gridN.front());
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    typedef trie_based::TrieBased<trie_based::NodeCount<int>,int> sample_type;
    std::shared_ptr<sample_type> sample = std::make_shared<sample_type>();

    std::vector<int> temp(gridN.size());
    for(size_t i = 0; i != Nsamples;)
    {
        for(size_t j = 0; j != gridN.size(); j++)
        {
            temp[j] = static_cast<int>(std::round(ureal01(generator)*(gridN[j] - 1.0)));
        }
        if(!sample->search(temp))
        {
            sample->insert(temp);
            i++;
        }
    }

    std::vector<std::vector<float> > values01;

    std::vector<float> temp1(gridN.size());
    std::vector<float> temp2(temp1.size());

    for(size_t i = 0; i != Nrolls; ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        values01.push_back(temp1);
    }

    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 0.0;
        values01.push_back(temp1);
    }
    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 1.0;
        values01.push_back(temp1);
    }

    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = 0.0;
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        values01.push_back(temp1);
    }
    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = 1.0;
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        values01.push_back(temp1);
    }

    for(size_t i = 0; i != temp1.size() - 1; ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = 0.0;
        }
        values01.push_back(temp1);
    }
    for(size_t i = 0; i != temp1.size() - 1; ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = 1.0;
        }
        values01.push_back(temp1);
    }

    std::vector<std::vector<float> > sampled(values01.size(), std::vector<float>(values01.front().size()));

    timer::Timer time_cpp11;

    empirical_quantile::ImplicitQuantile<int, float> quant_impl(lb, ub, gridN);
    quant_impl.set_sample_shared(sample);
    time_cpp11.reset();
    for(size_t i = 0; i != values01.size(); i++)
        quant_impl.transform(values01[i], sampled[i]);
    double first = time_cpp11.elapsed_seconds()/double(sampled.size());
    for(size_t i = 0; i != sampled.size(); i++)
    {
        for(size_t j = 0; j != sampled[i].size(); j++)
        {
            if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
            {
                std::cout << "beyond bounds" << std::endl;
                std::cout << sampled[i][j] << std::endl;
            }
        }
    }

    empirical_quantile::ImplicitQuantileSorted<int, float> quant_impls(lb, ub, gridN);
    quant_impls.set_sample_shared(sample);
    time_cpp11.reset();
    for(size_t i = 0; i != values01.size(); i++)
        quant_impls.transform(values01[i], sampled[i]);
    double second = time_cpp11.elapsed_seconds()/double(sampled.size());
    for(size_t i = 0; i != sampled.size(); i++)
    {
        for(size_t j = 0; j != sampled[i].size(); j++)
        {
            if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
            {
                std::cout << "beyond bounds" << std::endl;
                std::cout << sampled[i][j] << std::endl;
            }
        }
    }
    return std::make_pair(first, second);
}

void grid_test_Nd()
{
    for(size_t g_size = 1000; g_size < 10000 + 1; g_size+=1000)
    {
//        std::cout << g_size << std::endl;

        size_t N = 100;
        std::vector<size_t> g(N);
        for(size_t i = 0; i != N; i++)
        {
            g[i] = g_size;
        }
        std::vector<float> lb(N, -10);
        std::vector<float> ub(N, 10);
        auto rez = test_Nd_time(g, lb, ub, 500000, 1e5);
        std::cout << g_size << '\t' << std::scientific << rez.first << '\t' << rez.second << std::endl;
    }
}

void dim_test_Nd()
{
    for(size_t dim_size = 10; dim_size < 500 + 1; dim_size+=10)
    {
        size_t N = dim_size;
        std::vector<size_t> g(N, 1000);
        std::vector<float> lb(N, -10);
        std::vector<float> ub(N, 10);
        auto rez = test_Nd_time(g, lb, ub, 100000, 1e3);
        std::cout << dim_size << '\t' << std::scientific << rez.first << '\t' << rez.second << std::endl;
    }
}







std::vector<double> test_2d_time(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nsamples, size_t Nrolls)
{
    std::mt19937_64 generator;
    generator.seed(1 + gridN.size() + gridN.front());
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    typedef trie_based::TrieBased<trie_based::NodeCount<int>,int> sample_type;
    std::shared_ptr<sample_type> sample = std::make_shared<sample_type>();
    
    std::vector<std::vector<int>> sample_int;
    std::vector<int> temp(gridN.size());
    for(size_t i = 0; i != Nsamples;)
    {
        for(size_t j = 0; j != gridN.size(); j++)
        {
            temp[j] = static_cast<int>(std::round(ureal01(generator)*(gridN[j] - 1.0)));
        }
        if(!sample->search(temp))
        {
            sample->insert(temp);
            sample_int.push_back(temp);
            i++;
        }
    }

    std::vector<std::vector<float> > values01;

    std::vector<float> temp1(gridN.size());
    std::vector<float> temp2(temp1.size());

    for(size_t i = 0; i != Nrolls; ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        values01.push_back(temp1);
    }

    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 0.0;
        values01.push_back(temp1);
    }
    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 1.0;
        values01.push_back(temp1);
    }

    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = 0.0;
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        values01.push_back(temp1);
    }
    for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = 1.0;
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        values01.push_back(temp1);
    }

    for(size_t i = 0; i != temp1.size() - 1; ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = 0.0;
        }
        values01.push_back(temp1);
    }
    for(size_t i = 0; i != temp1.size() - 1; ++i)
    {
        for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        for(size_t j = i + 1; j < temp1.size(); j++)
        {
            temp1[j] = 1.0;
        }
        values01.push_back(temp1);
    }

    std::vector<std::vector<float> > sampled(values01.size(), std::vector<float>(values01.front().size()));

    timer::Timer time_all_trans, full_time;
    std::vector<double> result;
    
//    empirical_quantile::ExplicitQuantile<int, float> quant_expl(lb, ub, gridN);
//    quant_expl.set_sample(sample_int);
//    full_time.reset();
//    time_all_trans.reset();
//    for(size_t i = 0; i != values01.size(); i++)
//        quant_expl.transform(values01[i], sampled[i]);
//    result.push_back(time_all_trans.elapsed_seconds());
//    result.push_back(result.back()/double(sampled.size()));
//    result.push_back(full_time.elapsed_seconds());
//    for(size_t i = 0; i != sampled.size(); i++)
//    {
//        for(size_t j = 0; j != sampled[i].size(); j++)
//        {
//            if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
//            {
//                std::cout << "beyond bounds" << std::endl;
//                std::cout << sampled[i][j] << std::endl;
//            }
//        }
//    }

    empirical_quantile::ImplicitQuantile<int, float> quant_impl(lb, ub, gridN);
    full_time.reset();
    quant_impl.set_sample_shared(sample);
    time_all_trans.reset();
    for(size_t i = 0; i != values01.size(); i++)
        quant_impl.transform(values01[i], sampled[i]);
    result.push_back(time_all_trans.elapsed_seconds());
    result.push_back(result.back()/double(sampled.size()));
    result.push_back(full_time.elapsed_seconds());
    for(size_t i = 0; i != sampled.size(); i++)
    {
        for(size_t j = 0; j != sampled[i].size(); j++)
        {
            if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
            {
                std::cout << "beyond bounds" << std::endl;
                std::cout << sampled[i][j] << std::endl;
            }
        }
    }

    empirical_quantile::ImplicitQuantileSorted<int, float> quant_impls(lb, ub, gridN);
    full_time.reset();
    quant_impls.set_sample_shared(sample);
    time_all_trans.reset();
    for(size_t i = 0; i != values01.size(); i++)
        quant_impls.transform(values01[i], sampled[i]);
    result.push_back(time_all_trans.elapsed_seconds());
    result.push_back(result.back()/double(sampled.size()));
    result.push_back(full_time.elapsed_seconds());
    for(size_t i = 0; i != sampled.size(); i++)
    {
        for(size_t j = 0; j != sampled[i].size(); j++)
        {
            if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
            {
                std::cout << "beyond bounds" << std::endl;
                std::cout << sampled[i][j] << std::endl;
            }
        }
    }
    return result;
}



void grid_test_2d()
{
    size_t N = 2;
    std::vector<size_t> g(N, 1440);
    std::vector<float> lb(N, -acos(-1.0));
    std::vector<float> ub(N, acos(-1.0));
    auto rez = test_2d_time(g, lb, ub, 1e6, 1e4);
    for(const auto &i : rez)
        std::cout << std::scientific << i << '\t';
    std::cout << std::endl;
}
