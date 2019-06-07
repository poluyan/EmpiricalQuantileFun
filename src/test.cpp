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
#include <chrono>
#include <thread>

std::vector<std::vector<double>> centers = { {3, 3}, {-5, 0}, {0, -5} };
double threeExp(double x, double y)
{
    double rez = 0;
    rez += std::exp(-(std::pow(x - centers[0][0], 2.0) + std::pow(y - centers[0][1], 2.0))*0.75);
    rez += std::exp(-(std::pow(x - centers[1][0], 2.0) + std::pow(y - centers[1][1], 2.0))*0.5)*0.75;
    rez += std::exp(-(std::pow(x - centers[2][0], 2.0) + std::pow(y - centers[2][1], 2.0))*0.25)*0.5;
    return rez;
}

double threeExp1d(double x)
{
    double rez = 0;
    rez += std::exp(-(std::pow(x - centers[0][0], 2.0))*0.75);
    rez += std::exp(-(std::pow(x - centers[1][0], 2.0))*0.5)*0.75;
    rez += std::exp(-(std::pow(x - centers[2][0], 2.0))*0.25)*0.5;
    return rez;
}

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
    quant.set_sample_and_fill_count(sample);

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
    quant.set_sample_and_fill_count(sample);

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

void implicit_quantile_graph_sorted(float lb, float ub, std::vector<size_t> gridn, size_t nrolls)
{
    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    std::vector<std::vector<float> > sampled;
    std::vector<std::vector<float> > values01;

    empirical_quantile::GraphQuantile <int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);

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
    data_io::write_default2d("maps/sampled_implicit_graph_sorted.dat", sampled, 5);
}


void test_1d1()
{
    std::vector<size_t> grid_number = {6};

    std::vector<std::vector<int>> sample_implicit = {{0},{1},{2},{3},{4},{5}};

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(-2, 4, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class(-2, 4, grid_number, sample_implicit, 1e3);
    implicit_quantile_class_sorted(-2, 4, grid_number, sample_implicit, 1e3);
}

void test_1d5()
{
    std::vector<size_t> grid_number = {6};

    std::vector<std::vector<int>> sample_implicit = {{3}};

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

    std::vector<std::vector<int>> sample_implicit = {{4},{1},{3}};
    //sample_implicit.push_back(std::vector{8});

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(-2, 4, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class(-2, 4, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class_sorted(-2, 4, grid_number, sample_implicit, 1e+3);
}


void test_2d1()
{
    std::vector<size_t> grid_number = {9, 10};

    std::vector<std::vector<int>> sample_implicit =
    {
        {2,6},

        {3,2},
        {3,3},
        {3,5},
        {3,6},
        {3,7},

        {4,5},
        {4,6},
        {4,7},

        {5,3},
        {5,4},
        {5,5},
        {5,6},
        {5,7},

        {6,3},
        {6,4}
    };

/// multivariate quantile function [0,1]^n -> [-3,3]^n nrolls = 2e3
    explicit_quantile(-3, 3, grid_number, sample_implicit, 2e3);
    implicit_quantile_class(-3, 3, grid_number, sample_implicit, 2e3);
    implicit_quantile_class_sorted(-3, 3, grid_number, sample_implicit, 2e3);
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

void test_1d_func()
{
    std::vector<size_t> gridn = {10};
    std::vector<std::vector<int>> sample_implicit;


//    auto sample = std::make_shared<trie_based::TrieBased<trie_based::NodeCount<int>,int>>();
//    sample->set_dimension(gridn.size());
    std::vector<std::vector<int>> sample;

    std::vector<double> dx(gridn.size());
    std::vector<std::vector<double>> grids(gridn.size());
    for(size_t i = 0; i != grids.size(); i++)
    {
        std::vector<double> grid(gridn[i] + 1);
        double startp = -10.0;
        double endp = 10.0;
        double es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/double(gridn[i]);
        }
        grids[i] = grid;
        dx[i] = es/(double(gridn[i])*2);
    }

    for(size_t i = 0; i != gridn[0]; i++)
    {
        //std::cout << grids[0][i] + dx[0] << '\t' << grids[1][j] + dx[1] << std::endl;
        size_t v = size_t(100.0*threeExp1d(grids[0][i] + dx[0]));
        if(v)
        {
            std::vector<int> a = {int(i)};
            //std::cout << "add " << i << '\t' << j << std::endl;
//                sample->insert(a, v);
            for(size_t k = 0; k != v; k++)
                sample.push_back(a);
        }
    }
    /*std::vector<int> a = {5,5};
    //sample->insert(a, 4);
    sample.push_back(a);
    sample.push_back(a);
    sample.push_back(a);
    sample.push_back(a);
    sample.push_back(a);

    a[0] = 4;
    a[1] = 5;
    //sample->insert(a, 2);
    sample.push_back(a);
    sample.push_back(a);
    sample.push_back(a);

    a[0] = 5;
    a[1] = 4;
    //sample->insert(a, 2);
    sample.push_back(a);*/



//    std::cout << sample->root->count << std::endl;
    std::cout << sample.size() << std::endl;


    size_t nrolls = 1e+5;

    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    std::vector<std::vector<float> > sampled;
    std::vector<std::vector<float> > values01;

    float lb = -10.0, ub = 10.0;

//    empirical_quantile::ImplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
//    quant.set_sample_shared_and_fill_count(sample);
//    quant.set_sample_shared(sample);

    empirical_quantile::ExplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
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
    data_io::write_default2d("maps/sampled_implicit.dat", sampled, 5);
}

void test_2d_func()
{
    std::vector<size_t> gridn = {100, 100};
//    std::vector<size_t> gridn = {50, 50};

//    auto sample = std::make_shared<trie_based::Trie<trie_based::NodeCount<int>,int>>();
    auto sample = std::make_shared<trie_based::TrieBased<trie_based::NodeCount<int>,int>>();
    sample->set_dimension(gridn.size());
//    std::vector<std::vector<int>> sample;

    std::vector<double> dx(gridn.size());
    std::vector<std::vector<double>> grids(gridn.size());
    for(size_t i = 0; i != grids.size(); i++)
    {
        std::vector<double> grid(gridn[i] + 1);
        double startp = -2.0;
        double endp = 2.0;
        double es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/double(gridn[i]);
        }
        grids[i] = grid;
        dx[i] = es/(double(gridn[i])*2);
    }

    for(size_t i = 0; i != gridn[0]; i++)
    {
        for(size_t j = 0; j != gridn[1]; j++)
        {
            //std::cout << grids[0][i] + dx[0] << '\t' << grids[1][j] + dx[1] << std::endl;
            size_t v = size_t(1000.0*threeExp(grids[0][i] + dx[0], grids[1][j] + dx[1]));
            if(v)
            {
                std::vector<int> a = {int(i), int(j)};
                //std::cout << "add " << i << '\t' << j << std::endl;
                sample->insert(a);
//                sample->insert(a, v);
//                for(size_t k = 0; k != v; k++)
//                    sample.push_back(a);
            }
        }
    }

    /*std::vector<int> a = {5,5};
    //sample->insert(a, 4);
    sample.push_back(a);
    sample.push_back(a);
    sample.push_back(a);
    sample.push_back(a);
    sample.push_back(a);

    a[0] = 4;
    a[1] = 5;
    //sample->insert(a, 2);
    sample.push_back(a);
    sample.push_back(a);
    sample.push_back(a);

    a[0] = 5;
    a[1] = 4;
    //sample->insert(a, 2);
    sample.push_back(a);*/


    size_t nrolls = 1e+5;

    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    std::vector<std::vector<float> > sampled;
    std::vector<std::vector<float> > values01;

    float lb = -2.0, ub = 2.0;

    empirical_quantile::ImplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
//    empirical_quantile::ImplicitTrieQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
    quant.set_sample_shared_and_fill_count(sample);
//    quant.set_sample_shared(sample);

//    empirical_quantile::ExplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
//    quant.set_sample(sample);

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
    data_io::write_default2d("maps/sampled_implicit.dat", sampled, 5);
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
    sample_implicit =
    {
        {1,0,0}, /* baa */
        {2,0,0}, /* caa */
        {4,0,0}, /* eaa */
        {0,2,0}, /* aca */
        {4,4,0}, /* eea */
        {4,3,0}, /* eda */
        {3,3,0}, /* dda */

        {0,0,1}, /* aab */

        {3,0,2}, /* dac */
        {0,3,2}, /* adc */

        {0,3,3}, /* add */

        {2,0,4}, /* cae */
        {2,1,4}, /* cbe */
        {2,2,4} /* cce */
    };

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(0, 5, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class(0, 5, grid_number, sample_implicit, 1e+3);
    implicit_quantile_class_sorted(0, 5, grid_number, sample_implicit, 1e+3);
    implicit_quantile_graph_sorted(0, 5, grid_number, 1e+3);
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
    sample_implicit =
    {
        {0,0,0}, /* aaa */
        {0,0,0}, /* aab */
        {0,0,0}, /* aac */

        {0,2,0}, /* aba */
        {4,4,0}, /* abb */
        {4,3,0}, /* abc */

        {3,3,0}, /* aca */
        {3,3,0}, /* acb */
        {3,3,0}, /* acc */

        {3,3,0}, /* baa */
        {3,3,0}, /* bab */
        {3,3,0}, /* bac */

        {3,3,0}, /* bba */
        {3,3,0}, /* bbb */
        {3,3,0}, /* bbc */

        {3,3,0}, /* bca */
        {3,3,0}, /* bcb */
        {3,3,0}, /* bcc */

        {3,3,0}, /* caa */
        {3,3,0}, /* cab */
        {0,0,1}, /* cac */

        {0,0,1}, /* cba */
        {0,0,1}, /* cbb */
        {0,0,1}, /* cbc */

        {0,0,1}, /* cca */
        {0,0,1}, /* ccb */
        {0,0,1} /* ccc */
    };
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

    std::vector<std::vector<int>> sample_implicit =
    {
        {0},{2},{4},{11},{12},{16},{19},{20},{24},{25},{26},{27},{28},{30},{31},{33},{34},{36},{39},
        {41},{44},{46},{47},{48},{50},{51},{54},{56},{61},{62},{63},{64},{66},{68},{71},{77},{78},{79},{80},
        {81},{82},{84},{85},{87},{88},{93},{94},{95},{96},{99},{101},{103},{104},{106},{107},{110},{112},
        {114},{118},{121},{123},{126},{129},{130},{133},{137},{138},{140},{143},{146},{150},{151},{153},
        {154},{157},{162},{163},{164},{166}
    };

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

    std::vector<std::vector<int>> sample_implicit =
    {
        {0},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{15},{16},{22},{24},{25},{29},{32},{33},{35},
        {38},{39},{41},{43},{48},{51},{54},{55},{56},{58},{59},{62},{65},{69},{73},{76},{77},{78},{80},
        {81},{82},{83},{84},{85},{89},{93},{96},{98},{99},{101},{103},{104},{105},{107},{109},{112},{115},
        {117},{119},{121},{124},{126},{127},{129},{137},{138},{139},{142},{143},{144},{146},{147},{149},{153},
        {154},{155},{156},{158},{162},{163},{164},{165},{166},{171},{172},{174},{176},{177},{179},{184},{185},
        {189},{190},{193},{194},{195},{197},{200},{201},{202},{204},{208},{209},{210},{211},{213},{215},{216},
        {218},{222},{223},{224},{225},{226},{228},{231},{234},{239},{240},{246},{247},{249},{251},{252},{254},
        {255},{257},{258},{259},{262},{268},{269},{271},{275},{277},{280},{281},{283},{284},{286},{287},{289},
        {290},{292},{293},{298},{299},{300},{304},{306},{311},{315},{320},{324},{326},{327},{331},{332},{334},
        {337},{338},{341},{342},{344},{346},{347},{348},{349},{352},{353},{354},{356},{357},{358},{359},{360},
        {361},{364},{366},{367},{369},{370},{373},{376},{378},{379},{380},{381},{387},{388},{389},{391},{393},
        {394},{395},{396},{397},{398},{399}
    };

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
    quant_impl.set_sample_shared_and_fill_count(sample);
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
    quant_impls.set_sample_shared_and_fill_count(sample);
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
    quant_impl.set_sample_shared_and_fill_count(sample);
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
    quant_impls.set_sample_shared_and_fill_count(sample);
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
    quant_impl.set_sample_shared_and_fill_count(sample);
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
    quant_impls.set_sample_shared_and_fill_count(sample);
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

std::vector<double> sample_size_procedure(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nsamples, size_t Nrolls, size_t seed_append)
{
    std::mt19937_64 generator;
    generator.seed(1 + gridN.size() + gridN.front() + Nsamples + seed_append);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    typedef trie_based::TrieBased<trie_based::NodeCount<std::uint8_t>,std::uint8_t> sample_type;
    std::shared_ptr<sample_type> sample = std::make_shared<sample_type>();

    std::vector<std::vector<std::uint8_t>> sample_int;
    std::vector<std::uint8_t> temp(gridN.size());
    for(size_t i = 0; i != Nsamples;)
    {
        for(size_t j = 0; j != gridN.size(); j++)
        {
            temp[j] = static_cast<std::uint8_t>(std::round(ureal01(generator)*(gridN[j] - 1.0)));
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

    std::vector<std::vector<float> > sampled(values01.size(), std::vector<float>(values01.front().size()));

    timer::Timer time_all_trans, full_time;
    std::vector<double> result;

    empirical_quantile::ExplicitQuantile<std::uint8_t, float> quant_expl(lb, ub, gridN);
    quant_expl.set_sample(sample_int);
    full_time.reset();
    time_all_trans.reset();
    for(size_t i = 0; i != values01.size(); i++)
        quant_expl.transform(values01[i], sampled[i]);
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

    empirical_quantile::ImplicitQuantile<std::uint8_t, float> quant_impl(lb, ub, gridN);
    full_time.reset();
    quant_impl.set_sample_shared_and_fill_count(sample);
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

    empirical_quantile::ImplicitQuantileSorted<std::uint8_t, float> quant_impls(lb, ub, gridN);
    full_time.reset();
    quant_impls.set_sample_shared_and_fill_count(sample);
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

void sample_size_test(size_t dim)
{
    size_t tries = 30;
    for(size_t sample_size = 1e5; sample_size < 1e6 + 1; sample_size+=1e5)
    {
        std::cout << sample_size << '\t';
        std::vector<double> times;
        for(size_t i = 0; i != tries; i++)
        {
            std::vector<size_t> g(dim, 100);
            std::vector<float> lb(dim, -1.0);
            std::vector<float> ub(dim, 1.0);
            auto rez = sample_size_procedure(g, lb, ub, sample_size, 1e1, i);
            times.resize(rez.size());
            for(size_t j = 0; j != times.size(); j++)
            {
                times[j] += rez[j];
            }
        }
        for(size_t j = 0; j != times.size(); j++)
        {
            times[j] /= double(tries);
        }
        for(const auto &i : times)
            std::cout << std::scientific << i << '\t';
        std::cout << std::endl;
    }
}






///


std::vector<double> worst_space(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nrolls, size_t seed_append)
{
    std::mt19937_64 generator;
    generator.seed(1 + gridN.size() + gridN.front() + seed_append);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    typedef trie_based::TrieBased<trie_based::NodeCount<std::uint8_t>,std::uint8_t> sample_type;
    std::shared_ptr<sample_type> sample = std::make_shared<sample_type>();

    std::vector<std::vector<std::uint8_t>> variable_values(gridN.size());
    for(size_t i = 0; i != gridN.size(); i++)
    {
        variable_values[i].resize(gridN[i]);
        for(size_t j = 0; j != gridN[i]; j++)
        {
            variable_values[i][j] = j;
        }
    }
    std::vector<std::vector<std::uint8_t>> sample_int = iterate(variable_values);

    for(size_t i = 0; i != sample_int.size(); i++)
    {
        sample->insert(sample_int[i]);
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

    std::vector<std::vector<float> > sampled(values01.size(), std::vector<float>(values01.front().size()));

    timer::Timer time_all_trans, full_time;
    std::vector<double> result;

    empirical_quantile::ExplicitQuantile<std::uint8_t, float> quant_expl(lb, ub, gridN);
    quant_expl.set_sample(sample_int);
    full_time.reset();
    time_all_trans.reset();
    for(size_t i = 0; i != 10; i++)
        quant_expl.transform(values01[i], sampled[i]);
    result.push_back(time_all_trans.elapsed_seconds());
    result.push_back(result.back()/double(10));
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

    empirical_quantile::ImplicitQuantile<std::uint8_t, float> quant_impl(lb, ub, gridN);
    full_time.reset();
    quant_impl.set_sample_shared_and_fill_count(sample);
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

    empirical_quantile::ImplicitQuantileSorted<std::uint8_t, float> quant_impls(lb, ub, gridN);
    full_time.reset();
    quant_impls.set_sample_shared_and_fill_count(sample);
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


void worst_space_test_dim()
{
    size_t tries = 50, nrolls = 1e5, grid_size = 10;
    for(size_t dim = 1; dim != 7; dim++)
    {
        std::cout << dim << '\t';
        std::vector<double> times;
        for(size_t i = 0; i != tries; i++)
        {
            std::vector<size_t> g(dim, grid_size);
            std::vector<float> lb(dim, -1.0);
            std::vector<float> ub(dim, 1.0);
            auto rez = worst_space(g, lb, ub, nrolls, i);

            times.resize(rez.size());
            for(size_t j = 0; j != times.size(); j++)
            {
                times[j] += rez[j];
            }
        }
        for(size_t j = 0; j != times.size(); j++)
        {
            times[j] /= double(tries);
        }
        for(const auto &i : times)
            std::cout << std::scientific << i << '\t';
        std::vector<size_t> g(dim, grid_size);
        std::cout << std::accumulate(g.begin(), g.end(), 1, std::multiplies<size_t>());
        std::cout << '\t' << grid_size << '^' << dim << std::endl;
    }
}

void worst_space_test_grid()
{
    size_t dim = 4, tries = 50, nrolls = 1e5;
    for(size_t grid_size = 1; grid_size != 33; grid_size++)
    {
        std::cout << grid_size << '\t';
        std::vector<double> times;
        for(size_t i = 0; i != tries; i++)
        {
            std::vector<size_t> g(dim, grid_size);
            std::vector<float> lb(dim, -1.0);
            std::vector<float> ub(dim, 1.0);
            auto rez = worst_space(g, lb, ub, nrolls, i);

            times.resize(rez.size());
            for(size_t j = 0; j != times.size(); j++)
            {
                times[j] += rez[j];
            }
        }
        for(size_t j = 0; j != times.size(); j++)
        {
            times[j] /= double(tries);
        }
        for(const auto &i : times)
            std::cout << std::scientific << i << '\t';
        std::vector<size_t> g(dim, grid_size);
        std::cout << std::accumulate(g.begin(), g.end(), 1, std::multiplies<size_t>());
        std::cout << '\t' << grid_size << '^' << dim << std::endl;
    }
}

void worst_space_test_grid_2d()
{
    size_t dim = 2, tries = 50, nrolls = 1e6;
    for(size_t grid_size = 1; grid_size != 100 + 1; grid_size++)
    {
        std::cout << grid_size << '\t';
        std::vector<double> times;
        for(size_t i = 0; i != tries; i++)
        {
            std::vector<size_t> g(dim, grid_size);
            std::vector<float> lb(dim, -1.0);
            std::vector<float> ub(dim, 1.0);
            auto rez = worst_space(g, lb, ub, nrolls, i);

            times.resize(rez.size());
            for(size_t j = 0; j != times.size(); j++)
            {
                times[j] += rez[j];
            }
        }
        for(size_t j = 0; j != times.size(); j++)
        {
            times[j] /= double(tries);
        }
        for(const auto &i : times)
            std::cout << std::scientific << i << '\t';
        std::vector<size_t> g(dim, grid_size);
        std::cout << std::accumulate(g.begin(), g.end(), 1, std::multiplies<size_t>());
        std::cout << '\t' << grid_size << '^' << dim << std::endl;
    }
}





void worst_space_f(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, bool impl)
{
    std::vector<std::vector<std::uint8_t>> variable_values(gridN.size());
    for(size_t i = 0; i != gridN.size(); i++)
    {
        variable_values[i].resize(gridN[i]);
        for(size_t j = 0; j != gridN[i]; j++)
        {
            variable_values[i][j] = j;
        }
    }

    if(impl)
    {
        typedef trie_based::TrieBased<trie_based::NodeCount<std::uint8_t>,std::uint8_t> sample_type;
        std::shared_ptr<sample_type> sample = std::make_shared<sample_type>();

        iterate_trie(variable_values, sample);

        empirical_quantile::ImplicitQuantile<std::uint8_t, float> quant_impl(lb, ub, gridN);
        quant_impl.set_sample_shared_and_fill_count(sample);

        std::cout << "hit it" << std::endl;
        std::chrono::seconds dura(1000);
        std::this_thread::sleep_for(dura);
    }
    else
    {
        std::vector<std::vector<std::uint8_t>> sample_int = iterate(variable_values);
        empirical_quantile::ExplicitQuantile<std::uint8_t, float> quant_expl(lb, ub, gridN);
//        quant_expl.set_sample(sample_int);
//        sample_int.clear();
//        sample_int.shrink_to_fit();

        std::cout << "hit it" << std::endl;
        std::chrono::seconds dura(1000);
        std::this_thread::sleep_for(dura);
    }

//    empirical_quantile::ImplicitQuantileSorted<std::uint8_t, float> quant_impls(lb, ub, gridN);
//    quant_impls.set_sample_shared(sample);
}

void worst_space_d(std::vector<size_t> gridN, std::vector<double> lb, std::vector<double> ub, bool impl)
{
    std::vector<std::vector<int>> variable_values(gridN.size());
    for(size_t i = 0; i != gridN.size(); i++)
    {
        variable_values[i].resize(gridN[i]);
        for(size_t j = 0; j != gridN[i]; j++)
        {
            variable_values[i][j] = j;
        }
    }

    if(impl)
    {
        typedef trie_based::TrieBased<trie_based::NodeCount<int>,int> sample_type;
        std::shared_ptr<sample_type> sample = std::make_shared<sample_type>();

        iterate_trie(variable_values, sample);

        empirical_quantile::ImplicitQuantile<int, double> quant_impl(lb, ub, gridN);
        quant_impl.set_sample_shared_and_fill_count(sample);

        std::cout << "hit it" << std::endl;
        std::chrono::seconds dura(1000);
        std::this_thread::sleep_for(dura);
    }
    else
    {
        std::vector<std::vector<int>> sample_int = iterate(variable_values);
        empirical_quantile::ExplicitQuantile<int, double> quant_expl(lb, ub, gridN);
//        quant_expl.set_sample(sample_int);
//        sample_int.clear();
//        sample_int.shrink_to_fit();

        std::cout << "hit it" << std::endl;
        std::chrono::seconds dura(1000);
        std::this_thread::sleep_for(dura);
    }

//    empirical_quantile::ImplicitQuantileSorted<std::uint8_t, float> quant_impls(lb, ub, gridN);
//    quant_impls.set_sample_shared(sample);
}


void worst_space_check()
{
    size_t dim = 8, grid_size = 8;

    std::vector<size_t> g(dim, grid_size);
    std::vector<double> lb(dim, -1.0);
    std::vector<double> ub(dim, 1.0);
    worst_space_d(g, lb, ub, true);
}
