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
#include "test1d.h"
#include "test2d.h"
#include "test3d.h"
#include "data_io.h"
#include "test_kde.h"
#include "kquantile.h"
#include "mveqf.h"

void test_mveqf()
{
    std::vector<std::vector<double>> sample =
    {
        {-2.442222e-01, 1.137655e+00},
        {5.242222e-01, 9.542629e-01},
        {6.066667e-02, 8.547957e-01},
        {-1.788889e-01, 8.454707e-01},
        {-1.788889e-01, 5.750444e-01},
        {4.246667e-01, 5.470693e-01},
        {7.357778e-01, 4.849023e-01},
        {7.171111e-01, 3.854352e-01},
        {8.866667e-02, 3.698934e-01},
        {5.397778e-01, 3.450266e-01},
        {5.864444e-01, 3.201599e-01},
        {-7.917778e-01, 3.108348e-01},
        {-1.221111e+00, 2.828597e-01},
        {-9.940000e-01, 2.486679e-01},
        {-4.868889e-01, 2.082593e-01},
        {-5.428889e-01, 2.020426e-01},
        {-3.748889e-01, 2.020426e-01},
        {3.997778e-01, 1.616341e-01},
        {-2.411111e-01, 1.460924e-01},
        {1.104444e-01, 1.398757e-01},
        {-1.260000e-01, 1.056838e-01},
        {-1.650444e+00, 9.635879e-02},
        {-2.815556e-01, 9.014210e-02},
        {4.822222e-02, 8.703375e-02},
        {1.042222e-01, 4.351687e-02},
        {1.088889e-02, 1.554174e-02},
        {-2.597778e-01, -4.973357e-02},
        {-1.099778e+00, -7.460036e-02},
        {9.193333e-01, -9.635879e-02},
        {1.197778e-01, -1.802842e-01},
        {6.984444e-01, -1.865009e-01},
        {5.802222e-01, -2.238011e-01},
        {-1.851111e-01, -2.300178e-01},
        {9.815556e-01, -2.517762e-01},
        {-4.184444e-01, -2.735346e-01},
        {7.000000e-02, -2.766430e-01},
        {-2.131111e-01, -2.828597e-01},
        {3.873333e-01, -2.828597e-01},
        {7.482222e-01, -3.232682e-01},
        {2.722222e-01, -3.761101e-01},
        {1.011111e-01, -4.258437e-01},
        {-7.700000e-01, -4.320604e-01},
        {-4.200000e-02, -4.320604e-01},
        {-4.666667e-03, -4.973357e-01},
        {1.042222e-01, -5.222025e-01},
        {1.944444e-01, -5.346359e-01},
        {3.157778e-01, -5.563943e-01},
        {5.242222e-01, -5.626110e-01},
        {5.273333e-01, -6.993783e-01},
        {5.397778e-01, -9.325044e-01}
    };
    data_io::write_default2d("maps/sample2d.dat", sample, 5);

    mveqf::MVEQF<int, double> qf;
    qf.set_kernel_type(0);
    std::vector<size_t> empty;
    qf.set_sample_and_bounds(sample, std::vector<double>(2, -2), std::vector<double>(2, 2), empty, 500, 1e-8, 100000);

    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<double> ureal01(0.0,1.0);
    size_t nrolls = 500;
    std::vector<std::vector<double>> sampled(nrolls, std::vector<double>(2));
    for(size_t i = 0; i != nrolls; i++)
    {
        std::vector<double> t = {ureal01(generator), ureal01(generator)};
        sampled[i] = qf.transform(t);
    }
    data_io::write_default2d("maps/sampled2d.dat", sampled, 5);


    nrolls = 50000;
    sampled = std::vector<std::vector<double>>(nrolls, std::vector<double>(2));
    auto vals01 = sampled;
    for(auto & i : vals01)
    {
        for(auto & j : i)
        {
            j = ureal01(generator);
        }
    }

    const auto nthreads = std::thread::hardware_concurrency();
    auto first = sampled.begin();
    auto last = sampled.end();
    const auto size = last - first;
    const auto size_per_thread = size / nthreads;

    std::vector<std::future<void>> futures;
    for(unsigned int i = 0; i < nthreads - 1; i++)
    {
        futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, &vals01, &sampled, &qf]()
        {
            for(auto it = start; it != start + size_per_thread; ++it)
            {
                qf.transform(vals01[std::distance(sampled.begin(), it)], *it);
            }
        }));
    }
    futures.emplace_back(
        std::async([start = first + (nthreads - 1) * size_per_thread, last, &vals01, &sampled, &qf]()
    {
        for(auto it = start; it != last; ++it)
        {
            qf.transform(vals01[std::distance(sampled.begin(), it)], *it);
        }
    }));

    for(auto &&future : futures)
    {
        if(future.valid())
        {
            future.get();
        }
        else
        {
            throw std::runtime_error("Something going wrong.");
        }
    }
    data_io::write_default2d("maps/sampled2d_1m.dat", sampled, 5);
}


void test_mveqf2()
{
    std::vector<std::vector<double>> sample =
    {
        {5.000000e+01,	5.000000e+01},
        {-5.000000e+01,	5.000000e+01},
        {5.000000e+01,	-5.000000e+01},
        {-5.000000e+01,	-5.000000e+01},
        {-3.603057e+01,	7.704624e+01},
        {-3.210936e+01,	-3.876422e+01},
        {-5.673265e+00,	-6.626629e+01},
        {-9.162707e+01,	1.064726e+01},
        {-5.508124e+01,	-9.601210e+01},
        {-9.279164e+01,	-2.018437e-01}
    };
    data_io::write_default2d("maps/sample2d.dat", sample, 5);

//    std::vector<size_t> count = {1, 25, 84, 35, 101, 36, 77, 54, 76, 62};
//    for(size_t i = 0; i != count.size(); i++)
//    {
//        for(size_t j = 0; j < count[i] - 1; j++)
//        {
//            sample.push_back(sample[i]);
//        }
//    }

    mveqf::MVEQF<int, double> qf;
    qf.set_kernel_type(0);
    std::vector<size_t> empty;// = {1, 25, 84, 35, 101, 36, 77, 54, 76, 62};
    qf.set_sample_and_bounds(sample, std::vector<double>(2, -100), std::vector<double>(2, 100), empty, 500, 1e-8, 100000);

    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<double> ureal01(0.0,1.0);
    size_t nrolls = 500;
    std::vector<std::vector<double>> sampled(nrolls, std::vector<double>(2));
    for(size_t i = 0; i != nrolls; i++)
    {
        std::vector<double> t = {ureal01(generator), ureal01(generator)};
        sampled[i] = qf.transform(t);
    }
    data_io::write_default2d("maps/sampled2d.dat", sampled, 5);


    nrolls = 50000;
    sampled = std::vector<std::vector<double>>(nrolls, std::vector<double>(2));
    auto vals01 = sampled;
    for(auto & i : vals01)
    {
        for(auto & j : i)
        {
            j = ureal01(generator);
        }
    }

    const auto nthreads = std::thread::hardware_concurrency();
    auto first = sampled.begin();
    auto last = sampled.end();
    const auto size = last - first;
    const auto size_per_thread = size / nthreads;

    std::vector<std::future<void>> futures;
    for(unsigned int i = 0; i < nthreads - 1; i++)
    {
        futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, &vals01, &sampled, &qf]()
        {
            for(auto it = start; it != start + size_per_thread; ++it)
            {
                qf.transform(vals01[std::distance(sampled.begin(), it)], *it);
            }
        }));
    }
    futures.emplace_back(
        std::async([start = first + (nthreads - 1) * size_per_thread, last, &vals01, &sampled, &qf]()
    {
        for(auto it = start; it != last; ++it)
        {
            qf.transform(vals01[std::distance(sampled.begin(), it)], *it);
        }
    }));

    for(auto &&future : futures)
    {
        if(future.valid())
        {
            future.get();
        }
        else
        {
            throw std::runtime_error("Something going wrong.");
        }
    }
    data_io::write_default2d("maps/sampled2d_1m.dat", sampled, 5);
}

int main()
{
    timer::Timer time;
//
//    std::cout << tt.front() == tt.back() << std::endl;

//    simple_empirical_1d();
//    simple1d_example();

//    test_1d1();
//    test_1d2();
//    test_1d3();
//    test_1d4();
//    test_1d5();


//    test_2d1();
//    test_2d2();
//    test_3d1();
//    test_3d2();

//    test_1d_func();
//    test_2d_func();

//    kquantile::test_1d_kquantile();
//    test_2d_kquantile<double>();
//    test_3d_kquantile<double>();

//    test_1d_uniform_vs_nonuniform();
//    test_2d_uniform_vs_nonuniform();

//    test_grid_10d();

    /// 2-dimensional test

//    std::vector<size_t> g = {13, 20};
//    //std::vector<size_t> g(N, 10);
//    std::vector<float> lb = {-3, -2};
//    std::vector<float> ub = {1, 4};
//    test_Nd(g, lb, ub, 100, 1e5);

    /// N-dimensional test

//    std::mt19937_64 generator;
//    generator.seed(1);
//
//    std::uniform_int_distribution<size_t> dis(2, 20);
//    size_t N = 10;
//    std::vector<size_t> g(N);
//    for(size_t i = 0; i != N; i++)
//    {
//        g[i] = dis(generator);
//        std::cout << g[i] << '\t';
//    }
//    std::cout << std::endl;
//    std::vector<float> lb(N, -10);
//    std::vector<float> ub(N, 10);
//    test_Nd(g, lb, ub, 50000, 1e2);


//    size_t N = 2;
//    std::vector<size_t> g(N);
//    for(size_t i = 0; i != N; i++)
//    {
//        g[i] = 10000;
//    }
//    std::cout << std::endl;
//    std::vector<float> lb(N, -10);
//    std::vector<float> ub(N, 10);
//    test_Nd(g, lb, ub, 10000, 1e5);

    /// grid and dim test
//    grid_test_Nd();
//    dim_test_Nd();


//    size_t N = 1;
//    std::vector<size_t> g(N);
//    for(size_t i = 0; i != N; i++)
//    {
//        g[i] = 100;
//    }
//    std::cout << std::endl;
//    std::vector<float> lb(N, -10);
//    std::vector<float> ub(N, 10);
//    test_Nd(g, lb, ub, 50, 1e5);

//    grid_test_2d();
//    sample_size_test(10);

//    worst_space_test_dim();
//    worst_space_test_grid();
//    worst_space_test_grid_2d();

//    worst_space_check();
    //std::vector<std::vector<double>> x(100,std::vector<double>(10));


    /// kde
//    mveqf::kde::test1d();
//    mveqf::kde::test1d_1();
//    mveqf::kde::test2d();
//    mveqf::kde::test2d_2();

    test_mveqf();

    ///
    std::cout << time.elapsed_seconds() << std::endl;
}

// // trie based layer test
//    trie_based::TrieBased<trie_based::NodeCount<int>,int>  sample;
//    sample.set_dimension(3);
//    std::vector<std::vector<int>> in_sample;
//    in_sample =
//    {
//        {1,0,0}, /* baa */
//        {2,0,0}, /* caa */
//        {4,0,0}, /* eaa */
//        {0,2,0}, /* aca */
//        {4,4,0}, /* eea */
//        {4,3,0}, /* eda */
//        {3,3,0}, /* dda */
//
//        {0,0,1}, /* aab */
//
//        {3,0,2}, /* dac */
//        {0,3,2}, /* adc */
//
//        {0,3,3}, /* add */
//
//        {2,0,4}, /* cae */
//        {2,1,4}, /* cbe */
//        {2,2,4} /* cce */
//    };
//    for(const auto & i : in_sample)
//        sample.insert(i);
//    sample.fill_tree_count();
//    //sort();
//    auto rez = sample.get_layer_count();
//    std::cout << rez.size() << std::endl;
//    std::cout << sample.last_layer.size() << std::endl;
//
//    std::cout << std::endl;
//
//    for(const auto &i : rez)
//    {
//        std::cout << i.first << std::endl;
//        for(const auto &j : i.second)
//            std::cout << j << ' ';
//        std::cout << std::endl;
//    }
