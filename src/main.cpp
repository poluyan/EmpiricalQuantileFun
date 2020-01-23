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
#include "kde.h"
#include "kquantile.h"

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
    mveqf::test_2d_kquantile<double>();
//    mveqf::test_3d_kquantile<float>();
    
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
//    kde::test1d();
//    kde::test2d();
    
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
