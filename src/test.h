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

#ifndef TEST_H
#define TEST_H

#include "timer.h"
#include "trie_based.h"
#include "quantile.h"

void explicit_quantile(float lb, float ub, std::vector<size_t> gridn, std::vector<std::vector<int> > &sample, size_t nrolls);
void implicit_quantile_class(float lb, float ub, std::vector<size_t> gridn, std::vector<std::vector<int> > &sample, size_t nrolls);
void implicit_quantile_class_sorted(float lb, float ub, std::vector<size_t> gridn, std::vector<std::vector<int> > &sample, size_t nrolls);

void test_1d1();
void test_1d2();
void test_1d3();
void test_1d4();
void test_1d5();

void test_2d1();
void test_2d2();

void test_3d1();
void test_3d2();

void test_grid_10d();

void test_Nd(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nsamples, size_t Nrolls);


std::pair<double, double> test_Nd_time(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nsamples, size_t Nrolls);
void grid_test_Nd();
void dim_test_Nd();

void grid_test_2d();

void sample_size_test(size_t dim);



template <typename T>
bool increase(const std::vector<std::vector<T>> &v, std::vector<size_t> &it)
{
    for(size_t i = 0, size = it.size(); i != size; i++)
    {
        const size_t index = size - 1 - i;
        ++it[index];
        if(it[index] == v[index].size())
        {
            it[index] = 0;
        }
        else
        {
            return true;
        }
    }
    return false;
}

template <typename T>
std::vector<T> get_line(const std::vector<std::vector<T>> &v, std::vector<size_t> &it)
{
    std::vector<T> rez(v.size());
    for(size_t i = 0, size = v.size(); i != size; i++)
    {
        rez[i] = v[i][it[i]];
    }
    return rez;
}

template <typename T>
std::vector<std::vector<T>> iterate(const std::vector<std::vector<T>> &v)
{
    std::vector<size_t> it(v.size(), 0);
    std::vector<std::vector<T>> values;
    do
    {
        values.push_back(get_line(v, it));
    }
    while(increase(v, it));
    return values;
}

void worst_space_test_dim();
void worst_space_test_grid();

#endif
