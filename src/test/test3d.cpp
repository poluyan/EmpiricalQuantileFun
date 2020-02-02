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
#include "test3d.h"
#include <vector>

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
    explicit_quantile(0, 5, grid_number, sample_implicit, 1e+4);
    implicit_quantile_class(0, 5, grid_number, sample_implicit, 1e+4);
    implicit_quantile_class_sorted(0, 5, grid_number, sample_implicit, 1e+4);
    implicit_quantile_class_sorted_interp(0, 5, grid_number, sample_implicit, 1e+4);
//    implicit_quantile_graph_sorted(0, 5, grid_number, 1e+3);
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
        {0,0,1}, /* aab */
        {0,0,2}, /* aac */

        {0,1,0}, /* aba */
        {0,1,1}, /* abb */
        {0,1,2}, /* abc */

        {0,2,0}, /* aca */
        {0,2,1}, /* acb */
        {0,2,2}, /* acc */

        {1,0,0}, /* baa */
        {1,0,1}, /* bab */
        {1,0,2}, /* bac */

        {1,1,0}, /* bba */
        {1,1,1}, /* bbb */
        {1,1,2}, /* bbc */

        {1,2,0}, /* bca */
        {1,2,1}, /* bcb */
        {1,2,2}, /* bcc */

        {2,0,0}, /* caa */
        {2,0,1}, /* cab */
        {2,0,2}, /* cac */

        {2,1,0}, /* cba */
        {2,1,1}, /* cbb */
        {2,1,2}, /* cbc */

        {2,2,0}, /* cca */
        {2,2,1}, /* ccb */
        {2,2,2} /* ccc */
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
    implicit_quantile_class_sorted_interp(0, 3, grid_number, sample_implicit, 1e+3);
}