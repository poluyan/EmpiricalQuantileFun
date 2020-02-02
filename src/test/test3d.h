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
#ifndef TEST3D_H
#define TEST3D_H


#include "../kquantile.h"

void test_3d1();
void test_3d2();


template <typename T>
void test_3d_kquantile()
{
    
    const size_t kt = 4;//4
    const size_t dim = 3;
    const size_t N = 500;
    const size_t nrolls = 10000;
    const T threshold = 1e-11;
    const size_t multi = 1e+9;


    std::vector<size_t> gridn(dim, N);

    auto sample = std::make_shared<std::vector<std::vector<T>>>();
    std::vector<std::vector<T>> ss;

    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<T> ureal01(0.0,1.0);
    /*std::normal_distribution<T> norm1(0.0,0.5);
    std::normal_distribution<T> norm2(0.0,1.0);

    std::vector<T> p1 = {2.5, 2.5, 2.5};
    std::vector<T> p2 = {7.5, 7.5, 7.5};

    for(size_t i = 0; i != 1e4; i++)
    {
        if(ureal01(generator) < 0.5)
        {
            std::vector<T> t = {p1[0] + norm1(generator), p1[1] + norm1(generator), p1[2] + norm1(generator)};
            sample->push_back(t);
            ss.push_back(t);
        }
        else
        {
            std::vector<T> t = {p2[0] + norm2(generator), p2[1] + norm2(generator), p1[2] + norm2(generator)};
            sample->push_back(t);
            ss.push_back(t);
        }
    }*/

    std::ifstream gridIn;
    gridIn.open("maps/Trie/3d/123.dat");
    if(!gridIn.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }
    std::vector<T> lbv(dim, std::numeric_limits<T>::max());
    std::vector<T> ubv(dim, std::numeric_limits<T>::min());
    std::vector<T> t(dim);
    while(!gridIn.eof())
    {
        gridIn >> t[0];
        gridIn >> t[1];
        gridIn >> t[2];

        for(size_t i = 0; i != t.size(); i++)
        {
            if(t[i] < lbv[i])
                lbv[i] = t[i];
            if(t[i] > ubv[i])
                ubv[i] = t[i];
        }

        sample->push_back(t);
        ss.push_back(t);
        if(sample->size() > nrolls)
            break;
    }
//    sample->pop_back();
//    ss.pop_back();
    gridIn.close();
//
    std::cout << lbv[0] << '\t' << lbv[1] << '\t' << lbv[2] << std::endl;
    std::cout << ubv[0] << '\t' << ubv[1] << '\t' << ubv[2] << std::endl;


    mveqf::kde::KDE<T> obj1;
    obj1.set_dimension(dim);
    obj1.set_kernel_type(kt);
    obj1.set_sample_shared(sample);
    std::vector<std::vector<T>> output;
    for(size_t i = 0; i != ss.size(); i++)
    {
        auto t = ss[i];
        t.push_back(obj1.pdf(t));
        output.push_back(t);
    }
    data_io::write_default2d("maps/Trie/3d/sampled3d_init.dat", output, 5);



    mveqf::kde::KDE<T> obj;
    obj.set_dimension(dim);
    obj.set_kernel_type(kt);
    obj.set_sample_shared(sample);


    ///
    std::vector<std::vector<T>> grids(gridn.size());
    std::vector<T> dx(gridn.size());

    const T lb = -15, ub = 15;
    for(size_t i = 0; i != grids.size(); i++)
    {
        size_t num_points = gridn[i];
        std::vector<T> onegrid(num_points);
        T es = ub - lb;
        for(size_t j = 0; j != onegrid.size(); j++)
        {
            onegrid[j] = lb + j*es/(num_points);
        }
        grids[i] = onegrid;
        dx[i] = es/(num_points*2);
    }

    // finding start dot over created grid
    std::vector<std::vector<int>> points;
    for(auto it = sample->begin(); it != sample->end(); ++it)
    {
        std::vector<int> startdot(dim);
        for(size_t i = 0; i != startdot.size(); i++)
        {
//            std::vector<T> val(grids[i].size());
//            for(size_t j = 0; j != val.size(); j++)
//            {
//                val[j] = grids[i][j] + dx[i];
//            }
//            auto pos1 = std::lower_bound(val.begin(), val.end(), (*it)[i]);
//            startdot[i] = std::distance(val.begin(), pos1) - 1;
            auto pos1 = std::lower_bound(grids[i].begin(), grids[i].end(), (*it)[i]);
            if(pos1 == grids[i].end())
                startdot[i] = grids[i].size() - 2;
            else if(pos1 == grids[i].begin())
                startdot[i] = 0;
            else
                startdot[i] = std::distance(grids[i].begin(), pos1) - 1;
        }
        points.push_back(startdot);
    }

    std::set<std::vector<int>> visited;

    std::function<T(const std::vector<T> &)> pdffunc = std::bind(&mveqf::kde::KDE<T>::pdf, std::ref(obj), std::placeholders::_1, 1.0);

    auto triess = std::make_shared<mveqf::trie_based::Trie<mveqf::trie_based::NodeCount<int>,int>>();
    triess->set_dimension(gridn.size());

    size_t counter = 0, fe_count = 0;
    std::cout << "FloodFill" << std::endl;
    mveqf::mvff::FloodFill_MultipleGrids_VonNeumann_trie<T>(grids, points, triess, dx, counter, fe_count, pdffunc, threshold, multi);
    std::cout << "FloodFill" << std::endl;

    std::vector<std::vector<T> > sampled;
    std::vector<std::vector<T> > values01;

    mveqf::ImplicitTrieKQuantile<int, T> quant(std::vector<T>(gridn.size(), lb), std::vector<T>(gridn.size(), ub), gridn, kt);
    quant.set_sample_shared(triess);

    timer::Timer time_cpp11;
    time_cpp11.reset();

    std::vector<T> temp1(gridn.size());
    std::vector<T> temp2(temp1.size());


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

    /*for(size_t i = 0; i != temp1.size(); ++i)
    {
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        temp1[i] = 0.0;
        quant.ktransform(temp1,temp2);
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
        quant.ktransform(temp1,temp2);
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
        quant.ktransform(temp1,temp2);
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
        quant.ktransform(temp1,temp2);
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
        quant.ktransform(temp1,temp2);
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
        quant.ktransform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
    }

    for(size_t i = 0; i != sampled.size(); i++)
    {
        for(size_t j = 0; j != sampled[i].size(); j++)
        {
            if(sampled[i][j] < lb || sampled[i][j] > ub)
            {
                std::cout << "out of bounds!" << std::endl;
            }
        }
    }*/


    sample->clear();
    output.clear();

    for(size_t i = 0; i != sampled.size(); i++)
    {
        sample->push_back(sampled[i]);
    }

    mveqf::kde::KDE<T> obj2;
    obj2.set_dimension(dim);
    obj2.set_kernel_type(kt);
    obj2.set_sample_shared(sample);
    for(size_t i = 0; i != sampled.size(); i++)
    {
        auto t = sampled[i];
        T res = obj2.pdf(t);
        sampled[i].push_back(res);
    }
    data_io::write_default2d("maps/Trie/3d/sampled3d.dat", sampled, 5);

}

#endif
