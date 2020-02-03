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
#include <test/test.h>
#include <utility/data_io.h>
#include <utility/timer.h>
#include <test/testNd.h>
#include <vector>
#include <random>
#include <chrono>
#include <thread>

namespace testNd
{


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

    typedef mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<int>,int> sample_type;
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

    mveqf::ExplicitQuantile<int, float> quant_expl(lb, ub, gridN);
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

    mveqf::ImplicitQuantile<int, float> quant_impl(lb, ub, gridN);
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

    mveqf::ImplicitQuantileSorted<int, float> quant_impls(lb, ub, gridN);
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

    typedef mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<int>,int> sample_type;
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

    mveqf::ImplicitQuantile<int, float> quant_impl(lb, ub, gridN);
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

    mveqf::ImplicitQuantileSorted<int, float> quant_impls(lb, ub, gridN);
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



std::vector<double> worst_space(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nrolls, size_t seed_append)
{
    std::mt19937_64 generator;
    generator.seed(1 + gridN.size() + gridN.front() + seed_append);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    typedef mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<std::uint8_t>,std::uint8_t> sample_type;
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

    mveqf::ExplicitQuantile<std::uint8_t, float> quant_expl(lb, ub, gridN);
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

    mveqf::ImplicitQuantile<std::uint8_t, float> quant_impl(lb, ub, gridN);
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

    mveqf::ImplicitQuantileSorted<std::uint8_t, float> quant_impls(lb, ub, gridN);
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
        typedef mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<std::uint8_t>,std::uint8_t> sample_type;
        std::shared_ptr<sample_type> sample = std::make_shared<sample_type>();

        iterate_trie(variable_values, sample);

        mveqf::ImplicitQuantile<std::uint8_t, float> quant_impl(lb, ub, gridN);
        quant_impl.set_sample_shared_and_fill_count(sample);

        std::cout << "hit it" << std::endl;
        std::chrono::seconds dura(1000);
        std::this_thread::sleep_for(dura);
    }
    else
    {
        std::vector<std::vector<std::uint8_t>> sample_int = iterate(variable_values);
        mveqf::ExplicitQuantile<std::uint8_t, float> quant_expl(lb, ub, gridN);
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
        typedef mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<int>,int> sample_type;
        std::shared_ptr<sample_type> sample = std::make_shared<sample_type>();

        iterate_trie(variable_values, sample);

        mveqf::ImplicitQuantile<int, double> quant_impl(lb, ub, gridN);
        quant_impl.set_sample_shared_and_fill_count(sample);

        std::cout << "hit it" << std::endl;
        std::chrono::seconds dura(1000);
        std::this_thread::sleep_for(dura);
    }
    else
    {
        std::vector<std::vector<int>> sample_int = iterate(variable_values);
        mveqf::ExplicitQuantile<int, double> quant_expl(lb, ub, gridN);
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





std::vector<double> sample_size_procedure(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nsamples, size_t Nrolls, size_t seed_append)
{
    std::mt19937_64 generator;
    generator.seed(1 + gridN.size() + gridN.front() + Nsamples + seed_append);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    typedef mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<std::uint8_t>,std::uint8_t> sample_type;
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

    mveqf::ExplicitQuantile<std::uint8_t, float> quant_expl(lb, ub, gridN);
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

    mveqf::ImplicitQuantile<std::uint8_t, float> quant_impl(lb, ub, gridN);
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

    mveqf::ImplicitQuantileSorted<std::uint8_t, float> quant_impls(lb, ub, gridN);
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




}