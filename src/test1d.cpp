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
#include "test1d.h"
#include "timer.h"
#include "data_io.h"
#include <vector>
#include <iostream>
#include <random>
#include "kquantile.h"

void test_1d1()
{
    std::vector<size_t> grid_number = {6};

    std::vector<std::vector<int>> sample_implicit = {{0},{1},{2},{3},{4},{5}};

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(-2, 4, grid_number, sample_implicit, 1e+3);
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


void test_1d5()
{
    std::vector<size_t> grid_number = {6};

    std::vector<std::vector<int>> sample_implicit = {{3}};

    /// multivariate quantile function [0,1]^n -> [-3,3]^n
    explicit_quantile(-2, 4, grid_number, sample_implicit, 1e3);
    implicit_quantile_class(-2, 4, grid_number, sample_implicit, 1e3);
    implicit_quantile_class_sorted(-2, 4, grid_number, sample_implicit, 1e3);
}


double threeExp1d(double x)
{
    std::vector<std::vector<double>> centers = { {3, 3}, {-5, 0}, {0, -5} };
    double rez = 0;
    rez += std::exp(-(std::pow(x - centers[0][0], 2.0))*0.75);
    rez += std::exp(-(std::pow(x - centers[1][0], 2.0))*0.5)*0.75;
    rez += std::exp(-(std::pow(x - centers[2][0], 2.0))*0.25)*0.5;
    return rez;
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

    mveqf::ExplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
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


void test_1d_uniform_vs_nonuniform()
{
    std::vector<size_t> gridn = {12};

    std::vector<std::vector<int>> sample;
    
//    sample.push_back(std::vector<int>{1});
//    sample.push_back(std::vector<int>{2});
//    sample.push_back(std::vector<int>{3});
//    sample.push_back(std::vector<int>{7});
//    sample.push_back(std::vector<int>{8});
//    sample.push_back(std::vector<int>{9});
//    sample.push_back(std::vector<int>{10});
    
    
    // nonuniform
    sample.push_back(std::vector<int>{1});
    sample.push_back(std::vector<int>{1});
    sample.push_back(std::vector<int>{1});
    sample.push_back(std::vector<int>{2});
    sample.push_back(std::vector<int>{3});
    sample.push_back(std::vector<int>{3});
    sample.push_back(std::vector<int>{7});
    sample.push_back(std::vector<int>{7});
    sample.push_back(std::vector<int>{7});
    sample.push_back(std::vector<int>{8});
    sample.push_back(std::vector<int>{8});
    sample.push_back(std::vector<int>{8});
    sample.push_back(std::vector<int>{8});
    sample.push_back(std::vector<int>{9});
    sample.push_back(std::vector<int>{9});
    sample.push_back(std::vector<int>{10});
    
    typedef mveqf::trie_based::Trie<mveqf::trie_based::NodeCount<int>,int> trie_type;
    std::shared_ptr<trie_type> trie_sample = std::make_shared<trie_type>();
    trie_sample->insert(std::vector<int>{1},3);
    trie_sample->insert(std::vector<int>{2},1);
    trie_sample->insert(std::vector<int>{3},2);
    trie_sample->insert(std::vector<int>{7},3);
    trie_sample->insert(std::vector<int>{8},4);
    trie_sample->insert(std::vector<int>{9},2);
    trie_sample->insert(std::vector<int>{10},1);
    
    std::cout << sample.size() << std::endl;

    size_t nrolls = 1e+5;

    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    std::vector<std::vector<float> > sampled;
    std::vector<std::vector<float> > values01;

    float lb = -2.0, ub = 4.0;

//    empirical_quantile::ImplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
//    quant.set_sample_and_fill_count(sample);
    
//    quant.set_sample_shared_and_fill_count(shared_sample);
//    quant.set_sample_shared(shared_sample);

//    empirical_quantile::ExplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
//    quant.set_sample(sample);

    mveqf::ImplicitTrieQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
    quant.set_sample_shared(trie_sample);

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






void test_1d_kquantile()
{
    std::vector<size_t> gridn = {12};

    typedef mveqf::trie_based::Trie<mveqf::trie_based::NodeCount<int>,int> trie_type;
    std::shared_ptr<trie_type> trie_sample = std::make_shared<trie_type>();
//    trie_sample->insert(std::vector<int> {2},1);
//    trie_sample->insert(std::vector<int> {3},2);
//    trie_sample->insert(std::vector<int> {4},3);
//    trie_sample->insert(std::vector<int> {5},2);
//    trie_sample->insert(std::vector<int> {6},1);

    trie_sample->insert(std::vector<int> {3},1);
    trie_sample->insert(std::vector<int> {8},1);

//    trie_sample->insert(std::vector<int> {6},1);

    size_t nrolls = 1e+3;

    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<float> ureal01(0.0,1.0);

    std::vector<std::vector<float> > sampled;
    std::vector<std::vector<float> > values01;

    float lb = -2.0, ub = 4.0;

//    empirical_quantile::ImplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
//    quant.set_sample_and_fill_count(sample);

//    quant.set_sample_shared_and_fill_count(shared_sample);
//    quant.set_sample_shared(shared_sample);

//    empirical_quantile::ExplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
//    quant.set_sample(sample);

    mveqf::ImplicitTrieKQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn, 0);
    quant.set_sample_shared(trie_sample);

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

    ///


//    auto kde_sample = std::make_shared<std::vector<std::vector<float>>>();
//
//
//    kde_sample->
//
//
//
//    kde::KDE<float> obj;
//    obj.set_dimension(1);
//    obj.set_kernel_type(0);
//    obj.set_sample_shared(kde_sample);
//
//    float es = ub - lb;

    /*const size_t n = 1e3;
    std::vector<float> pdf(n), cdf(n);
    for(size_t i = 0; i != pdf.size(); i++)
    {
        pdf[i] = obj.pdf(std::vector<float> {lb + i*es/(n - 1)});
        cdf[i] = obj.cdf(std::vector<float> {lb + i*es/(n - 1)});
    }
    data_io::write_default1d("maps/kde/kpdf.dat", pdf, 1, 5);
    data_io::write_default1d("maps/kde/kcdf.dat", cdf, 1, 5);*/
}