/**************************************************************************

   Copyright © 2020 Sergey Poluyan <svpoluyan@gmail.com>

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

#include "kde.h"
#include "data_io.h"
#include <random>

namespace mveqf
{

namespace kde
{

void test1d()
{
    auto sample = std::make_shared<std::vector<std::vector<double>>>();

    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<double> ureal01(0.0,1.0);
    std::normal_distribution<double> norm(0.0,1.0);

    /// first
    for(size_t i = 0; i != 1e3; i++)
        sample->push_back(std::vector<double> {norm(generator)});

    /// second
//    sample->push_back(std::vector<double>{1.0});

    /// third
//    sample->push_back(std::vector<double>{-2.0});
//    sample->push_back(std::vector<double>{2.0});


    mveqf::kde::KDE<double> obj;
    obj.set_dimension(1);
    obj.set_kernel_type(0);
    obj.set_sample_shared(sample);

    const double lb = -5;
    const double ub = 5;
    const double es = ub - lb;

    const size_t n = 1e3;
    std::vector<double> pdf(n), cdf(n);
    for(size_t i = 0; i != pdf.size(); i++)
    {
        pdf[i] = obj.pdf(std::vector<double> {lb + i*es/(n - 1)});
        cdf[i] = obj.cdf(std::vector<double> {lb + i*es/(n - 1)});
    }
    std::cout << obj.pdf(std::vector<double> {1.0}) << std::endl;
    data_io::write_default1d("maps/kde/kpdf.dat", pdf, 1, 5);
    data_io::write_default1d("maps/kde/kcdf.dat", cdf, 1, 5);
}

void test1d_1()
{
    auto sample = std::make_shared<std::vector<std::vector<double>>>();

    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<double> ureal01(0.0,1.0);
    std::normal_distribution<double> norm(0.0,1.0);

    /// first
//    for(size_t i = 0; i != 1e3; i++)
//        sample->push_back(std::vector<double> {norm(generator)});

    /// second
//    sample->push_back(std::vector<double>{1.0});

    /// third
    sample->push_back(std::vector<double>{-2.0});
    sample->push_back(std::vector<double>{2.0});
    auto count = std::make_shared<std::vector<size_t>>(std::vector<size_t>{2, 1});


    mveqf::kde::KDE<double> obj;
    obj.set_dimension(1);
    obj.set_kernel_type(4);
//    obj.set_sample_shared(sample);
    obj.set_sample_shared(sample, count);

    const double lb = -5;
    const double ub = 5;
    const double es = ub - lb;

    const size_t n = 1e3;
    std::vector<double> pdf(n), cdf(n);
    for(size_t i = 0; i != pdf.size(); i++)
    {
        pdf[i] = obj.pdf(std::vector<double> {lb + i*es/(n - 1)});
        cdf[i] = obj.cdf(std::vector<double> {lb + i*es/(n - 1)});
    }
    std::cout << obj.pdf(std::vector<double> {1.0}) << std::endl;
    data_io::write_default1d("maps/kde/kpdf.dat", pdf, 1, 5);
    data_io::write_default1d("maps/kde/kcdf.dat", cdf, 1, 5);
}

void test2d()
{
    auto sample = std::make_shared<std::vector<std::vector<double>>>();

    ///first
//    sample->push_back(std::vector<double>{0.0, 0.0});

    /// second
    // 1125, 1126
    std::vector<std::vector<double> > dots =
    {
        {484,197},
        {731,256},
        {582,288},
        {505,291},
        {505,378},
        {699,387},
        {799,407},
        {793,439},
        {591,444},
        {736,452},
        {751,460},
        {308,463},
        {170,472},
        {243,483},
        {406,496},
        {388,498},
        {442,498},
        {691,511},
        {485,516},
        {598,518},
        {522,529},
        {32,532},
        {472,534},
        {578,535},
        {596,549},
        {566,558},
        {479,579},
        {209,587},
        {858,594},
        {601,621},
        {787,623},
        {749,635},
        {503,637},
        {878,644},
        {428,651},
        {585,652},
        {494,654},
        {687,654},
        {803,667},
        {650,684},
        {595,700},
        {315,702},
        {549,702},
        {561,723},
        {596,731},
        {625,735},
        {664,742},
        {731,744},
        {732,788},
        {736,863}
    };
    for(const auto & i : dots)
    {
        std::vector<double> t = {-1.75 + 3.5*(i.front())/double(1125), -1.75 + 3.5*(double(1126) - i.back())/double(1126)};
        sample->push_back(t);
    }

    KDE<double> obj;
    obj.set_dimension(2);
    obj.set_kernel_type(0);
    obj.set_sample_shared(sample);

    double lb = -2.0;
    double ub = 2.0;
    double es = ub - lb;

    size_t n = 1000;

    std::vector<std::vector<double> > pdf(n, std::vector<double>(n));
    auto cdf = pdf;
    for(size_t i = 0; i != pdf.size(); i++)
    {
        double a = lb + i*es/(n - 1);
        for(size_t j = 0; j != pdf[i].size(); j++)
        {
            std::vector<double> temp = {a, lb + j*es/(n - 1)};
            cdf[i][j] = obj.cdf(temp);
            pdf[i][j] = obj.pdf(temp);
        }
    }
    data_io::write_default2d("maps/kde/kpdf2d.dat", pdf, 5);
    data_io::write_default2d("maps/kde/kcdf2d.dat", cdf, 5);
}

void test2d_2()
{
    auto sample = std::make_shared<std::vector<std::vector<double>>>();

    ///first
//    sample->push_back(std::vector<double>{-0.5, -0.25});
//    sample->push_back(std::vector<double>{-0.5, -0.25});
//    sample->push_back(std::vector<double>{-0.5, -0.25});
//    sample->push_back(std::vector<double>{0.0, 0.75});
//    sample->push_back(std::vector<double>{0.0, 0.75});
//    sample->push_back(std::vector<double>{0.25, 0.0});
    
    sample->push_back(std::vector<double>{-0.5, -0.25});
    sample->push_back(std::vector<double>{0.0, 0.75});
    sample->push_back(std::vector<double>{0.25, 0.0});
    auto count = std::make_shared<std::vector<size_t>>(std::vector<size_t>{3, 2, 1});

    KDE<double> obj;
    obj.set_dimension(2);
    obj.set_kernel_type(0);
//    obj.set_sample_shared(sample);
    obj.set_sample_shared(sample, count);
    double lb = -2.0;
    double ub = 2.0;
    double es = ub - lb;

    size_t n = 1000;

    std::vector<std::vector<double> > pdf(n, std::vector<double>(n));
    auto cdf = pdf;
    for(size_t i = 0; i != pdf.size(); i++)
    {
        double a = lb + i*es/(n - 1);
        for(size_t j = 0; j != pdf[i].size(); j++)
        {
            std::vector<double> temp = {a, lb + j*es/(n - 1)};
            cdf[i][j] = obj.cdf(temp);
            pdf[i][j] = obj.pdf(temp);
        }
    }
    data_io::write_default2d("maps/kde/kpdf2d.dat", pdf, 5);
    data_io::write_default2d("maps/kde/kcdf2d.dat", cdf, 5);
}

}

}
