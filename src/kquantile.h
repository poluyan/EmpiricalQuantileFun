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
#ifndef KQUANTILE_H
#define KQUANTILE_H

#include "quantile.h"
#include "kde.h"
#include "mvff.h"

namespace mveqf
{

namespace kquantile
{

template <typename T, typename U>
class Qkde : public kde::Kernels<U>
{
protected:
    using kde::Kernels<U>::compute_pdf;
    using kde::Kernels<U>::compute_cdf;

    using kde::Kernels<U>::gaussian_pdf;
    using kde::Kernels<U>::gaussian_cdf;

    using kde::Kernels<U>::epanechnikov_pdf;
    using kde::Kernels<U>::epanechnikov_cdf;

    using kde::Kernels<U>::uniform_pdf;
    using kde::Kernels<U>::uniform_cdf;

    using kde::Kernels<U>::biweight_pdf;
    using kde::Kernels<U>::biweight_cdf;

    using kde::Kernels<U>::triweight_pdf;
    using kde::Kernels<U>::triweight_cdf;

    size_t count;
    U sum, ssum, min, max, bandwidth;

    size_t kernel_type; // 0 - gauss, 1 - epanechnikov, 2 - uniform, 3 - biweight, 4 - triweight

    void calculate_bandwidth(trie_based::NodeCount<T> *layer)
    {
        U x  = sum/count;
        U y = ssum/count;
        U sigma = std::sqrt(y - std::pow(x, 2.0));
        bandwidth = sigma*(std::pow((3.0*count/4.0), (-1.0/5.0)));
        if(layer->children.size() == 1)
        {
            bandwidth = 1.0;
        }
    }
public:
    Qkde():sum(0), ssum(0), min(std::numeric_limits<U>::max()), max(std::numeric_limits<U>::min()) {}
    void set_kernel_type(size_t kt)
    {
        if(kt >= 0 && kt < 5)
            kernel_type = kt;
        else
            throw std::logic_error("kernel type");
    }
    void set_sample(trie_based::NodeCount<T> *layer, size_t ind, const std::vector<std::vector<U>> &grids, const std::vector<U> &dx)
    {
        for(size_t i = 0; i != layer->children.size(); i++)
        {
            U sample = grids[ind][layer->children[i]->index] + dx[ind];
//            for(size_t j = 0; j != layer->children[i]->count; j++)
//            {
//                sum += sample;
//                ssum += std::pow(sample, 2.0);
//            }

            sum += sample*layer->children[i]->count;
            ssum += std::pow(sample, 2.0)*layer->children[i]->count;

            min = sample < min ? sample : min;
            max = sample > max ? sample : max;
        }
        count = layer->count;
        calculate_bandwidth(layer);
    }
    U cdf(U x, trie_based::NodeCount<T> *layer, size_t ind, const std::vector<std::vector<U>> &grids, const std::vector<U> &dx) const
    {
//        if(layer->children.size() < 100)
//        {
        U d = 0;
        for(size_t i = 0; i != layer->children.size(); i++)
        {
            d += compute_cdf(kernel_type, x, grids[ind][layer->children[i]->index] + dx[ind], bandwidth)*layer->children[i]->count;
        }
        return d/layer->count;
//        }
//        else
//        {
//            const auto first = layer->children.begin();
//            const auto last = layer->children.end();
//            U res = 0.0;
//            const auto size = last - first;
//            const auto nthreads = std::thread::hardware_concurrency();
//            const auto size_per_thread = size / nthreads;
//
//            std::vector<std::future<U>> futures;
//            for(unsigned int i = 0; i < nthreads - 1; i++)
//            {
//                futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, &grids, &dx, &layer, this, x, ind]()
//                {
//                    U d = 0;
//                    for(auto it = start; it != start + size_per_thread; ++it)
//                    {
//                        size_t i = std::distance(layer->children.begin(), it);
//                        d += compute_cdf(this->kernel_type, x, grids[ind][layer->children[i]->index] + dx[ind], this->bandwidth)*layer->children[i]->count;
//                    }
//                    return d;
//                }));
//            }
//            futures.emplace_back(
//                std::async([start = first + (nthreads - 1) * size_per_thread, last, &grids, &dx, &layer, this, x, ind]()
//            {
//                U d = 0;
//                for(auto it = start; it != last; ++it)
//                {
//                    size_t i = std::distance(layer->children.begin(), it);
//                    d += compute_cdf(this->kernel_type, x, grids[ind][layer->children[i]->index] + dx[ind], this->bandwidth)*layer->children[i]->count;
//                }
//                return d;
//            }));
//
//            for(auto &&future : futures)
//            {
//                if(future.valid())
//                {
//                    res += future.get();
//                }
//                else
//                {
//                    throw std::runtime_error("Something going wrong.");
//                }
//            }
//            return res/layer->count;
//        }







//        for(size_t i = 0; i != layer->children.size(); i++)
//        {
//            for(size_t j = 0; j != layer->children[i]->count; j++)
//            {
//                d += compute_cdf(kernel_type, x, grids[ind][layer->children[i]->index] + dx[ind], bandwidth);
//            }
//        }

//        for(size_t i = 0; i != layer->children.size(); i++)
//        {
//            U t = compute_cdf(kernel_type, x, grids[ind][layer->children[i]->index] + dx[ind], bandwidth);
//            for(size_t j = 0; j != layer->children[i]->count; j++)
//                d += t;
//        }

//        return d/layer->count;
    }
};

}

template <typename T, typename U>
class ImplicitTrieKQuantile : public mveqf::ImplicitQuantile<T, U>
{
protected:
    using mveqf::ImplicitQuantile<T, U>::grids;
    using mveqf::ImplicitQuantile<T, U>::dx;

    typedef trie_based::Trie<trie_based::NodeCount<T>,T> trie_type;
    std::shared_ptr<trie_type> sample;

    using mveqf::ImplicitQuantile<T, U>::count_less;
    using mveqf::ImplicitQuantile<T, U>::quantile_transform;

    size_t kernel_type;

    U kquantile_transform(trie_based::NodeCount<T> *layer, size_t ind, U val01) const;
public:
    ImplicitTrieKQuantile();
    ImplicitTrieKQuantile(std::vector<U> in_lb, std::vector<U> in_ub, std::vector<size_t> in_gridn, size_t kt);
    void set_kernel_type(size_t kt);
    void set_sample_shared(std::shared_ptr<trie_type> in_sample);
    void transform(const std::vector<U>& in01, std::vector<U>& out) const;
    std::vector<U> transform(const std::vector<U>& in01) const;
};
template <typename T, typename U>
ImplicitTrieKQuantile<T, U>::ImplicitTrieKQuantile() : kernel_type(0)
{
}
template <typename T, typename U>
ImplicitTrieKQuantile<T, U>::ImplicitTrieKQuantile(std::vector<U> in_lb,
        std::vector<U> in_ub,
        std::vector<size_t> in_gridn,
        size_t kt) : mveqf::ImplicitQuantile<T, U>(in_lb, in_ub, in_gridn), kernel_type(kt)
{
}
template <typename T, typename U>
void ImplicitTrieKQuantile<T, U>::set_sample_shared(std::shared_ptr<trie_type> in_sample)
{
    sample = std::move(in_sample);
}
template <typename T, typename U>
void ImplicitTrieKQuantile<T, U>::set_kernel_type(size_t kt)
{
    kernel_type = kt;
}

template <typename T, typename U>
void ImplicitTrieKQuantile<T, U>::transform(const std::vector<U>& in01, std::vector<U>& out) const
{
    auto p = sample->root.get();
    for(size_t i = 0, k; i != in01.size(); i++)
    {
        out[i] = kquantile_transform(p, i, in01[i]);
        k = quantile_transform(p, i, in01[i]).first;
        p = p->children[k].get();
    }
}

template <typename T, typename U>
std::vector<U> ImplicitTrieKQuantile<T, U>::transform(const std::vector<U>& in01) const
{
    std::vector<U> out(grids.size());
    transform(in01, out);
    return out;
}

template <typename T, typename U>
U ImplicitTrieKQuantile<T, U>::kquantile_transform(trie_based::NodeCount<T> *layer, size_t ind, U val01) const
{

//    timer::Timer time;


    kquantile::Qkde<T, U> obj;
    obj.set_kernel_type(kernel_type);
    obj.set_sample(layer, ind, grids, dx);

//    std::cout << std::fixed << '\t' << time.elapsed_seconds() << std::endl;
//    time.reset();

    const U lower_bound = grids[ind].front();
    const U upper_bound = grids[ind].back();

    const U es = upper_bound - lower_bound;
    const U min = obj.cdf(lower_bound, layer, ind, grids, dx);
    const U max = obj.cdf(upper_bound, layer, ind, grids, dx);

    const size_t n = 100000;
    size_t c = n - 1, step;
    U f1, f2, p1, p2;
    size_t f = 0, i = 0;
    while(c > 0)
    {
        i = f;
        step = c / 2;
        i = i + step;

        p1 = lower_bound + es * i / (n - 1);
        f1 = (obj.cdf(p1, layer, ind, grids, dx) - min)/(max - min);
        if(f1 < val01)
        {
            ++i;
            f = i;
            c -= step + 1;
        }
        else
            c = step;
    }

    p2 = lower_bound + es * (i + 1) / (n - 1);
    f2 = (obj.cdf(p2, layer, ind, grids, dx) - min)/(max - min);

    U check = p1 + (val01 - f1)*(p2 - p1)/(f2 - f1);
    if(!std::isfinite(check))
    {
        size_t max_ind = 0;
        for(size_t j = 0; j != n - 1; j++)
        {
            p2 = lower_bound + es * (j + 1) / (n - 1);
            f2 = (obj.cdf(p2, layer, ind, grids, dx) - min)/(max - min);
            if(f2 < 1.0)
                max_ind = j;
            if(f2 > val01)
            {
                p1 = lower_bound + es * j / (n - 1);
                f1 = (obj.cdf(p1, layer, ind, grids, dx) - min)/(max - min);
                return p1 + (val01 - f1)*(p2 - p1)/(f2 - f1);
            }
        }
        std::cout << "fff" << std::endl;// std::fixed << '\t' << time.elapsed_seconds() << std::endl;
//
//        for(size_t j = 0; j != layer->children.size(); j++)
//        {
//            std::cout << layer->children[i]->count << std::endl;
//        }
//

//        std::cin.get();
        return lower_bound + es * max_ind / (n - 1);
    }

//    std::cout << std::fixed << '\t' << time.elapsed_seconds() << std::endl;

    return p1 + (val01 - f1)*(p2 - p1)/(f2 - f1);
}


void test_1d_kquantile()
{
    std::vector<size_t> gridn = {12};

    typedef trie_based::Trie<trie_based::NodeCount<int>,int> trie_type;
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

template<class InputIt, class T>
void parallel_step(InputIt first,
                   InputIt last,
                   std::vector<std::vector<T>> &pdf,
                   const T lb,
                   const T ub,
                   std::function<T(const std::vector<T> &)> f)
{
    T es = ub - lb;
    for(auto it = first; it != last; ++it)
    {
        size_t i = std::distance(pdf.begin(), it);
        T a = lb + i*es/(pdf.size() - 1);
        for(size_t j = 0; j != it->size(); j++)
        {
            std::vector<T> temp = {a, lb + j*es/(pdf.size() - 1)};
            pdf[i][j] = f(temp);
        }
    }
}

template<class InputIt, class T>
void parallel_pdf(InputIt first,
                  InputIt last,
                  std::vector<std::vector<T>> &pdf,
                  const T lb,
                  const T ub,
                  std::function<T(const std::vector<T> &)> f)
{
    const auto size = last - first;
    const auto nthreads = std::thread::hardware_concurrency();
    const auto size_per_thread = size / nthreads;

    std::vector<std::future<void>> futures;
    for(unsigned int i = 0; i < nthreads - 1; i++)
    {
        futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, &pdf, lb, ub, f]()
        {
            T es = ub - lb;
            for(auto it = start; it != start + size_per_thread; ++it)
            {
                size_t i = std::distance(pdf.begin(), it);
                T a = lb + i*es/(pdf.size() - 1);
                for(size_t j = 0; j != it->size(); j++)
                {
                    std::vector<T> temp = {a, lb + j*es/(pdf.size() - 1)};
                    pdf[i][j] = f(temp);
                }
            }
        }));
    }
    futures.emplace_back(
        std::async([start = first + (nthreads - 1) * size_per_thread, last, &pdf, lb, ub, f]()
    {
        T es = ub - lb;
        for(auto it = start; it != last; ++it)
        {
            size_t i = std::distance(pdf.begin(), it);
            T a = lb + i*es/(pdf.size() - 1);
            for(size_t j = 0; j != it->size(); j++)
            {
                std::vector<T> temp = {a, lb + j*es/(pdf.size() - 1)};
                pdf[i][j] = f(temp);
            }
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
}

template <typename T>
void print_pdf(std::string fname, std::function<T(const std::vector<T> &)> f, const T lb, const T ub, const size_t n)
{
    std::vector<std::vector<T> > pdf(n, std::vector<T>(n));
    //parallel_step(pdf.begin(), pdf.end(), pdf, lb, ub, f);
    parallel_pdf(pdf.begin(), pdf.end(), pdf, lb, ub, f);
    data_io::write_default2d(fname, pdf, 5);
}

template <typename T>
void test_2d_kquantile()
{
    const size_t kt = 0;
    const size_t dim = 2;
    const size_t N = 1000;
    const size_t nrolls = 10000;
    const T threshold = 1e-2;
    const size_t multi = 1000000;
    std::vector<size_t> gridn(dim, N);

//    auto sample = std::make_shared<trie_based::TrieBased<trie_based::NodeCount<int>,int>>();
//    sample->set_dimension(gridn.size());
//    sample->insert(std::vector<int>{5,5});
//    sample->insert(std::vector<int>{4,5});
//    sample->insert(std::vector<int>{5,4});
//    sample->insert(std::vector<int>{1,3});



    auto sample = std::make_shared<trie_based::Trie<trie_based::NodeCount<int>,int>>();
    sample->set_dimension(gridn.size());
    // 10x10
    /*sample->insert(std::vector<int>{3,2}, 1);
    sample->insert(std::vector<int>{2,3}, 1);
    sample->insert(std::vector<int>{3,4}, 1);
    sample->insert(std::vector<int>{4,3}, 1);
    sample->insert(std::vector<int>{3,3}, 2);*/

    ///


    auto kdesample = std::make_shared<std::vector<std::vector<T>>>();

    ///first
//    sample->push_back(std::vector<T>{0.0, 0.0});

    /// second
    // 1125, 1126
    /*std::vector<std::vector<size_t> > dots =
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
        std::vector<T> t = {-1.75 + 3.5*(i.front())/T(1125), -1.75 + 3.5*(T(1126) - i.back())/T(1126)};
        std::cout << std::scientific << t.front() << '\t' << t.back() << std::endl;
        kdesample->push_back(t);
    }*/
    
    std::ifstream gridIn;
    gridIn.open("maps/Rpics/sampled_implicit.dat");
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

        for(size_t i = 0; i != t.size(); i++)
        {
            if(t[i] < lbv[i])
                lbv[i] = t[i];
            if(t[i] > ubv[i])
                ubv[i] = t[i];
        }

        kdesample->push_back(t);
        if(kdesample->size() > 10000)
            break;
    }
    gridIn.close();
    

    const T lb = -2, ub = 2;

    kde::KDE<T> obj;
    obj.set_dimension(dim);
    obj.set_kernel_type(kt);
    obj.set_sample_shared(kdesample);

    std::function<T(const std::vector<T> &)> fff = std::bind(&kde::KDE<T>::pdf, std::ref(obj), std::placeholders::_1);

    print_pdf<T>("maps/Trie/kpdf2d_init.dat", fff, lb, ub, 1000);

    ///

    std::vector<std::vector<T>> grids(gridn.size());
    std::vector<T> dx(gridn.size());

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
    for(auto it = kdesample->begin(); it != kdesample->end(); ++it)
    {
        std::vector<int> startdot(dim);
        for(size_t i = 0; i != startdot.size(); i++)
        {
            std::vector<T> val(grids[i].size());
            for(size_t j = 0; j != val.size(); j++)
            {
                val[j] = grids[i][j] + dx[i];
            }
            auto pos1 = std::lower_bound(val.begin(), val.end(), (*it)[i]);
            startdot[i] = std::distance(val.begin(), pos1) - 1;
        }
//        std::cout << startdot.front() << '\t' << startdot.back() << std::endl;
//        std::cin.get();
        points.push_back(startdot);
    }

    std::set<std::vector<int>> visited;

    //std::vector<std::vector<T> > samples;

    std::function<T(const std::vector<T> &)> pdffunc = std::bind(&kde::KDE<T>::pdf, std::ref(obj), std::placeholders::_1);
//    std::cout << pdffunc(std::vector<float>{-1.0,1.0}) << std::endl;
//    std::cout << obj.pdf(std::vector<float>{-1.0,1.0}) << std::endl;
//
//    std::cin.get();

    size_t counter = 0, fe_count = 0;
    std::cout << "FloodFill start" << std::endl;
    mveqf::mvff::FloodFill_MultipleGrids_VonNeumann_trie<T>(grids, points, sample, dx, counter, fe_count, pdffunc, threshold, multi);
    std::cout << "FloodFill end" << std::endl;

    //std::cout << samples.size() << std::endl;
    //std::cin.get();
    ///

    std::mt19937_64 generator;
    generator.seed(1);
    std::uniform_real_distribution<T> ureal01(0.0,1.0);

    std::vector<std::vector<T> > sampled;
    std::vector<std::vector<T> > values01;

    mveqf::ImplicitTrieKQuantile<int, T> quant(std::vector<T>(gridn.size(), lb), std::vector<T>(gridn.size(), ub), gridn, kt);
    quant.set_sample_shared(sample);

    timer::Timer time_cpp11;
    time_cpp11.reset();

    std::vector<T> temp1(gridn.size());
    std::vector<T> temp2(temp1.size());

    std::cout << "Transform start" << std::endl;
    for(size_t i = 0; i != nrolls; ++i)
    {
//        timer::Timer t;
        for(size_t j = 0; j != temp1.size(); j++)
        {
            temp1[j] = ureal01(generator);
        }
        quant.transform(temp1,temp2);
        values01.push_back(temp1);
        sampled.push_back(temp2);
//        std::cout << std::fixed << t.elapsed_seconds() << std::endl;
    }
    std::cout << "Transform end" << std::endl;
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

    data_io::write_default2d("maps/sampled_implicit.dat", sampled, 5);

    std::shared_ptr<std::vector<std::vector<T> >> ss = std::make_shared<std::vector<std::vector<T>>>();
    for(const auto & i : sampled)
        ss->push_back(i);

    kde::KDE<T> obj1;
    obj1.set_dimension(dim);
    obj1.set_kernel_type(kt);
    obj1.set_sample_shared(ss);
    std::function<T(const std::vector<T> &)> ff = std::bind(&kde::KDE<T>::pdf, std::ref(obj1), std::placeholders::_1);

    print_pdf<T>("maps/Trie/kpdf2d.dat", ff, lb, ub, 1000);
}

template <typename T>
void test_3d_kquantile()
{
    const size_t kt = 3;
    const size_t dim = 3;
    const size_t N = 200;
    const size_t nrolls = 10000;
    const T threshold = 1e-8;
    const size_t multi = 10000;

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


    kde::KDE<T> obj1;
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



    kde::KDE<T> obj;
    obj.set_dimension(dim);
    obj.set_kernel_type(kt);
    obj.set_sample_shared(sample);


    ///
    std::vector<std::vector<T>> grids(gridn.size());
    std::vector<T> dx(gridn.size());

    const T lb = -31, ub = 31;
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
            std::vector<T> val(grids[i].size());
            for(size_t j = 0; j != val.size(); j++)
            {
                val[j] = grids[i][j] + dx[i];
            }
            auto pos1 = std::lower_bound(val.begin(), val.end(), (*it)[i]);
            startdot[i] = std::distance(val.begin(), pos1) - 1;
        }
        points.push_back(startdot);
    }

    std::set<std::vector<int>> visited;

    std::function<T(const std::vector<T> &)> pdffunc = std::bind(&kde::KDE<T>::pdf, std::ref(obj), std::placeholders::_1);

    auto triess = std::make_shared<trie_based::Trie<trie_based::NodeCount<int>,int>>();
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

    kde::KDE<T> obj2;
    obj2.set_dimension(dim);
    obj2.set_kernel_type(kt);
    obj2.set_sample_shared(sample);
    for(size_t i = 0; i != sampled.size(); i++)
    {
        auto t = sampled[i];
        sampled[i].push_back(obj2.pdf(t));
    }
    data_io::write_default2d("maps/Trie/3d/sampled3d.dat", sampled, 5);

}

}

#endif
