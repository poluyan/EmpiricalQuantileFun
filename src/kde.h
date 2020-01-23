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

#ifndef KDE_H
#define KDE_H

#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <future>
//#include "fmath.hpp"

namespace mveqf
{

namespace kde
{

template <typename T>
class Kernels
{
protected:
    const T pi = std::acos(-1.0);
public:
    Kernels() {}
    inline T gaussian_cdf(T x, T mu, T sigma) const
    {
        // 0.5*(1.0 + std::erf((x - mu)/(sigma*std::sqrt(2.0))))
        // Abramowitz Stegun normal CDF
        T t = 1.0/(1.0 + 0.2316419*std::abs(x - mu)/sigma);
        T y = t*(0.319381530 + t*(-0.356563782 + t*(1.781477937 + t*(-1.821255978 + t*1.330274429))));
        if(x >= mu)
        {
            return 1.0 - gaussian_pdf(x, mu, sigma)*y*sigma;
        }
        else
        {
            return gaussian_pdf(x, mu, sigma)*y*sigma;
        }
    }
    inline T epanechnikov_cdf(T x, T mu, T sigma) const
    {
        T z = (x - mu)/sigma;
        if(z < -1.0)
            return 0.0;
        if(z >  1.0)
            return 1.0;
        return 0.25*(2.0 + 3.0*z - std::pow(z, 3.0));
    }
    inline T biweight_cdf(T x, T mu, T sigma) const
    {
        T z = (x - mu)/sigma;
        if(z < -1.0)
            return 0.0;
        if(z >  1.0)
            return 1.0;
        return 0.9375*(z -  2.0*std::pow(z, 3.0)/3.0 + 0.2*std::pow(z, 5.0)) + 0.5;
    }
    inline T triweight_cdf(T x, T mu, T sigma) const
    {
        T z = (x - mu)/sigma;
        if(z < -1.0)
            return 0.0;
        if(z >  1.0)
            return 1.0;
        return 1.09375*(z - std::pow(z, 3.0) + 0.6*std::pow(z, 5.0) - std::pow(z, 7.0)/7.0) + 0.5;
    }
    inline T uniform_cdf(T x, T mu, T sigma) const
    {
        if(x < mu - 0.5*sigma)
            return 0.0;
        if(x > mu + 0.5*sigma)
            return 1.0;
        return (x-mu)/sigma + 0.5;
    }
    inline T gaussian_pdf(T x, T mu, T sigma) const
    {
        return std::exp(-0.5*std::pow((x - mu)/sigma, 2.0))/(sigma*std::sqrt(2.0*pi));

//        T y = (x - mu)/sigma;
//        return fmath::expd(-0.5*y*y)/(sigma*std::sqrt(2.0*pi));
    }
    inline T epanechnikov_pdf(T x, T mu, T sigma) const
    {
        T z = (x - mu)/sigma;
        return std::abs(z) > 1.0 ? 0.0 : 0.75*(1.0-std::pow(z, 2.0))/sigma;
    }
    inline T biweight_pdf(T x, T mu, T sigma) const
    {
        T z = (x - mu)/sigma;
        return std::abs(z) > 1.0 ? 0.0 : 0.9375*std::pow(1.0 - std::pow(z, 2.0), 2.0);
    }
    inline T triweight_pdf(T x, T mu, T sigma) const
    {
        T z = (x - mu)/sigma;
        return std::abs(z) > 1.0 ? 0.0 : 1.09375*std::pow(1.0 - std::pow(z, 2.0), 3.0);
    }
    inline T uniform_pdf(T x, T mu, T sigma) const
    {
        return (x < mu - 0.5*sigma || x > mu + 0.5*sigma) ? 0.0 : 1.0/sigma;
    }
    inline T compute_pdf(size_t kt, T x, T mu, T sigma) const
    {
        switch(kt)
        {
            case 1:
                return epanechnikov_pdf(x, mu, sigma);
            case 2:
                return uniform_pdf(x, mu, sigma);
            case 3:
                return biweight_pdf(x, mu, sigma);
            case 4:
                return triweight_pdf(x, mu, sigma);
            default:
                return gaussian_pdf(x, mu, sigma);
        }
    }
    inline T compute_cdf(size_t kt, T x, T mu, T sigma) const
    {
        switch(kt)
        {
            case 1:
                return epanechnikov_cdf(x, mu, sigma);
            case 2:
                return uniform_cdf(x, mu, sigma);
            case 3:
                return biweight_cdf(x, mu, sigma);
            case 4:
                return triweight_cdf(x, mu, sigma);
            default:
                return gaussian_cdf(x, mu, sigma);
        }
    }
};

template <typename T>
class KDE : public kde::Kernels<T>
{
protected:
    using kde::Kernels<T>::compute_pdf;
    using kde::Kernels<T>::compute_cdf;

    using kde::Kernels<T>::gaussian_pdf;
    using kde::Kernels<T>::gaussian_cdf;

    using kde::Kernels<T>::epanechnikov_pdf;
    using kde::Kernels<T>::epanechnikov_cdf;

    using kde::Kernels<T>::uniform_pdf;
    using kde::Kernels<T>::uniform_cdf;

    using kde::Kernels<T>::biweight_pdf;
    using kde::Kernels<T>::biweight_cdf;

    using kde::Kernels<T>::triweight_pdf;
    using kde::Kernels<T>::triweight_cdf;

    typedef std::vector<std::vector<T>> sample_type;
public:
    KDE() {}
    //KDE(const KDE&) = delete;
    KDE(const KDE& f)
    {
        this->dimension = f.dimension;
        this->kernel_type = f.dimension;
        set_sample_shared(f.sample);
    }
    KDE& operator=(const KDE&) = delete;

    void set_kernel_type(size_t kt)
    {
        if(kt >= 0 && kt < 5)
            kernel_type = kt;
        else
            throw std::logic_error("kernel type");
    }
    void set_dimension(size_t dim)
    {
        dimension = dim;
        sum = std::vector<T>(dimension, 0);
        ssum = std::vector<T>(dimension, 0);
        min = std::vector<T>(dimension, std::numeric_limits<T>::max());
        max = std::vector<T>(dimension, std::numeric_limits<T>::min());
        bandwidth = std::vector<T>(dimension);
    }
    void set_sample_shared(std::shared_ptr<sample_type> in_sample)
    {
        sample = std::move(in_sample);

        for(auto it = sample->begin(); it != sample->end(); ++it)
        {
            if((*it).size() != dimension)
                throw std::logic_error("in_sample->front().size() != dimension");

            for(size_t j = 0; j != dimension; j++)
            {
                sum[j] += (*it)[j];
                ssum[j] += std::pow((*it)[j], 2.0);
                min[j] = (*it)[j] < min[j] ? (*it)[j] : min[j];
                max[j] = (*it)[j] > max[j] ? (*it)[j] : max[j];
            }
        }

        count = sample->size();
        calculate_bandwidth();
    }
    T pdf(const std::vector<T> &x) const
    {
        if(sample->size() < 1000000)
        {
            T res = 0.0;
            for(auto it = sample->begin(); it != sample->end(); ++it)
            {
                T t = 1.0;
                for(size_t i = 0; i != dimension; i++)
                {
                    t *= compute_pdf(kernel_type, x[i], (*it)[i], bandwidth[i]);
                }
                res += t;
            }
            return res/count;
        }
        else
        {
            return parallel_pdf(sample->begin(), sample->end(), x);
        }
    }
    T cdf(const std::vector<T> &x) const
    {
        if(sample->size() < 1000000)
        {
            T res = 0.0;
            for(auto it = sample->begin(); it != sample->end(); ++it)
            {
                T t = 1.0;
                for(size_t j = 0; j != (*it).size(); j++)
                {
                    t *= compute_cdf(kernel_type, x[j], (*it)[j], bandwidth[j]);
                }
                res += t;
            }
            return res/count;
        }
        else
        {
            return parallel_cdf(sample->begin(), sample->end(), x);
        }
    }
protected:
    template<typename InputIt>
    T pdf_it(InputIt first, InputIt last, const std::vector<T> &x) const
    {
        T res = 0.0;
        for(auto it = first; it != last; ++it)
        {
            T t = 1.0;
            for(size_t i = 0; i != dimension; i++)
            {
                t *= compute_pdf(kernel_type, x[i], (*it)[i], bandwidth[i]);
            }
            res += t;
        }
        return res;
    }
    template<typename InputIt>
    T cdf(InputIt first, InputIt last, const std::vector<T> &x) const
    {
        T res = 0.0;
        for(auto it = first; it != last; ++it)
        {
            T t = 1.0;
            for(size_t i = 0; i != dimension; i++)
            {
                t *= compute_cdf(kernel_type, x[i], (*it)[i], bandwidth[i]);
            }
            res += t;
        }
        return res;
    }
    template<class InputIt>
    T parallel_pdf(InputIt first, InputIt last, const std::vector<T> &x) const
    {
        T res = 0.0;
        const auto size = last - first;
        const auto nthreads = std::thread::hardware_concurrency();
        const auto size_per_thread = size / nthreads;

        std::vector<std::future<T>> futures;
        for(unsigned int i = 0; i < nthreads - 1; i++)
        {
            futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, x, this]()
            {
                return this->pdf_it(start, start + size_per_thread, x);
            }));
        }
        futures.emplace_back(
            std::async([start = first + (nthreads - 1) * size_per_thread, last, x, this]()
        {
            return this->pdf_it(start, last, x);
        }));

        for(auto &&future : futures)
        {
            if(future.valid())
            {
                res += future.get();
            }
            else
            {
                throw std::runtime_error("Something going wrong.");
            }
        }
        return res/count;
    }

    template<class InputIt>
    T parallel_cdf(InputIt first, InputIt last, const std::vector<T> &x) const
    {
        T res = 0.0;
        const auto size = last - first;
        const auto nthreads = std::thread::hardware_concurrency();
        const auto size_per_thread = size / nthreads;

        std::vector<std::future<T>> futures;
        for(unsigned int i = 0; i < nthreads - 1; i++)
        {
            futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, x, this]()
            {
                return this->cdf(start, start + size_per_thread, x);
            }));
        }
        futures.emplace_back(
            std::async([start = first + (nthreads - 1) * size_per_thread, last, x, this]()
        {
            return this->cdf(start, last, x);
        }));

        for(auto &&future : futures)
        {
            if(future.valid())
            {
                res += future.get();
            }
            else
            {
                throw std::runtime_error("Something going wrong.");
            }
        }
        return res/count;
    }

    size_t dimension;
    size_t kernel_type; // 0 - gaussian, 1 - epanechnikov, 2 - uniform, 3 - biweight, 4 - triweight
    size_t count;
    std::vector<T> sum, ssum, min, max, bandwidth;
    std::shared_ptr<sample_type> sample;

    void calculate_bandwidth()
    {
        for(size_t i = 0; i != dimension; i++)
        {
            T x = sum[i]/count;
            T y = ssum[i]/count;
            T sigma = std::sqrt(y - std::pow(x, 2.0));
            bandwidth[i] = sigma*(std::pow((3.0*count/4.0), (-1.0/5.0)));
        }
        if(sample->size() == 1)
        {
            for(size_t i = 0; i != dimension; i++)
            {
                bandwidth[i] = 1.0;
            }
        }
    }
};

void test1d();
void test2d();
}

}

#endif
