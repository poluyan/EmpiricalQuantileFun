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

#include <quantile.h>
#include <kde.h>
#include <mvff.h>
//#include <utility/timer.h>
//#include <utility/data_io.h>
#include <random>
#include <ctime>

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

			size_t count;
			U sum, ssum, min, max, bandwidth;

			std::vector<U> sample;

			size_t kernel_type; // 0 - gauss, 1 - epanechnikov, 2 - uniform, 3 - biweight, 4 - triweight

			void calculate_bandwidth(trie_based::NodeCount<T> *layer)
			{
				U x = sum/count;
				U y = ssum/count;
				U sigma = std::sqrt(y - std::pow(x, 2.0));
				bandwidth = sigma*(std::pow((3.0*count/4.0), (-1.0/5.0)));
//        std::cout << bandwidth << std::endl;
				if(layer->children.size() == 1)
				{
					bandwidth = 1.0;
				}
			}

			U inverse(size_t ind, const std::vector<std::vector<U>> &grids, const U mu, U val01, U in_bandwidth, const U lambda) const
			{
				const U lower_bound = grids[ind].front();
				const U upper_bound = grids[ind].back();

				const U es = upper_bound - lower_bound;
				const U min = compute_cdf(kernel_type, lower_bound, mu, in_bandwidth, lambda);
				const U max = compute_cdf(kernel_type, upper_bound, mu, in_bandwidth, lambda);

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
					f1 = (compute_cdf(kernel_type, p1, mu, in_bandwidth, lambda) - min)/(max - min);
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
				f2 = (compute_cdf(kernel_type, p2, mu, in_bandwidth, lambda) - min)/(max - min);

				U check = p1 + (val01 - f1)*(p2 - p1)/(f2 - f1);
				if(!std::isfinite(check))
				{
					size_t max_ind = 0;
					for(size_t j = 0; j != n - 1; j++)
					{
						p2 = lower_bound + es * (j + 1) / (n - 1);
						f2 = (compute_cdf(kernel_type, p2, mu, in_bandwidth, lambda)  - min)/(max - min);
						if(f2 < 1.0)
							max_ind = j;
						if(f2 > val01)
						{
							p1 = lower_bound + es * j / (n - 1);
							f1 = (compute_cdf(kernel_type, p1, mu, in_bandwidth, lambda)  - min)/(max - min);
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
				if(!std::isnormal(check))
				{
//        std::cout << check << std::endl;

//        std::cout.precision(20);
//        std::cout << std::scientific << std::numeric_limits<float>::denorm_min() << std::endl;
//        std::cout << std::scientific << std::numeric_limits<float>::lowest() << std::endl;
//        std::cout << std::scientific << p1 << " + (" << val01 << " - " << f1 << " )*( " << p2 << " - " << p1 << " )/( " << f2 << " - " << f1 << " ) " << std::endl;
//        std::cout << "fff" << std::endl;

//        U t1 = (val01 - f1)*(p2 - p1);
//        U t2 = t1/(f2 - f1);
//        U check = p1 + t2;
//        if(!std::isnormal(check))
//            std::cout << "ff" << std::endl;




					size_t max_ind = 0;
					for(size_t j = 0; j != n - 1; j++)
					{
						p2 = lower_bound + es * (j + 1) / (n - 1);
						f2 = (compute_cdf(kernel_type, p2, mu, in_bandwidth, lambda)  - min)/(max - min);
						if(f2 < 1.0)
							max_ind = j;
						if(f2 > val01)
						{
							p1 = lower_bound + es * j / (n - 1);
							f1 = (compute_cdf(kernel_type, p1, mu, in_bandwidth, lambda)  - min)/(max - min);
							return p1 + (val01 - f1)*(p2 - p1)/(f2 - f1);
						}
					}
					check = lower_bound + es * max_ind / (n - 1);
					if(std::isnormal(check))
					{
						return check;
					}
					else
					{
						if(max_ind > 0)
							return lower_bound + es * (max_ind - 1) / (n - 1);
						else
							return lower_bound + es * (max_ind + 1) / (n - 1);
					}
				}

//    std::cout << std::fixed << '\t' << time.elapsed_seconds() << std::endl;

				return p1 + (val01 - f1)*(p2 - p1)/(f2 - f1);
			}

		public:
			Qkde():sum(0), ssum(0), min(std::numeric_limits<U>::max()), max(std::numeric_limits<U>::min()) {}
			void set_kernel_type(size_t kt)
			{
//        if(kt >= 0 && kt < 5)
				kernel_type = kt;
//        else
//            throw std::logic_error("kernel type");
			}
			void set_sample(trie_based::NodeCount<T> *layer, size_t ind, const std::vector<std::vector<U>> &grids, const std::vector<U> &dx, U in_bandwidth, U lambda)
			{
				count = 0;
				size_t min_c = layer->count + 1, cc = 0;
				for(size_t i = 0; i != layer->children.size(); i++)
				{
					if(layer->children[i]->count < min_c)
					{
						min_c = layer->children[i]->count;
					}
				}
				min_c = min_c - 1;
				for(size_t i = 0; i != layer->children.size(); i++)
				{
					cc += layer->children[i]->count - min_c;
				}
				sample.resize(cc);

//        std::cout << sample.size() << std::endl;

//        sample.resize(layer->count);
				for(size_t i = 0; i != layer->children.size(); i++)
				{
//            U init_point = grids[ind][layer->children[i]->index] + dx[ind];
//            for(size_t j = 0; j != layer->children[i]->count; j++)
//            {
//                sum += sample;
//                ssum += std::pow(sample, 2.0);
//            }
//            size_t c = 0;
					for(size_t j = 1; j != layer->children[i]->count - min_c + 1; j++)
					{
						U temp = j/U(layer->children[i]->count - min_c + 1);
						//U inverse(size_t ind, const std::vector<std::vector<U>> &grids, const U mu, U val01, U in_bandwidth, const U lambda) const
						U rez = inverse(ind, grids, grids[ind][layer->children[i]->index] + dx[ind], temp, in_bandwidth, lambda);
						//U rez = grids[ind][layer->children[i]->index] + dx[ind] + myErfInv2(temp);
//                std::cout << temp << std::endl;
//                std::cin.get();
						if(rez > grids[ind].front() && rez < grids[ind].back())
						{
							sample[count] = rez;
							++count;
						}
						sum += rez;
						ssum += std::pow(rez, 2.0);

						min = rez < min ? rez : min;
						max = rez > max ? rez : max;
					}

//            sum += sample*layer->children[i]->count;
//            ssum += std::pow(sample, 2.0)*layer->children[i]->count;

//            min = sample < min ? sample : min;
//            max = sample > max ? sample : max;
				}

//        std::cout << sample.size() << std::endl;

				sample.resize(count);
				//count = layer->count;
				calculate_bandwidth(layer);
			}
			U cdf(U x, /*trie_based::NodeCount<T> *layer, size_t ind, const std::vector<std::vector<U>> &grids, const std::vector<U> &dx,*/ const U lambda = 1.0) const
			{
//        if(layer->children.size() < 100)
//        {

//        U d = 0;
//        for(size_t i = 0; i != layer->children.size(); i++)
//        {
//            d += compute_cdf(kernel_type, x, grids[ind][layer->children[i]->index] + dx[ind], bandwidth, lambda)*layer->children[i]->count;
//        }
//        return d/layer->count;

//        std::cout << sample.size() << std::endl;
				U d = 0;
				for(size_t i = 0; i != sample.size(); i++)
				{
					d += compute_cdf(kernel_type, x, sample[i], bandwidth, lambda);
				}
				return d/count;

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
//		using mveqf::ImplicitQuantile<T, U>::grids;
		std::vector<std::vector<U>> grids;
		using mveqf::ImplicitQuantile<T, U>::grid_number;
		using mveqf::ImplicitQuantile<T, U>::dx;
		using mveqf::ImplicitQuantile<T, U>::lb;
		using mveqf::ImplicitQuantile<T, U>::ub;

		typedef trie_based::Trie<trie_based::NodeCount<T>,T> trie_type;
		std::shared_ptr<trie_type> sample;

		using mveqf::ImplicitQuantile<T, U>::count_less;
		using mveqf::ImplicitQuantile<T, U>::quantile_transform;

		size_t kernel_type;
		std::vector<U> bandwidth;

		U kquantile_transform(trie_based::NodeCount<T> *layer, size_t ind, U val01, U in_bandwidth, const U lambda) const;
	public:
		ImplicitTrieKQuantile();
		ImplicitTrieKQuantile(std::vector<U> in_lb, std::vector<U> in_ub, std::vector<size_t> in_gridn, size_t kt);
		void set_kernel_type(size_t kt);
		void set_sample_shared(std::shared_ptr<trie_type> in_sample);
		void set_bandwidth(std::vector<U> in_bandwidth);
		void transform(const std::vector<U>& in01, std::vector<U>& out/*, const U lambda = 1.0*/) const;
		std::vector<U> transform(const std::vector<U>& in01/*, const U lambda = 1.0*/) const;
		std::vector<std::vector<U>> get_grid() const;
		std::vector<U> get_dx() const;
	};
	template <typename T, typename U>
	ImplicitTrieKQuantile<T, U>::ImplicitTrieKQuantile() : kernel_type(0)
	{
		grids.resize(grid_number.size());
		for(size_t i = 0; i != grids.size(); i++)
		{
			std::vector<U> grid(grid_number[i] + 1);
			U startp = lb[i];
			U endp = ub[i];
			U es = endp - startp;
			for(size_t j = 0; j != grid.size(); j++)
			{
				grid[j] = startp + j*es/U(grid_number[i]);
			}
			grids[i] = grid;
		}
	}
	template <typename T, typename U>
	ImplicitTrieKQuantile<T, U>::ImplicitTrieKQuantile(std::vector<U> in_lb,
	    std::vector<U> in_ub,
	    std::vector<size_t> in_gridn,
	    size_t kt) : mveqf::ImplicitQuantile<T, U>(in_lb, in_ub, in_gridn), kernel_type(kt)
	{
		grids.resize(grid_number.size());
		for(size_t i = 0; i != grids.size(); i++)
		{
			std::vector<U> grid(grid_number[i] + 1);
			U startp = lb[i];
			U endp = ub[i];
			U es = endp - startp;
			for(size_t j = 0; j != grid.size(); j++)
			{
				grid[j] = startp + j*es/U(grid_number[i]);
			}
			grids[i] = grid;
		}
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
	void ImplicitTrieKQuantile<T, U>::set_bandwidth(std::vector<U> in_bandwidth)
	{
		bandwidth = in_bandwidth;
	}

	template <typename T, typename U>
	void ImplicitTrieKQuantile<T, U>::transform(const std::vector<U>& in01, std::vector<U>& out/*, const U lambda*/) const
	{
		auto p = sample->root.get();
		for(size_t i = 0, k; i != in01.size(); i++)
		{
//        out[i] = kquantile_transform(p, i, in01[i], bandwidth[i], lambda);
//        k = quantile_transform(p, i, in01[i]).first;
			std::tie(k, out[i]) = quantile_transform(p, i, in01[i]);

//        k = 0;
//        U min_distance = std::abs(out[i] - grids[i][p->children.front()->index] + dx[i]);
//        for(size_t j = 1; j < p->children.size(); j++)
//        {
//            U temp = std::abs(out[i] - grids[i][p->children[j]->index] + dx[i]);
//            if(temp < min_distance)
//            {
//                k = j;
//                min_distance = temp;
//            }
//        }
			p = p->children[k].get();
		}
	}

	template <typename T, typename U>
	std::vector<U> ImplicitTrieKQuantile<T, U>::transform(const std::vector<U>& in01/*, const U lambda*/) const
	{
		std::vector<U> out(grid_number.size());
		transform(in01, out/*, lambda*/);
		return out;
	}

	template <typename T, typename U>
	U ImplicitTrieKQuantile<T, U>::kquantile_transform(trie_based::NodeCount<T> *layer, size_t ind, U val01, U in_bandwidth, const U lambda) const
	{

//    timer::Timer time;


		kquantile::Qkde<T, U> obj;
		obj.set_kernel_type(kernel_type);
		obj.set_sample(layer, ind, grids, dx, in_bandwidth, lambda);

//    std::cout << std::fixed << '\t' << time.elapsed_seconds() << std::endl;
//    time.reset();

		const U lower_bound = grids[ind].front();
		const U upper_bound = grids[ind].back();

		const U es = upper_bound - lower_bound;
		const U min = obj.cdf(lower_bound,/* layer, ind, grids, dx,*/ lambda);
		const U max = obj.cdf(upper_bound,/* layer, ind, grids, dx,*/ lambda);

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
			f1 = (obj.cdf(p1, /*layer, ind, grids, dx,*/ lambda) - min)/(max - min);
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
		f2 = (obj.cdf(p2, /*layer, ind, grids, dx,*/ lambda) - min)/(max - min);

		U check = p1 + (val01 - f1)*(p2 - p1)/(f2 - f1);
		if(!std::isfinite(check))
		{
			size_t max_ind = 0;
			for(size_t j = 0; j != n - 1; j++)
			{
				p2 = lower_bound + es * (j + 1) / (n - 1);
				f2 = (obj.cdf(p2, /*layer, ind, grids, dx,*/ lambda) - min)/(max - min);
				if(f2 < 1.0)
					max_ind = j;
				if(f2 > val01)
				{
					p1 = lower_bound + es * j / (n - 1);
					f1 = (obj.cdf(p1, /*layer, ind, grids, dx,*/ lambda) - min)/(max - min);
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
		if(!std::isnormal(check))
		{
//        std::cout << check << std::endl;

//        std::cout.precision(20);
//        std::cout << std::scientific << std::numeric_limits<float>::denorm_min() << std::endl;
//        std::cout << std::scientific << std::numeric_limits<float>::lowest() << std::endl;
//        std::cout << std::scientific << p1 << " + (" << val01 << " - " << f1 << " )*( " << p2 << " - " << p1 << " )/( " << f2 << " - " << f1 << " ) " << std::endl;
//        std::cout << "fff" << std::endl;

//        U t1 = (val01 - f1)*(p2 - p1);
//        U t2 = t1/(f2 - f1);
//        U check = p1 + t2;
//        if(!std::isnormal(check))
//            std::cout << "ff" << std::endl;




			size_t max_ind = 0;
			for(size_t j = 0; j != n - 1; j++)
			{
				p2 = lower_bound + es * (j + 1) / (n - 1);
				f2 = (obj.cdf(p2, /*layer, ind, grids, dx,*/ lambda) - min)/(max - min);
				if(f2 < 1.0)
					max_ind = j;
				if(f2 > val01)
				{
					p1 = lower_bound + es * j / (n - 1);
					f1 = (obj.cdf(p1, /*layer, ind, grids, dx,*/ lambda) - min)/(max - min);
					return p1 + (val01 - f1)*(p2 - p1)/(f2 - f1);
				}
			}
			check = lower_bound + es * max_ind / (n - 1);
			if(std::isnormal(check))
			{
				return check;
			}
			else
			{
				if(max_ind > 0)
					return lower_bound + es * (max_ind - 1) / (n - 1);
				else
					return lower_bound + es * (max_ind + 1) / (n - 1);
			}
		}

//    std::cout << std::fixed << '\t' << time.elapsed_seconds() << std::endl;

		return p1 + (val01 - f1)*(p2 - p1)/(f2 - f1);
	}

	template <typename T, typename U>
	std::vector<std::vector<U>> ImplicitTrieKQuantile<T, U>::get_grid() const
	{
		return grids;
	}

	template <typename T, typename U>
	std::vector<U> ImplicitTrieKQuantile<T, U>::get_dx() const
	{
		return dx;
	}

}

#endif
