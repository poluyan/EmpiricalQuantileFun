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
#ifndef SMOT_H
#define SMOT_H

#include <mveqf/trie_based.h>
#include <mveqf/explicit.h>
#include <mveqf/implicit.h>

#include <iostream>
#include <string>
#include <unordered_map>

namespace mveqf
{
	namespace ot
	{
		class Hasher
		{
		public:
			size_t operator()(const std::vector<size_t>& key) const
			{
				std::size_t seed = key.size();
				for(auto& i : key)
				{
					seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
				}
				return seed;
			}
		};
		class EqualFn
		{
		public:
			bool operator()(const std::vector<size_t>& l, const std::vector<size_t>& r) const
			{
				for(size_t i = 0; i != l.size(); i++)
				{
					if(l[i] != r[i])
						return false;
				}
				return true;
			}
		};
		template <typename TFloat>
		class SDOT
		{
			std::shared_ptr<mveqf::Quantile<size_t, TFloat>> qf;
			const std::vector<std::vector<TFloat>>& init_dots;
			const std::vector<size_t>& weights;
			std::unordered_map<std::vector<size_t>, size_t, Hasher, EqualFn> index;
			std::vector<TFloat> lb;
			std::vector<TFloat> ub;
		public:
			SDOT(const std::vector<std::vector<TFloat>> &in_dots, const std::vector<size_t> &in_weights) : init_dots(in_dots), weights(in_weights)
			{
				lb = std::vector<TFloat>(init_dots.front().size(), std::numeric_limits<TFloat>::max());
				ub = std::vector<TFloat>(init_dots.front().size(), std::numeric_limits<TFloat>::min());

				if(std::adjacent_find(weights.begin(), weights.end(), std::not_equal_to<>()) == weights.end() && weights.front() == 1)
					initialize(false);
				else
					initialize(true);
			}

			void initialize(bool weights_flag)
			{
				for(size_t i = 0; i != init_dots.size(); i++)
				{
					for(size_t j = 0; j != init_dots[i].size(); j++)
					{
						if(init_dots[i][j] < lb[j])
							lb[j] = init_dots[i][j];
						if(init_dots[i][j] > ub[j])
							ub[j] = init_dots[i][j];
					}
				}

				if(weights_flag)
				{
					//std::make_shared<mveqf::ImplicitTrieQuantile<size_t, TFloat>>()
					qf = std::make_shared<mveqf::ExplicitQuantile<size_t, TFloat>>();
				}
				else
				{
					qf = std::make_shared<mveqf::ImplicitQuantile<size_t, TFloat>>();
				}
				qf->set_grid_from_sample(lb, ub, init_dots);
				qf->set_sample(init_dots, weights);
				auto grid_number = qf->get_grid_number();
				for(size_t i = 0; i != init_dots.size(); i++)
				{
					std::vector<size_t> point(lb.size());
					for(size_t j = 0; j != point.size(); j++)
					{
						point[j] = qf->get_the_closest_grid_node_to_the_value(lb[j], ub[j], grid_number[j], init_dots[i][j]);
					}
					index[point] = i;
				}
			}

			SDOT(const std::vector<std::vector<TFloat>> &&) = delete;  // prevents rvalue binding
//			void set_sample(const std::vector<std::vector<TFloat>> &sample, const std::vector<size_t> &weights);
//			void set_sample(const std::vector<std::vector<TFloat>> &sample, const std::vector<size_t> &weights, std::vector<TFloat> lb, std::vector<TFloat> ub);
			std::vector<size_t> transform(const std::vector<TFloat> &val01)
			{
				std::vector<size_t> res(val01.size());
				qf->transform(val01, res);
				//return init_dots[index[res]];
				return res;
			}
			std::vector<TFloat> transform_float(const std::vector<TFloat> &val01)
			{
				std::vector<TFloat> res(val01.size());
				qf->transform(val01, res);
				return res;
			}
			size_t transform_index(const std::vector<TFloat> &val01)
			{
				std::vector<size_t> res(val01.size());
				qf->transform(val01, res);
				return index[res];
			}
		};
	}
}

#endif
