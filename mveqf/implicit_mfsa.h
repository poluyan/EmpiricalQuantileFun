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
#ifndef IMPLICIT_MFSA_H
#define IMPLICIT_MFSA_H

#include <mveqf/quantile.h>
#include <mveqf/mfsa.h>

namespace mveqf
{
	template <typename TIndex, typename TFloat>
	class ImplicitQuantileMFSA : public Quantile<TIndex, TFloat>
	{
	protected:
		typedef mveqf::mfsa::MFSA<TIndex> sample_type;
		std::shared_ptr<sample_type> sample;

		//using Quantile<TIndex, TFloat>::grids;
		using Quantile<TIndex, TFloat>::grid_number;
		using Quantile<TIndex, TFloat>::dx;
		using Quantile<TIndex, TFloat>::lb;
		using Quantile<TIndex, TFloat>::ub;

		using Quantile<TIndex, TFloat>::get_grid_value;

		std::pair<size_t, size_t> count_less(mfsa::Node<TIndex> *layer, const size_t &r) const;
		std::pair<size_t, TFloat> quantile_transform(mfsa::Node<TIndex> *layer, size_t ind, TFloat val01) const;
	public:
		ImplicitQuantileMFSA() = default;
		ImplicitQuantileMFSA(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		ImplicitQuantileMFSA(const ImplicitQuantileMFSA&) = delete;
		ImplicitQuantileMFSA& operator=(const ImplicitQuantileMFSA&) = delete;
		void set_sample(const std::vector<std::vector<TIndex>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights) override;
		void set_sample_and_fill_count(const std::vector<std::vector<TIndex>> &in_sample);
		void set_sample_shared_and_fill_count(std::shared_ptr<sample_type> in_sample);
		void set_sample_shared(std::shared_ptr<sample_type> in_sample);
		void transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const override;
		void transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const override;
		size_t get_node_count() const;
		size_t get_link_count() const;
		using Quantile<TIndex, TFloat>::get_the_closest_grid_node_to_the_value;
		using Quantile<TIndex, TFloat>::get_real_node_values;
		~ImplicitQuantileMFSA();
	};

	template <typename TIndex, typename TFloat>
	ImplicitQuantileMFSA<TIndex, TFloat>::ImplicitQuantileMFSA(std::vector<TFloat> in_lb,
	    std::vector<TFloat> in_ub,
	    std::vector<size_t> in_gridn) : Quantile<TIndex, TFloat>(in_lb, in_ub, in_gridn)
	{}

	template <typename TIndex, typename TFloat>
	ImplicitQuantileMFSA<TIndex, TFloat>::~ImplicitQuantileMFSA()
	{}

	template <typename TIndex, typename TFloat>
	size_t ImplicitQuantileMFSA<TIndex, TFloat>::get_node_count() const
	{
		return sample->get_node_count();
	}

	template <typename TIndex, typename TFloat>
	size_t ImplicitQuantileMFSA<TIndex, TFloat>::get_link_count() const
	{
		return sample->get_link_count();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantileMFSA<TIndex, TFloat>::set_sample_and_fill_count(const std::vector<std::vector<TIndex>> &in_sample)
	{
		sample = std::make_shared<sample_type>();
		sample->set_dimension(grid_number.size());
		for(const auto & i : in_sample)
			sample->insert(i);
		sample->fill_tree_count();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantileMFSA<TIndex, TFloat>::set_sample(const std::vector<std::vector<TIndex>> &in_sample)
	{
		set_sample_and_fill_count(in_sample);
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantileMFSA<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample)
	{
		sample = std::make_shared<sample_type>();
		sample->set_dimension(grid_number.size());
		for(size_t i = 0; i != in_sample.size(); ++i)
		{
			std::vector<TIndex> temp(in_sample[i].size());
			for(size_t j = 0; j != in_sample[i].size(); ++j)
			{
				temp[j] = get_the_closest_grid_node_to_the_value(lb[j], ub[j], grid_number[j], in_sample[i][j]);
			}
//			for(size_t j = 0; j != temp.size(); j++)
//			{
//				std::cout << temp[j] << ' ';
//			}
//			std::cout << std::endl;
			sample->insert(temp);
		}
		sample->fill_tree_count();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantileMFSA<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights)
	{
		set_sample(in_sample);
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantileMFSA<TIndex, TFloat>::set_sample_shared_and_fill_count(std::shared_ptr<sample_type> in_sample)
	{
		sample = std::move(in_sample);
		sample->fill_tree_count();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitQuantileMFSA<TIndex, TFloat>::set_sample_shared(std::shared_ptr<sample_type> in_sample)
	{
		sample = std::move(in_sample);
	}

	template <typename TIndex, typename TFloat>
	std::pair<size_t, size_t> ImplicitQuantileMFSA<TIndex, TFloat>::count_less(mfsa::Node<TIndex> *layer, const size_t &r) const
	{
		std::pair<size_t, size_t> res;
		for(const auto &i : layer->children)
		{
			size_t j = static_cast<size_t>(i.first);
			if(j < r + 1)
			{
				res.second += i.second->count;
				if(j < r)
				{
					res.first += i.second->count;
				}
			}
		}
		return res;
	}
	template <typename TIndex, typename TFloat>
	void ImplicitQuantileMFSA<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const
	{
		auto p = sample->root.get();
		for(size_t i = 0, k; i != in01.size(); ++i)
		{
			std::tie(k, out[i]) = quantile_transform(p, i, in01[i]);
			p = p->children[k].second.get();
		}
	}
	template <typename TIndex, typename TFloat>
	void ImplicitQuantileMFSA<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const
	{
		auto p = sample->root.get();
		for(size_t i = 0; i != in01.size(); ++i)
		{
			//std::tie(k, out[i]) = quantile_transform(p, i, in01[i]);
			auto [k, result] = quantile_transform(p, i, in01[i]);
			out[i] = p->children[k].first;
			p = p->children[k].second.get();
		}
	}
	template <typename TIndex, typename TFloat>
	std::pair<size_t, TFloat> ImplicitQuantileMFSA<TIndex, TFloat>::quantile_transform(mveqf::mfsa::Node<TIndex> *layer, size_t ind, TFloat val01) const
	{
		size_t m = 0, count = grid_number[ind], step, a = 0, b = 0;
		TFloat x = 0.0, y = 0.0, p = static_cast<TFloat>(layer->count);
		//auto first = grids[ind].begin();
		//auto it = grids[ind].begin();
		size_t it = 0, first = 0;

		while(count > 0)
		{
			it = first;
			step = count / 2;
			it += step;
			m = it;
			//std::advance(it, step);
			//m = std::distance(grids[ind].begin(), it);

			std::tie(a, b) = count_less(layer, m);
			x = static_cast<TFloat>(a)/p;

			if(x < val01)
			{
				y = static_cast<TFloat>(b)/p;
				if(val01 < y)
					break;

				first = ++it;
				count -= step + 1;
			}
			else
				count = step;
		}
		if(count == 0)
		{
			y = static_cast<TFloat>(b)/p;
		}
		if(a == b)
		{
			if(a == 0)
			{
				auto min_val_it = std::min_element(layer->children.begin(), layer->children.end(),
				                                   [](const auto &l,
				                                      const auto &r)
				{
					return l.first < r.first;
				});
				size_t min_ind = std::distance(layer->children.begin(), min_val_it);
				//return std::make_pair(min_ind, grids[ind][layer->children[min_ind]->index] + 2.0*val01*dx[ind]);
				return std::make_pair(min_ind, get_grid_value(ind, layer->children[min_ind].first) + 2.0*val01*dx[ind]);
			}
			if(a == layer->count)
			{
				auto max_val_it = std::max_element(layer->children.begin(), layer->children.end(),
				                                   [](const auto &l,
				                                      const auto &r)
				{
					return l.first < r.first;
				});
				size_t max_ind = std::distance(layer->children.begin(), max_val_it);
				//return std::make_pair(max_ind, grids[ind][layer->children[max_ind]->index] + 2.0*val01*dx[ind]);
				return std::make_pair(max_ind, get_grid_value(ind, layer->children[max_ind].first) + 2.0*val01*dx[ind]);
			}
			int diff = std::numeric_limits<int>::max();
			size_t index = 0;
			int min_ind = static_cast<int>(layer->children[index].first);
			for(size_t i = 1; i != layer->children.size(); ++i)
			{
				int t = static_cast<int>(layer->children[i].first);
				int curr = std::abs(t - static_cast<int>(m));
				if(diff > curr)
				{
					diff = curr;
					index = i;
					min_ind = t;
				}
				else if(diff == curr)
				{
					if(min_ind > t)
					{
						min_ind = t;
						index = i;
					}
				}
			}
			//return std::make_pair(index, grids[ind][layer->children[index]->index] + 2.0*val01*dx[ind]);
			return std::make_pair(index, get_grid_value(ind, layer->children[index].first) + 2.0*val01*dx[ind]);
		}
		size_t index = 0;
		TIndex target = static_cast<TIndex>(m);
		for(size_t j = 1; j < layer->children.size(); j++)
		{
			if(layer->children[j].first == target)
			{
				index = j;
				break;
			}
		}
		//return std::make_pair(index, grids[ind][m] + (val01 - x) * (grids[ind][m + 1] - grids[ind][m]) / (y - x));
		return std::make_pair(index, get_grid_value(ind, m) + (val01 - x) * (get_grid_value(ind, m + 1) - get_grid_value(ind, m)) / (y - x));
	}
}

#endif
