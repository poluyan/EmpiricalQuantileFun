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
#ifndef EXPLICIT_H
#define EXPLICIT_H

#include <quantile.h>

namespace mveqf
{
	template <typename TIndex, typename TFloat>
	class ExplicitQuantile : public Quantile<TIndex, TFloat>
	{
	protected:
		//using Quantile<TIndex, TFloat>::grids;
		using Quantile<TIndex, TFloat>::grid_number;

		using Quantile<TIndex, TFloat>::lb;
		using Quantile<TIndex, TFloat>::ub;
		using Quantile<TIndex, TFloat>::dx;

		using Quantile<TIndex, TFloat>::get_grid_value;
		using Quantile<TIndex, TFloat>::get_the_closest_grid_node_to_the_value;

		typedef std::vector<std::vector<TFloat>> sample_type;
		std::shared_ptr<sample_type> sample;

		size_t count_less(const std::vector<TFloat> &layer, TFloat target) const;
		std::pair<size_t, TFloat> quantile_transform(const std::vector<TFloat> &layer, size_t ind, TFloat val01) const;

		size_t get_lower_bound(size_t ind, const TFloat &value) const;
	public:
		ExplicitQuantile() = default;
		ExplicitQuantile(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		ExplicitQuantile(const ExplicitQuantile&) = delete;
		ExplicitQuantile& operator=(const ExplicitQuantile&) = delete;
		using Quantile<TIndex, TFloat>::set_grid_and_gridn;
		void set_sample(const std::vector<std::vector<TIndex>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample) override;
		void set_sample_shared(std::shared_ptr<sample_type> in_sample);
		void transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const override;
		void transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const override;
	};

	template <typename TIndex, typename TFloat>
	ExplicitQuantile<TIndex, TFloat>::ExplicitQuantile(std::vector<TFloat> in_lb,
	    std::vector<TFloat> in_ub,
	    std::vector<size_t> in_gridn): Quantile<TIndex, TFloat>(in_lb, in_ub, in_gridn)
	{
	}

	template <typename TIndex, typename TFloat>
	void ExplicitQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TIndex>> &in_sample)
	{
		sample = std::make_shared< std::vector<std::vector<TFloat>> >();
		for(size_t i = 0; i != in_sample.size(); ++i)
		{
			std::vector<TFloat> temp;
			for(size_t j = 0; j != in_sample[i].size(); ++j)
			{
				//temp.push_back(grids[j][in_sample[i][j]] + dx[j]);
				temp.push_back(get_grid_value(j, in_sample[i][j]) + dx[j]);
			}
			sample->push_back(temp);
		}
	}

	template <typename TIndex, typename TFloat>
	void ExplicitQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample)
	{
		sample = std::make_shared< std::vector<std::vector<TFloat>> >();
		for(size_t i = 0; i != in_sample.size(); ++i)
		{
			std::vector<TFloat> temp;
			for(size_t j = 0; j != in_sample[i].size(); ++j)
			{
				TIndex index = get_the_closest_grid_node_to_the_value(lb[j], ub[j], grid_number[j], in_sample[i][j]);
				temp.push_back(get_grid_value(j, index) + dx[j]);
			}
			sample->push_back(temp);
		}
	}

	template <typename TIndex, typename TFloat>
	void ExplicitQuantile<TIndex, TFloat>::set_sample_shared(std::shared_ptr<sample_type> in_sample)
	{
		sample = std::move(in_sample);
	}

	template <typename TIndex, typename TFloat>
	size_t ExplicitQuantile<TIndex, TFloat>::get_lower_bound(size_t ind, const TFloat &value) const
	{
		size_t it, first = 0;
		int count = grid_number[ind] + 1, step;

		while(count > 0)
		{
			it = first;
			step = count / 2;
			it += step;

			if(get_grid_value(ind, it) < value)
			{
				first = ++it;
				count -= step + 1;
			}
			else
				count = step;
		}
		return first;
	}

	template <typename TIndex, typename TFloat>
	void ExplicitQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const
	{
		std::vector<size_t> m(grid_number.size());
		for(size_t i = 0, g = in01.size(); i != g; i++)
		{
			std::vector<TFloat> row(sample->size());
			size_t index = 0;
			for(size_t j = 0, n = sample->size(); j != n; j++)
			{
				bool flag = true;
				for(size_t k = 0; k != i; k++)
				{
					//if(!((*sample)[j][k] > grids[k][m[k]] && (*sample)[j][k] < grids[k][m[k] + 1]))
					if(!((*sample)[j][k] > get_grid_value(k, m[k]) && (*sample)[j][k] < get_grid_value(k, m[k] + 1)))
					{
						flag = false;
						break;
					}
				}
				if(flag)
				{
					row[index] = (*sample)[j][i];
					++index;
				}
			}
			row.resize(index);

//			std::tie(m[i], out[i]) = quantile_transform(row, i, in01[i]);
			auto [k, res] = quantile_transform(row, i, in01[i]);
			out[i] = res;
			m[i] = k;
		}
	}

	template <typename TIndex, typename TFloat>
	void ExplicitQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const
	{
		std::vector<size_t> m(grid_number.size());
		for(size_t i = 0, g = in01.size(); i != g; i++)
		{
			std::vector<TFloat> row(sample->size());
			size_t index = 0;
			for(size_t j = 0, n = sample->size(); j != n; j++)
			{
				bool flag = true;
				for(size_t k = 0; k != i; k++)
				{
					//if(!((*sample)[j][k] > grids[k][m[k]] && (*sample)[j][k] < grids[k][m[k] + 1]))
					if(!((*sample)[j][k] > get_grid_value(k, m[k]) && (*sample)[j][k] < get_grid_value(k, m[k] + 1)))
					{
						flag = false;
						break;
					}
				}
				if(flag)
				{
					row[index] = (*sample)[j][i];
					++index;
				}
			}
			row.resize(index);

//			std::tie(m[i], out[i]) = quantile_transform(row, i, in01[i]);
			auto [k, res] = quantile_transform(row, i, in01[i]);
			out[i] = k;
			m[i] = k;
		}
	}

	template <typename TIndex, typename TFloat>
	size_t ExplicitQuantile<TIndex, TFloat>::count_less(const std::vector<TFloat> &layer, TFloat target) const
	{
		return std::count_if(layer.begin(), layer.end(), [&](const TFloat &v)
		{
			return v < target;
		});
	}

	template <typename TIndex, typename TFloat>
	std::pair<size_t, TFloat> ExplicitQuantile<TIndex, TFloat>::quantile_transform(const std::vector<TFloat> &layer, size_t ind, TFloat val01) const
	{
		size_t count = grid_number[ind], step, c1 = 0, c2 = 0, m = 0;
		TFloat f1 = 0.0, f2 = 0.0, n = layer.size();
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

			c1 = count_less(layer, get_grid_value(ind, m));
			f1 = c1/n;

			if(f1 < val01)
			{
				c2 = count_less(layer, get_grid_value(ind, m + 1));
				f2 = c2/n;

				if(val01 < f2)
					break;

				first = ++it;
				count -= step + 1;
			}
			else
				count = step;
		}

		if(count == 0)
		{
			c2 = count_less(layer, get_grid_value(ind, m + 1));
			f2 = c2/n;
		}

		if(c1 == c2)
		{
			if(c1 == 0)
			{
				auto min_val = *std::min_element(layer.begin(), layer.end()) - 2.0*dx[ind];
				//auto lb_min = std::lower_bound(grids[ind].begin(), grids[ind].end(), min_val);
				//size_t min_ind = std::distance(grids[ind].begin(), lb_min);
				size_t min_ind = get_lower_bound(ind, min_val);
				return std::make_pair(min_ind, get_grid_value(ind, min_ind) + 2.0*val01*dx[ind]);
			}
			if(c1 == layer.size())
			{
				auto max_val = *std::max_element(layer.begin(), layer.end()) - 2.0*dx[ind];
				//auto lb_max = std::lower_bound(grids[ind].begin(), grids[ind].end(), max_val);
				//size_t max_ind = std::distance(grids[ind].begin(), lb_max);
				size_t max_ind = get_lower_bound(ind, max_val);
				return std::make_pair(max_ind, get_grid_value(ind, max_ind) + 2.0*val01*dx[ind]);
			}

			TFloat target = get_grid_value(ind, m);
			TFloat diff = std::numeric_limits<TFloat>::max();
			size_t index = 0;
			for(size_t i = 1; i != layer.size(); ++i)
			{
				TFloat curr = std::abs(layer[i] - target);
				if(diff > curr)
				{
					diff = curr;
					index = i;
				}
			}

			//auto lb = std::lower_bound(grids[ind].begin(), grids[ind].end(), layer[index] - 2.0*dx[ind]);
			//size_t lb_ind = std::distance(grids[ind].begin(), lb);
			size_t lb_ind = get_lower_bound(ind, layer[index] - 2.0*dx[ind]);
			return std::make_pair(lb_ind, get_grid_value(ind, lb_ind) + 2.0*val01*dx[ind]);
		}
		//return std::make_pair(m, grids[ind][m] + (val01 - f1) * (grids[ind][m + 1] - grids[ind][m]) / (f2 - f1));
		return std::make_pair(m, get_grid_value(ind, m) + (val01 - f1) * (get_grid_value(ind, m + 1) - get_grid_value(ind, m)) / (f2 - f1));
	}
}

#endif
