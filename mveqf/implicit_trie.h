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
#ifndef IMPLICIT_TRIE_H
#define IMPLICIT_TRIE_H

#include <mveqf/trie_node.h>
#include <mveqf/implicit.h>

namespace mveqf
{
	template <typename TIndex, typename TFloat>
	class ImplicitTrieQuantile : public ImplicitQuantile<TIndex, TFloat>
	{
	protected:
//		using ImplicitQuantile<TIndex, TFloat>::grids;
//		using ImplicitQuantile<TIndex, TFloat>::dx;
		using ImplicitQuantile<TIndex, TFloat>::grid_number;
		using ImplicitQuantile<TIndex, TFloat>::lb;
		using ImplicitQuantile<TIndex, TFloat>::ub;

		typedef Trie<NodeCount<TIndex>,TIndex> trie_type;
		std::shared_ptr<trie_type> sample;

//		using ImplicitQuantile<TIndex, TFloat>::count_less;
		using ImplicitQuantile<TIndex, TFloat>::quantile_transform;
//		const std::vector<size_t> weights;
	public:
		ImplicitTrieQuantile() = default;
		ImplicitTrieQuantile(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		ImplicitTrieQuantile(const ImplicitTrieQuantile&) = delete;
		ImplicitTrieQuantile& operator=(const ImplicitTrieQuantile&) = delete;
		void set_sample_shared(std::shared_ptr<trie_type> in_sample);
		void set_sample(const std::vector<std::vector<TIndex>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights) override;
		void transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const override;
		void transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const override;
		using ImplicitQuantile<TIndex, TFloat>::get_the_closest_grid_node_to_the_value;
	};

	template <typename TIndex, typename TFloat>
	ImplicitTrieQuantile<TIndex, TFloat>::ImplicitTrieQuantile(std::vector<TFloat> in_lb,
	    std::vector<TFloat> in_ub,
	    std::vector<size_t> in_gridn) : ImplicitQuantile<TIndex, TFloat>(in_lb, in_ub, in_gridn)
	{
	}
	template <typename TIndex, typename TFloat>
	void ImplicitTrieQuantile<TIndex, TFloat>::set_sample_shared(std::shared_ptr<trie_type> in_sample)
	{
		sample = std::move(in_sample);
	}

	template <typename TIndex, typename TFloat>
	void ImplicitTrieQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const
	{
		auto p = sample->root;
		for(size_t i = 0, k; i != in01.size(); i++)
		{
			std::tie(k, out[i]) = quantile_transform(p, i, in01[i]);
			p = p->children[k];
		}
	}

	template <typename TIndex, typename TFloat>
	void ImplicitTrieQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const
	{
		auto p = sample->root;
		for(size_t i = 0; i != in01.size(); i++)
		{
			//std::tie(k, out[i]) = quantile_transform(p, i, in01[i]);
			auto [k, result] = quantile_transform(p, i, in01[i]);
			out[i] = p->children[k]->index;
			p = p->children[k];
		}
	}

	template <typename TIndex, typename TFloat>
	void ImplicitTrieQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TIndex>> &in_sample)
	{
		ImplicitQuantile<TIndex, TFloat>::set_sample(in_sample);
	}

	template <typename TIndex, typename TFloat>
	void ImplicitTrieQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample)
	{
		sample = std::make_shared<trie_type>();
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
			sample->insert(temp, 1);
		}
		//sample->fill_tree_count();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitTrieQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights)
	{
		sample = std::make_shared<trie_type>();
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
			sample->insert(temp, weights[i]);
		}
		//sample->fill_tree_count();
	}
}

#endif