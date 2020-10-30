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
#ifndef EXPERIMENTAL_H
#define EXPERIMENTAL_H

#include <mveqf/trie_node.h>
#include <mveqf/quantile.h>

namespace mveqf
{
	namespace trie_based
	{
		template <typename I>
		class TrieLayer
		{
		public:
			size_t dimension;
			std::vector< std::multimap<I,std::vector<I>> > layers;
			typename std::multimap<I,std::vector<I>>::iterator root;
			bool sorted;
		public:
			TrieLayer();
			TrieLayer(size_t dim);
			void set_dimension(size_t dim);
			size_t get_dimension() const;
			bool empty() const;
			void insert(const std::vector<I> &key);
			bool search(const std::vector<I> &key) const;
			void sort();
//    void print() const;
		};

		template <typename I>
		TrieLayer<I>::TrieLayer() {}

		template <typename I>
		TrieLayer<I>::TrieLayer(size_t dim): dimension(dim)
		{
			layers.resize(dimension);
			layers.front().insert(std::pair<I, std::vector<I>>(0, std::vector<I>()));
			root = layers.front().find(0);
		}
		template <typename I>
		void TrieLayer<I>::set_dimension(size_t dim)
		{
			dimension = dim;
			layers.resize(dimension);
			layers.front().insert(std::pair<I, std::vector<I>>(0, std::vector<I>()));
			root = layers.front().find(0);
		}
		template <typename I>
		size_t TrieLayer<I>::get_dimension() const
		{
			return dimension;
		}
		template <typename I>
		bool TrieLayer<I>::empty() const
		{
			return root->second.empty();
		}
		template <typename I>
		void TrieLayer<I>::insert(const std::vector<I> &key)
		{
			sorted = false;

			auto pos = std::find(root->second.begin(), root->second.end(), key.front());
			if(pos == root->second.end())
				root->second.push_back(key.front());

			for(size_t i = 0; i != dimension - 1; i++)
			{
				auto it = layers[i+1].find(key[i]);
				if(it == layers[i+1].end())
				{
					layers[i+1].insert(std::pair<I, std::vector<I>>(key[i],std::vector<I>(1, key[i+1])));
				}
				else
				{
					auto pos = std::find(it->second.begin(), it->second.end(), key[i + 1]);
					if(pos == it->second.end())
						it->second.push_back(key[i + 1]);
				}
			}
		}
		template <typename I>
		void TrieLayer<I>::sort()
		{
			for(size_t i = 0; i != layers.size(); i++)
			{
				for(auto it = layers[i].begin(); it != layers[i].end(); ++it)
				{
					std::sort(it->second.begin(), it->second.end());
				}
			}
			sorted = true;
		}
//template <typename I>
//void TrieLayer<I>::print() const
//{
//    for(size_t i = 0; i != layers.size(); i++)
//    {
//        for(auto it = layers[i].begin(); it != layers[i].end(); ++it)
//        {
//            std::cout << int(it->first) << " => ";
//            for(size_t j = 0; j != it->second.size(); j++)
//            {
//                std::cout << int(it->second[j]) << ' ';
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//    }
//}
		template <typename I>
		bool TrieLayer<I>::search(const std::vector<I> &key) const
		{
			if(sorted)
			{
				auto pos = std::lower_bound(root->second.begin(), root->second.end(), key.front());
				if(pos == root->second.end())
					return false;
				if(*pos != key.front())
					return false;

				for(size_t i = 0; i != dimension - 1; i++)
				{
					auto it = layers[i + 1].find(key[i]);
					if(it == layers[i + 1].end())
					{
						return false;
					}
					else
					{
						auto pos = std::lower_bound(it->second.begin(), it->second.end(), key[i + 1]);
						if(pos == it->second.end())
							return false;
						if(*pos != key[i + 1])
							return false;
					}
				}
				return true;
			}
			else
			{
				auto pos = std::find(root->second.begin(), root->second.end(), key.front());
				if(pos == root->second.end())
					return false;
				for(size_t i = 0; i != dimension - 1; i++)
				{
					auto it = layers[i + 1].find(key[i]);
					if(it == layers[i + 1].end())
					{
						return false;
					}
					else
					{
						auto pos = std::find(it->second.begin(), it->second.end(), key[i + 1]);
						if(pos == it->second.end())
							return false;
					}
				}
				return true;
			}
		}

		struct invect
		{
			char vname;
			int index;
			size_t count;
			invect(char _vname, int _index, int _count)
			{
				vname = _vname;
				index = _index;
				count = _count;
			}
		};

		template <typename I>
		class Graph
		{
		public:
			std::unordered_map< char, std::vector<invect> > layers;
			typename std::unordered_map<char, std::vector<invect>>::iterator root;
		public:
			Graph();
			void print() const;
		};

		template <typename I>
		Graph<I>::Graph()
		{
			std::vector<invect> y, a,b,c,d,e, f,g,h,j,k,l, z;
			y.push_back(invect('a',0,4));
			y.push_back(invect('b',1,1));
			y.push_back(invect('e',2,4));
			y.push_back(invect('d',3,2));
			y.push_back(invect('c',4,3));
			layers.insert(std::make_pair('y',y));

			a.push_back(invect('g',0,1));
			a.push_back(invect('h',2,1));
			a.push_back(invect('f',3,2));
			layers.insert(std::make_pair('a',a));

			b.push_back(invect('h',0,1));
			layers.insert(std::make_pair('b',b));

			c.push_back(invect('h',0,1));
			c.push_back(invect('h',3,1));
			c.push_back(invect('h',4,1));
			layers.insert(std::make_pair('c',c));

			d.push_back(invect('j',0,1));
			d.push_back(invect('h',3,1));
			layers.insert(std::make_pair('d',d));

			e.push_back(invect('k',0,2));
			e.push_back(invect('l',1,1));
			e.push_back(invect('l',2,1));
			layers.insert(std::make_pair('e',e));


			f.push_back(invect('z',2,1));
			f.push_back(invect('z',3,1));
			layers.insert(std::make_pair('f',f));

			g.push_back(invect('z',1,1));
			layers.insert(std::make_pair('g',g));

			h.push_back(invect('z',0,1));
			layers.insert(std::make_pair('h',h));

			j.push_back(invect('z',2,1));
			layers.insert(std::make_pair('j',j));

			k.push_back(invect('z',0,1));
			k.push_back(invect('z',4,1));
			layers.insert(std::make_pair('k',k));

			l.push_back(invect('z',4,1));
			layers.insert(std::make_pair('l',l));

			layers.insert(std::make_pair('z',z));

			root = layers.find('y');

			/*std::vector<invect> y, a,b,c,d,e, f,g,h,j,k,l, z;
			y.push_back(invect('a',0,1));
			y.push_back(invect('a',1,1));
			y.push_back(invect('a',2,1));
			y.push_back(invect('a',3,1));
			y.push_back(invect('a',4,1));
			layers.insert(std::make_pair('y',y));

			a.push_back(invect('b',0,1));
			a.push_back(invect('b',1,1));
			a.push_back(invect('b',2,1));
			a.push_back(invect('b',3,1));
			a.push_back(invect('b',4,1));
			layers.insert(std::make_pair('a',a));

			b.push_back(invect('z',0,1));
			b.push_back(invect('z',1,1));
			b.push_back(invect('z',2,1));
			b.push_back(invect('z',3,1));
			b.push_back(invect('z',4,1));
			layers.insert(std::make_pair('b',b));

			layers.insert(std::make_pair('z',z));

			root = layers.find('y');*/
		}


//template <typename I>
//void Graph<I>::print() const
//{
//    for(const auto & i : layers)
//    {
//
//        std::cout << i.first << " => ";
//        for(auto it = i.second.begin(); it != i.second.end(); ++it)
//        {
//            std::cout << it->vname << ' ' << it->index << ' ' << it->count << '\t';
//        }
//        std::cout << std::endl;
//    }
//}

	}


	template <typename TIndex, typename TFloat>
	class ImplicitGraphQuantile : public Quantile<TIndex, TFloat>
	{
	protected:
		//using Quantile<TIndex, TFloat>::grids;
		using Quantile<TIndex, TFloat>::grid_number;
		using Quantile<TIndex, TFloat>::dx;
		using Quantile<TIndex, TFloat>::lb;
		using Quantile<TIndex, TFloat>::ub;

		using Quantile<TIndex, TFloat>::get_grid_value;

		typedef trie_based::TrieLayer<TIndex> sample_type;
		std::shared_ptr<sample_type> sample;

		size_t count_less_binary(const std::vector<TIndex> &layer, TIndex target) const;
		std::pair<size_t, TFloat> quantile_transform(const std::vector<TIndex> &layer, const std::vector<size_t> &psum, size_t ind, TFloat val01) const;
	public:
		ImplicitGraphQuantile() = default;
		ImplicitGraphQuantile(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		using Quantile<TIndex, TFloat>::set_grid_and_gridn;
		void set_sample(const std::vector<std::vector<TIndex>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights) override;
		void set_sample(std::shared_ptr<sample_type> in_sample);
		void transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const override;
		void transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const override;
	};

	template <typename TIndex, typename TFloat>
	ImplicitGraphQuantile<TIndex, TFloat>::ImplicitGraphQuantile(std::vector<TFloat> in_lb,
	    std::vector<TFloat> in_ub,
	    std::vector<size_t> in_gridn) : Quantile<TIndex, TFloat>(in_lb, in_ub, in_gridn)
	{

	}
	template <typename TIndex, typename TFloat>
	void ImplicitGraphQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TIndex>> &in_sample)
	{
		sample = std::make_shared<sample_type>();
		sample->set_dimension(grid_number.size());
		for(const auto & i : in_sample)
			sample->insert(i);
		sample->sort();
	}
	template <typename TIndex, typename TFloat>
	void ImplicitGraphQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample)
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
			sample->insert(temp);
		}
		sample->sort();
	}

	template <typename TIndex, typename TFloat>
	void ImplicitGraphQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights)
	{
		set_sample(in_sample);
	}

	template <typename TIndex, typename TFloat>
	void ImplicitGraphQuantile<TIndex, TFloat>::set_sample(std::shared_ptr<sample_type> in_sample)
	{
		sample = std::move(in_sample);
		sample->sort();
	}
	template <typename TIndex, typename TFloat>
	void ImplicitGraphQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const
	{
//    // first step
//    std::vector<size_t> psum(sample->root->second.size() + 1, 0);
//    for(size_t j = 1, k = 0; j != sample->root->second.size(); ++j)
//    {
//        k += 1;
//        psum[j] = k;
//    }
//    psum[sample->root->second.size()] = sample->root->second.size();
//
//    auto rez = quantile_transform(sample->root->second, psum, 0, in01[0]);
//    out[0] = rez.second;
//    //p = p->children[rez.first].get();


		auto p = sample->root;
		for(size_t i = 0; i != in01.size(); ++i)
		{
			std::vector<size_t> psum(p->second.size() + 1, 0);
			for(size_t j = 1, k = 0; j != p->second.size(); ++j)
			{
				k += 1;
				psum[j] = k;
			}
			psum[p->second.size()] = p->second.size();

			auto rez = quantile_transform(p->second, psum, i, in01[i]);
			out[i] = rez.second;
			if(i < in01.size() - 1)
				p = sample->layers[i + 1].find(p->second[rez.first]);
		}
	}
	template <typename TIndex, typename TFloat>
	void ImplicitGraphQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const
	{
		auto p = sample->root;
		for(size_t i = 0; i != in01.size(); ++i)
		{
			std::vector<size_t> psum(p->second.size() + 1, 0);
			for(size_t j = 1, k = 0; j != p->second.size(); ++j)
			{
				k += 1;
				psum[j] = k;
			}
			psum[p->second.size()] = p->second.size();

			auto rez = quantile_transform(p->second, psum, i, in01[i]);
			out[i] = rez.first;
			if(i < in01.size() - 1)
				p = sample->layers[i + 1].find(p->second[rez.first]);
		}
	}
	template <typename TIndex, typename TFloat>
	std::pair<size_t, TFloat> ImplicitGraphQuantile<TIndex, TFloat>::quantile_transform(const std::vector<TIndex> &layer, const std::vector<size_t> &psum, size_t ind, TFloat val01) const
	{
		size_t m = 0, count = grid_number[ind], step, c1 = 0, c2 = 0;
		TFloat f1 = 0.0, f2 = 0.0, sample_size_u = static_cast<TFloat>(layer.size());
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

			c1 = psum[count_less_binary(layer, m)];
			f1 = c1/sample_size_u;

			if(f1 < val01)
			{
				c2 = psum[count_less_binary(layer, m + 1)];
				f2 = c2/sample_size_u;

				if(val01 < f2)
					break;

				first = ++it;
				count -= step + 1;
			}
			else
			{
				count = step;
			}
		}

		if(count == 0)
		{
			c2 = psum[count_less_binary(layer, m + 1)];
			f2 = c2/sample_size_u;
		}

		if(c1 == c2)
		{
			if(c1 == 0)
			{
				return std::make_pair(0, get_grid_value(ind, layer.front()) + 2.0*val01*dx[ind]);
			}
			if(c1 == layer.size())
			{
				return std::make_pair(layer.size() - 1, get_grid_value(ind, layer.back()) + 2.0*val01*dx[ind]);
			}

			TIndex target = m;
			auto pos = std::lower_bound(layer.begin(), layer.end(), target);
			size_t index = std::distance(layer.begin(), pos);

			if(index > 0)
			{
				int curr1 = std::abs(static_cast<int>(layer[index]) - static_cast<int>(m));
				int curr2 = std::abs(static_cast<int>(layer[index - 1]) - static_cast<int>(m));

				if(curr1 < curr2)
				{
					return std::make_pair(index, get_grid_value(ind, layer[index]) + 2.0*val01*dx[ind]);
				}
				else if(curr1 == curr2)
				{
					if(layer[index - 1] < layer[index])
						return std::make_pair(index - 1, get_grid_value(ind, layer[index - 1]) + 2.0*val01*dx[ind]);
					else
						return std::make_pair(index, get_grid_value(ind, layer[index]) + 2.0*val01*dx[ind]);
				}
				else
				{
					return std::make_pair(index - 1, get_grid_value(ind, layer[index - 1]) + 2.0*val01*dx[ind]);
				}

			}
			return std::make_pair(index, get_grid_value(ind, layer[index]) + 2.0*val01*dx[ind]);
		}
		TIndex target = m;
		auto pos = std::lower_bound(layer.begin(), layer.end(), target);
		TIndex index = std::distance(layer.begin(), pos);
		return std::make_pair(index, get_grid_value(ind, m) + (val01 - f1) * (get_grid_value(ind, m + 1) - get_grid_value(ind, m)) / (f2 - f1));
	}
	template <typename TIndex, typename TFloat>
	size_t ImplicitGraphQuantile<TIndex, TFloat>::count_less_binary(const std::vector<TIndex> &layer, TIndex target) const
	{
		auto lb = std::lower_bound(layer.begin(), layer.end(), target);
		size_t pos = std::distance(layer.begin(), lb);
		if(lb == layer.end())
			pos = layer.size(); // to psum! which is layer->children.size() + 1
		return pos; // to psum!
	}


// fsa

	template <typename TIndex, typename TFloat>
	class GraphQuantile : public Quantile<TIndex, TFloat>
	{
	protected:
//		using Quantile<TIndex, TFloat>::grids;
		using Quantile<TIndex, TFloat>::grid_number;
		using Quantile<TIndex, TFloat>::dx;

		using Quantile<TIndex, TFloat>::get_grid_value;

		typedef trie_based::Graph<TIndex> sample_type;
		std::shared_ptr<sample_type> sample;

		size_t count_less_binary(const std::vector<trie_based::invect> &layer, TIndex target) const;
		std::pair<size_t, TFloat> quantile_transform(const std::vector<trie_based::invect> &layer, const std::vector<size_t> &psum, size_t ind, TFloat val01) const;
	public:
		GraphQuantile() = default;
		GraphQuantile(std::vector<TFloat> in_lb, std::vector<TFloat> in_ub, std::vector<size_t> in_gridn);
		using Quantile<TIndex, TFloat>::set_grid_and_gridn;
		void set_sample(const std::vector<std::vector<TIndex>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample) override;
		void set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights) override;
		void transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const override;
		void transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const override;
	};

	template <typename TIndex, typename TFloat>
	GraphQuantile<TIndex, TFloat>::GraphQuantile(std::vector<TFloat> in_lb,
	    std::vector<TFloat> in_ub,
	    std::vector<size_t> in_gridn) : Quantile<TIndex, TFloat>(in_lb, in_ub, in_gridn)
	{
		sample = std::make_shared<sample_type>();
	}

	template <typename TIndex, typename TFloat>
	void GraphQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TFloat>& out) const
	{
		auto p = sample->root;
		for(size_t i = 0, k; i != in01.size(); ++i)
		{
			std::vector<size_t> psum(p->second.size() + 1, 0);
			for(size_t j = 1, m = 0; j != p->second.size(); ++j)
			{
				m += p->second[j-1].count;
				psum[j] = m;
				psum.back() += p->second[j-1].count;
			}
			psum.back() += p->second.back().count;
			std::tie(k, out[i]) = quantile_transform(p->second, psum, i, in01[i]);
			p = sample->layers.find(p->second[k].vname);
		}
	}
	template <typename TIndex, typename TFloat>
	void GraphQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TIndex>> &in_sample)
	{

	}
	template <typename TIndex, typename TFloat>
	void GraphQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample)
	{

	}
	template <typename TIndex, typename TFloat>
	void GraphQuantile<TIndex, TFloat>::set_sample(const std::vector<std::vector<TFloat>> &in_sample, const std::vector<size_t> &weights)
	{

	}
	template <typename TIndex, typename TFloat>
	void GraphQuantile<TIndex, TFloat>::transform(const std::vector<TFloat>& in01, std::vector<TIndex>& out) const
	{
		auto p = sample->root;
		for(size_t i = 0; i != in01.size(); ++i)
		{
			std::vector<size_t> psum(p->second.size() + 1, 0);
			for(size_t j = 1, m = 0; j != p->second.size(); ++j)
			{
				m += p->second[j-1].count;
				psum[j] = m;
				psum.back() += p->second[j-1].count;
			}
			psum.back() += p->second.back().count;
			//std::tie(k, out[i]) = quantile_transform(p->second, psum, i, in01[i]);
			auto [k, result] = quantile_transform(p->second, psum, i, in01[i]);
			p = sample->layers.find(p->second[k].vname);
		}
	}
	template <typename TIndex, typename TFloat>
	std::pair<size_t, TFloat> GraphQuantile<TIndex, TFloat>::quantile_transform(const std::vector<trie_based::invect> &layer, const std::vector<size_t> &psum, size_t ind, TFloat val01) const
	{
		size_t m = 0, count = grid_number[ind], step, c1 = 0, c2 = 0;
		TFloat f1 = 0.0, f2 = 0.0, sample_size_u = static_cast<TFloat>(psum.back());
		//auto first = grids[ind].begin();
		//auto it = grids[ind].begin();
		size_t it = 0, first = 0;

		while(count > 0)
		{
			it = first;
			step = count / 2;
			it += step;
			m = it;
//			std::advance(it, step);
//			m = std::distance(grids[ind].begin(), it);

			c1 = psum[count_less_binary(layer, m)];
			f1 = c1/sample_size_u;

			if(f1 < val01)
			{
				c2 = psum[count_less_binary(layer, m + 1)];
				f2 = c2/sample_size_u;

				if(val01 < f2)
					break;

				first = ++it;
				count -= step + 1;
			}
			else
			{
				count = step;
			}
		}

		if(count == 0)
		{
			c2 = psum[count_less_binary(layer, m + 1)];
			f2 = c2/sample_size_u;
		}

		if(c1 == c2)
		{
			if(c1 == 0)
			{
				return std::make_pair(0, get_grid_value(ind, layer.front().index) + 2.0*val01*dx[ind]);
			}
			if(c1 == psum.back())
			{
				return std::make_pair(layer.size() - 1, get_grid_value(ind, layer.back().index) + 2.0*val01*dx[ind]);
			}

			TIndex target = m;
			auto pos = std::lower_bound(layer.begin(), layer.end(), target,
			                            [](const trie_based::invect &l,
			                               const TIndex &r)
			{
				return l.index < r;
			});
			size_t index = std::distance(layer.begin(), pos);

			if(index > 0)
			{
				int curr1 = std::abs(static_cast<int>(layer[index].index) - static_cast<int>(m));
				int curr2 = std::abs(static_cast<int>(layer[index - 1].index) - static_cast<int>(m));

				if(curr1 < curr2)
				{
					return std::make_pair(index, get_grid_value(ind, layer[index].index) + 2.0*val01*dx[ind]);
				}
				else if(curr1 == curr2)
				{
					if(layer[index - 1].index < layer[index].index)
						return std::make_pair(index - 1, get_grid_value(ind, layer[index - 1].index) + 2.0*val01*dx[ind]);
					else
						return std::make_pair(index, get_grid_value(ind, layer[index].index) + 2.0*val01*dx[ind]);
				}
				else
				{
					return std::make_pair(index - 1, get_grid_value(ind, layer[index - 1].index) + 2.0*val01*dx[ind]);
				}

			}
			return std::make_pair(index, get_grid_value(ind, layer[index].index) + 2.0*val01*dx[ind]);
		}
		TIndex target = m;
		auto pos = std::lower_bound(layer.begin(), layer.end(), target,
		                            [](const trie_based::invect &l,
		                               const TIndex &r)
		{
			return l.index < r;
		});
		TIndex index = std::distance(layer.begin(), pos);
		return std::make_pair(index, get_grid_value(ind, m) + (val01 - f1) * (get_grid_value(ind, m + 1) - get_grid_value(ind, m)) / (f2 - f1));
	}
	template <typename TIndex, typename TFloat>
	size_t GraphQuantile<TIndex, TFloat>::count_less_binary(const std::vector<trie_based::invect> &layer, TIndex target) const
	{
		auto lb = std::lower_bound(layer.begin(), layer.end(), target,
		                           [](const trie_based::invect &l,
		                              const TIndex &r)
		{
			return l.index < r;
		});
		size_t pos = std::distance(layer.begin(), lb);
		if(lb == layer.end())
			pos = layer.size(); // to psum! which is layer->children.size() + 1
		return pos; // to psum!
	}

}


#endif
