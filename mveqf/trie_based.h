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
#ifndef TRIE_BASED_H
#define TRIE_BASED_H

#include <vector>
#include <map>
#include <algorithm>
#include <memory>
#include <mveqf/cstvect.h>
#include <mveqf/sample.h>

namespace mveqf
{
	template <typename TNode, typename TIndex>
	class TrieBased : public BaseSample<TIndex>
	{
	public:
		TNode *root;
		cst::vector<TNode*> last_layer;
		explicit TrieBased();
		explicit TrieBased(size_t dim);
		TrieBased(const TrieBased&) = delete;
		TrieBased& operator=(const TrieBased&) = delete;
		~TrieBased();
		void set_dimension(size_t dim) override;
		size_t get_dimension() const override;
		void insert(const std::vector<TIndex> &key) override;
		void insert(const std::vector<TIndex> &key, size_t number) override;
		bool search(const std::vector<TIndex> &key) const override;
		void fill_tree_count() override;
		bool empty() const;
		void remove_tree();
		size_t get_total_count() const;
		std::vector<TIndex> get_and_remove_last();
		size_t get_link_count() const override;
		size_t get_node_count() const override;
		std::map<size_t,std::vector<TIndex>> get_layer_count() const;
	protected:
		void layer_count(size_t current_layer, TNode *p, std::map<size_t,std::vector<TIndex>> &layers) const;
		void get_link_number(TNode *p, size_t &count) const;
		void get_node_number(TNode *p, size_t &count, size_t dim) const;
		void fill_tree_count(TNode *p);
		void get_number(TNode *p, size_t &count) const;
		void is_all_empty(TNode *p) const;
		void destroy(TNode *p);
		size_t dimension;
	};

	template <typename TNode, typename TIndex>
	void TrieBased<TNode,TIndex>::destroy(TNode *p)
	{
		if(!p->children.empty())
		{
			for(auto &i : p->children)
				destroy(i);
			delete p;
		}
	}

	template <typename TNode, typename TIndex>
	TrieBased<TNode,TIndex>::TrieBased() : root(new TNode()), dimension(0)
	{
	}
	template <typename TNode, typename TIndex>
	TrieBased<TNode,TIndex>::TrieBased(size_t dim) : root(new TNode()), dimension(dim)
	{
	}
	template <typename TNode, typename TIndex>
	TrieBased<TNode,TIndex>::~TrieBased()
	{
		for(TNode *i : root->children)
			destroy(i);
		for(TNode *i : last_layer)
			delete i;
		delete root;
	}

	template <typename TNode, typename TIndex>
	void TrieBased<TNode,TIndex>::set_dimension(size_t dim)
	{
		dimension = dim;
	}
	template <typename TNode, typename TIndex>
	size_t TrieBased<TNode,TIndex>::get_dimension() const
	{
		return dimension;
	}
	template <typename TNode, typename TIndex>
	bool TrieBased<TNode,TIndex>::empty() const
	{
		return root->children.empty();
	}
	template <typename TNode, typename TIndex>
	size_t TrieBased<TNode,TIndex>::get_total_count() const
	{
		return root->count;
	}
	template <typename TNode, typename TIndex>
	void TrieBased<TNode,TIndex>::get_node_number(TNode *p, size_t &count, size_t dim) const
	{
		if(dim > 1)
		{
			for(TNode *i : p->children)
			{
				++count;
				get_node_number(i, count, dim - 1);
			}
		}
	}
	template <typename TNode, typename TIndex>
	size_t TrieBased<TNode,TIndex>::get_node_count() const
	{
		size_t count = 0;
		get_node_number(root, count, dimension);
		return count + last_layer.size() + 1;
	}
	template <typename TNode, typename TIndex>
	void TrieBased<TNode,TIndex>::get_link_number(TNode *p, size_t &count) const
	{
		for(const auto &i : p->children)
		{
			++count;
			get_link_number(i, count);
		}
	}
	template <typename TNode, typename TIndex>
	size_t TrieBased<TNode,TIndex>::get_link_count() const
	{
		size_t count = 0;
		get_link_number(root, count);
		return count;
	}
	template <typename TNode, typename TIndex>
	void TrieBased<TNode,TIndex>::insert(const std::vector<TIndex> &key)
	{
		auto p = root;
		for(size_t i = 0; i != key.size() - 1; i++)
		{
			auto value = key[i];
			auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const auto &obj)
			{
				return obj->index == value;
			});
			if(it == p->children.end())
			{
				p->children.emplace_back(new TNode(value));
				p->children.shrink_to_fit();
				p = p->children.back();
			}
			else
			{
				p = p->children[std::distance(p->children.begin(), it)];
			}
		}
		auto value = key.back();
		auto it = std::find_if(last_layer.begin(), last_layer.end(), [&value](const auto &obj)
		{
			return obj->index == value;
		});
		size_t dist = 0;
		if(it == last_layer.end())
		{
			last_layer.emplace_back(new TNode(value));
			last_layer.shrink_to_fit();
			dist = last_layer.size() - 1;
		}
		else
		{
			dist = std::distance(last_layer.begin(), it);
		}

		auto iter = std::find_if(p->children.begin(), p->children.end(), [&value](const auto &obj)
		{
			return obj->index == value;
		});
		if(iter == p->children.end())
		{
//				auto t = last_layer[dist];
//				p->children.emplace_back(t);
			p->children.push_back(last_layer[dist]);
//				p->children.emplace_back(new TNode());
//				p->children.back() = last_layer[dist];
			p->children.shrink_to_fit();
		}
	}
	template <typename TNode, typename TIndex>
	void TrieBased<TNode,TIndex>::insert(const std::vector<TIndex> &key, size_t count)
	{
		auto p = root;
		for(size_t i = 0; i != key.size() - 1; i++)
		{
			p->count += count;
			auto value = key[i];
			auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const auto &obj)
			{
				return obj->index == value;
			});
			if(it == p->children.end())
			{
				auto t = new TNode(value);
				p->children.emplace_back(t);
				p->children.shrink_to_fit();
				p = p->children.back();
			}
			else
			{
				p = p->children[std::distance(p->children.begin(), it)];
			}
		}
		auto value = key.back();
		auto it = std::find_if(last_layer.begin(), last_layer.end(), [&value](const auto &obj)
		{
			return obj->index == value;
		});
		size_t dist = 0;
		if(it == last_layer.end())
		{
			auto t = new TNode(value);
			last_layer.emplace_back(t);
//        last_layer.back()->count += count;
			last_layer.back()->count = 1;
			last_layer.shrink_to_fit();
			dist = last_layer.size() - 1;
		}
		else
		{
			dist = std::distance(last_layer.begin(), it);
			last_layer[dist]->count = 1;
//        last_layer[dist]->count += count;
		}
		p->count += count;
		auto iter = std::find_if(p->children.begin(), p->children.end(), [&value](const auto &obj)
		{
			return obj->index == value;
		});
		if(iter == p->children.end())
		{
			p->children.emplace_back(last_layer[dist]);
			p->children.shrink_to_fit();
		}
	}
	template <typename TNode, typename TIndex>
	bool TrieBased<TNode,TIndex>::search(const std::vector<TIndex> &key) const
	{
		auto p = root;
		for(size_t i = 0; i != key.size(); i++)
		{
			auto value = key[i];
			auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const auto &obj)
			{
				return obj->index == value;
			});
			if(it == p->children.end())
			{
				return false;
			}
			else
			{
				p = p->children[std::distance(p->children.begin(), it)];
			}
		}
		return true;
	}
	template <typename TNode, typename TIndex>
	std::vector<TIndex> TrieBased<TNode,TIndex>::get_and_remove_last()
	{
		std::vector<TIndex> sample;
		std::vector<TIndex> back_size;

		auto p = root, rt = root;
		if(root == nullptr || p->children.empty())
			return sample;

		for(size_t i = 0; i != dimension; ++i)
		{
			p->count--;
			sample.push_back(p->children.back()->index);
			back_size.push_back(p->children.size());
			rt = p;
			p = p->children.back();
		}

		bool flag = true;
		for(const auto & i : back_size)
		{
			if(i > 1)
			{
				flag = false;
				break;
			}
		}
		if(flag)
		{
			auto t = root->children.back(), kt = root->children.back();
			root->children.pop_back();
			for(size_t k = 1; k < dimension; k++)
			{
				kt = t;
				t = t->children.back();
				delete kt;
			}
			return sample;
		}

		if(back_size.back() > 1)
		{
			rt->children.pop_back();
		}
		else
		{
			size_t j = 0;
			for(j = back_size.size() - 1; back_size[j] == 1 && j > 0; j--);
			p = root;
			for(size_t k = 0; k <= j; k++)
			{
				rt = p;
				p = p->children.back();
			}
			auto init = rt;
			rt = p;
			for(; j < dimension - 1; j++)
			{
				rt = p;
				p = p->children.back();
				delete rt;
			}
			init->children.pop_back();
		}
		return sample;
	}

	template <typename TNode, typename TIndex>
	void TrieBased<TNode,TIndex>::remove_tree()
	{
		auto p = root;
		while(!p->children.empty())
		{
			auto t = get_and_remove_last();
		}
	}
	template <typename TNode, typename TIndex>
	void TrieBased<TNode,TIndex>::fill_tree_count()
	{
		fill_tree_count(root);
		size_t count = 0;
		for(auto &i : root->children)
		{
			count += i->count;
		}
		root->count = count;
	}
	template <typename TNode, typename TIndex>
	void TrieBased<TNode,TIndex>::fill_tree_count(TNode *p)
	{
		for(auto &i : p->children)
		{
			size_t count = 0;
			get_number(i, count);
			i->count = count > 0 ? count : 1;
			fill_tree_count(i);
		}
	}
	template <typename TNode, typename TIndex>
	void TrieBased<TNode,TIndex>::get_number(TNode *p, size_t &count) const
	{
		for(auto &i : p->children)
		{
			if(i->children.empty())
				++count;
			get_number(i, count);
		}
	}
	template <typename TNode, typename TIndex>
	void TrieBased<TNode,TIndex>::layer_count(size_t current_layer, TNode *p, std::map<size_t, std::vector<TIndex>> &layers) const
	{
		auto it = layers.find(current_layer);
		if(it != layers.end())
		{
			it->second.push_back(p->children.size());
		}
		else
		{
			layers.insert(std::pair<size_t, std::vector<TIndex>>(current_layer, std::vector<TIndex>(1, p->children.size())));
		}
		for(const auto &i : p->children)
		{
			layer_count(current_layer + 1, i, layers);
		}
	}
	template <typename TNode, typename TIndex>
	std::map<size_t, std::vector<TIndex>> TrieBased<TNode,TIndex>::get_layer_count() const
	{
		std::map<size_t, std::vector<TIndex>> layers;
		size_t cur_layer = 0;
		layer_count(cur_layer, root, layers);
		layers.erase(std::prev(layers.end()));
		return layers;
	}
}
#endif
