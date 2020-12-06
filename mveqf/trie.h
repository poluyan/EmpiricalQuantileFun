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
#ifndef TRIE_H
#define TRIE_H

#include <memory>
#include <vector>
#include <algorithm>
#include <mveqf/trie_node.h>

namespace mveqf
{
	template <typename TNode, typename TIndex>
	struct BaseSample
	{
		virtual void set_dimension(size_t dim) = 0;
		virtual size_t get_dimension() const = 0;
		virtual void insert(const std::vector<TIndex> &key) = 0;
		virtual void insert(const std::vector<TIndex> &key, size_t number) = 0;
		virtual bool search(const std::vector<TIndex> &key) const = 0;
		virtual void fill_tree_count() = 0;
		virtual size_t get_link_count() const = 0;
		virtual size_t get_node_count() const = 0;
		virtual ~BaseSample() = default;
	};

	template <typename TNode, typename TIndex>
	class Trie : public BaseSample<TNode, TIndex>
	{
	protected:
		size_t dimension;
	public:
		TNode *root;
		Trie();
		Trie(size_t dim);
		void destroy(TNode *p);
		~Trie();
		void set_dimension(size_t dim) override;
		size_t get_dimension() const override;
		bool empty() const;
		void insert(const std::vector<TIndex> &key, size_t number) override;
		void insert(const std::vector<TIndex> &key) override;
		bool search(const std::vector<TIndex> &key) const override;
		void fill_tree_count() override;

		size_t get_link_count() const override;
		size_t get_node_count() const override;
	};
	template <typename TNode, typename TIndex>
	Trie<TNode,TIndex>::Trie() : dimension(0), root(new TNode())
	{
	}
	template <typename TNode, typename TIndex>
	Trie<TNode,TIndex>::Trie(size_t dim) : dimension(dim), root(new TNode())
	{
	}

	template <typename TNode, typename TIndex>
	void Trie<TNode,TIndex>::destroy(TNode *p)
	{
		for(auto &i : p->children)
			destroy(i);
		delete p;
	}

	template <typename TNode, typename TIndex>
	Trie<TNode,TIndex>::~Trie()
	{
		for(TNode *i : root->children)
			destroy(i);
		delete root;
	}
	template <typename TNode, typename TIndex>
	void Trie<TNode,TIndex>::set_dimension(size_t dim)
	{
		dimension = dim;
	}
	template <typename TNode, typename TIndex>
	size_t Trie<TNode,TIndex>::get_dimension() const
	{
		return dimension;
	}
	template <typename TNode, typename TIndex>
	bool Trie<TNode,TIndex>::empty() const
	{
		return root->children.empty();
	}
	template <typename TNode, typename TIndex>
	void Trie<TNode,TIndex>::insert(const std::vector<TIndex> &key, size_t count)
	{
		auto p = root;
		for(const auto &i : key)
		{
			p->count += count;
			auto it = std::find_if(p->children.begin(), p->children.end(), [&i](const auto &obj)
			{
				return obj->index == i;
			});
			if(it == p->children.end())
			{
				p->children.push_back(new TNode(i));
				p->children.shrink_to_fit();
				p = p->children.back();
			}
			else
			{
				p = p->children[std::distance(p->children.begin(), it)];
			}
		}
		p->count += count;
	}
	template <typename TNode, typename TIndex>
	void Trie<TNode,TIndex>::insert(const std::vector<TIndex> &key)
	{
		insert(key, 1);
	}
	template <typename TNode, typename TIndex>
	bool Trie<TNode,TIndex>::search(const std::vector<TIndex> &key) const
	{
		auto p = root;
		for(const auto &i : key)
		{
			auto it = std::find_if(p->children.begin(), p->children.end(), [&i](const auto &obj)
			{
				return obj->index == i;
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
	void Trie<TNode,TIndex>::fill_tree_count()
	{
	}
	template <typename TNode, typename TIndex>
	size_t Trie<TNode,TIndex>::get_link_count() const
	{
		return 0;
	}
	template <typename TNode, typename TIndex>
	size_t Trie<TNode,TIndex>::get_node_count() const
	{
		return 0;
	}
}

#endif
