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
//#include <unordered_map>
#include <algorithm>
#include <memory>

namespace mveqf
{
	namespace trie_based
	{
		template <typename TNode, typename TIndex>
		class TrieBased
		{
		public:
			std::shared_ptr<TNode> root;
			std::vector<std::shared_ptr<TNode>> last_layer;
			explicit TrieBased();
			explicit TrieBased(size_t dim);
			TrieBased(const TrieBased&) = delete;
			TrieBased& operator=(const TrieBased&) = delete;
			~TrieBased();
			void set_dimension(size_t dim);
			size_t get_dimension() const;
			void insert(const std::vector<TIndex> &key);
			void insert(const std::vector<TIndex> &key, size_t number);
			bool search(const std::vector<TIndex> &key) const;
			void fill_tree_count();
			bool empty() const;
			void remove_tree();
			size_t get_total_count() const;
			std::vector<TIndex> get_and_remove_last();
			size_t get_link_count() const;
			size_t get_node_count() const;
			std::map<size_t,std::vector<TIndex>> get_layer_count() const;
		protected:
			void layer_count(size_t current_layer, TNode *p, std::map<size_t,std::vector<TIndex>> &layers) const;
			void get_link_number(TNode *p, size_t &count) const;
			void get_node_number(TNode *p, size_t &count) const;
			void fill_tree_count(TNode *p);
			void get_number(TNode *p, size_t &count) const;
			void is_all_empty(TNode *p) const;
			void delete_last(int dim);
			size_t dimension;
		};
		template <typename TNode, typename TIndex>
		TrieBased<TNode,TIndex>::TrieBased() : root(std::make_shared<TNode>())
		{
			//root = std::make_shared<TNode>();
		}
		template <typename TNode, typename TIndex>
		TrieBased<TNode,TIndex>::TrieBased(size_t dim) : root(std::make_shared<TNode>()), dimension(dim)
		{
//			root = std::make_shared<TNode>();
		}
		template <typename TNode, typename TIndex>
		TrieBased<TNode,TIndex>::~TrieBased() {}
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
			size_t count = 0;
			for(auto &i : last_layer)
			{
				count += i.use_count();
			}
			return count - last_layer.size();
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::get_node_number(TNode *p, size_t &count) const
		{
			for(const auto &i : p->children)
			{
				++count;
				get_node_number(i.get(), count);
			}
		}
		template <typename TNode, typename TIndex>
		size_t TrieBased<TNode,TIndex>::get_node_count() const
		{
			size_t count = 0;
			get_node_number(root.get(), count);
			for(const auto &i : last_layer)
				count -= i.use_count();
			return count + 2*last_layer.size() + 1;
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::get_link_number(TNode *p, size_t &count) const
		{
			for(const auto &i : p->children)
			{
				++count;
				get_node_number(i.get(), count);
			}
		}
		template <typename TNode, typename TIndex>
		size_t TrieBased<TNode,TIndex>::get_link_count() const
		{
			size_t count = 0;
			get_node_number(root.get(), count);
//    for(const auto &i : last_layer)
//        count -= i.use_count();
			return count;
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::insert(const std::vector<TIndex> &key)
		{
			auto p = root.get();
			for(size_t i = 0; i != key.size() - 1; i++)
			{
				auto value = key[i];
				auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
				{
					return obj->index == value;
				});
				if(it == p->children.end())
				{
					p->children.emplace_back(std::make_shared<TNode>(value));
					p->children.shrink_to_fit();
					p = p->children.back().get();
				}
				else
				{
					p = p->children[std::distance(p->children.begin(), it)].get();
				}
			}
			auto value = key.back();
			auto it = std::find_if(last_layer.begin(), last_layer.end(), [&value](const std::shared_ptr<TNode> &obj)
			{
				return obj->index == value;
			});
			size_t dist = 0;
			if(it == last_layer.end())
			{
				last_layer.emplace_back(std::make_shared<TNode>(value));
				last_layer.shrink_to_fit();
				dist = last_layer.size() - 1;
			}
			else
			{
				dist = std::distance(last_layer.begin(), it);
			}

			auto iter = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
			{
				return obj->index == value;
			});
			if(iter == p->children.end())
			{
				std::shared_ptr<TNode> ptr(last_layer[dist]);
				p->children.emplace_back(ptr);
				p->children.shrink_to_fit();
			}
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::insert(const std::vector<TIndex> &key, size_t count)
		{
			auto p = root.get();
			for(size_t i = 0; i != key.size() - 1; i++)
			{
				p->count += count;
				auto value = key[i];
				auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
				{
					return obj->index == value;
				});
				if(it == p->children.end())
				{
					p->children.emplace_back(std::make_shared<TNode>(value));
					p->children.shrink_to_fit();
					p = p->children.back().get();
				}
				else
				{
					p = p->children[std::distance(p->children.begin(), it)].get();
				}
			}
			auto value = key.back();
			auto it = std::find_if(last_layer.begin(), last_layer.end(), [&value](const std::shared_ptr<TNode> &obj)
			{
				return obj->index == value;
			});
			size_t dist = 0;
			if(it == last_layer.end())
			{
				last_layer.emplace_back(std::make_shared<TNode>(value));
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
			auto iter = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
			{
				return obj->index == value;
			});
			if(iter == p->children.end())
			{
				std::shared_ptr<TNode> ptr(last_layer[dist]);
				p->children.emplace_back(ptr);
				p->children.shrink_to_fit();
			}
		}
		template <typename TNode, typename TIndex>
		bool TrieBased<TNode,TIndex>::search(const std::vector<TIndex> &key) const
		{
			auto p = root.get();
			for(size_t i = 0; i != key.size(); i++)
			{
				auto value = key[i];
				auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
				{
					return obj->index == value;
				});
				if(it == p->children.end())
				{
					return false;
				}
				else
				{
					p = p->children[std::distance(p->children.begin(), it)].get();
				}
			}
			return true;
		}
		template <typename TNode, typename TIndex>
		std::vector<TIndex> TrieBased<TNode,TIndex>::get_and_remove_last()
		{
			std::vector<TIndex> sample;

			auto p = root.get();
			if(p->children.empty())
				return sample;

			for(size_t i = 0; i != dimension; ++i)
			{
				p->count--;
				sample.push_back(p->children.back()->index);
				p = p->children.back().get();
			}

			size_t dim = sample.size() - 1;

			delete_last(dim);

			last_layer.erase(
			  std::remove_if(last_layer.begin(), last_layer.end(),
			                 [&sample](const std::shared_ptr<TNode> &obj)
			{
				if(obj->index != sample.back())
					return false;
				if(obj.use_count() == 1)
					return true;
				else
					return false;
			}),
			last_layer.end());

			return sample;
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::remove_tree()
		{
			auto p = root.get();
			while(!p->children.empty())
			{
				auto t = get_and_remove_last();
			}
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::delete_last(int dim)
		{
			if(dim < 0)
				return;

			auto p = root.get();
			if(p->children.empty())
				return;

			for(int i = 0; i != dim; ++i)
			{
				p = p->children.back().get();
			}
			p->children.pop_back();

			if(p->children.empty())
			{
				dim = dim - 1;
				delete_last(dim);
			}
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::fill_tree_count()
		{
			fill_tree_count(root.get());
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
				get_number(i.get(), count);
				i->count = count > 0 ? count : 1;
				fill_tree_count(i.get());
			}
		}
		template <typename TNode, typename TIndex>
		void TrieBased<TNode,TIndex>::get_number(TNode *p, size_t &count) const
		{
			for(auto &i : p->children)
			{
				if(i->children.empty())
					++count;
				get_number(i.get(), count);
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
				layer_count(current_layer + 1, i.get(), layers);
			}
		}
		template <typename TNode, typename TIndex>
		std::map<size_t, std::vector<TIndex>> TrieBased<TNode,TIndex>::get_layer_count() const
		{
			std::map<size_t, std::vector<TIndex>> layers;
			size_t cur_layer = 0;
			layer_count(cur_layer, root.get(), layers);
			layers.erase(std::prev(layers.end()));
			return layers;
		}

		template <typename TNode, typename TIndex>
		class TrieBasedInverse : public TrieBased<TNode, TIndex>
		{
		public:
			using TrieBased<TNode, TIndex>::root;
			using TrieBased<TNode, TIndex>::last_layer;

			TrieBasedInverse();
			TrieBasedInverse(size_t dim);
			void insert(const std::vector<TIndex> &key);
			void insert(const std::vector<TIndex> &key, size_t number);
		};

		template <typename TNode, typename TIndex>
		TrieBasedInverse<TNode,TIndex>::TrieBasedInverse() : TrieBased<TNode, TIndex>() { }

		template <typename TNode, typename TIndex>
		TrieBasedInverse<TNode,TIndex>::TrieBasedInverse(size_t dim) : TrieBased<TNode, TIndex>(dim) {	}

		template <typename TNode, typename TIndex>
		void TrieBasedInverse<TNode,TIndex>::insert(const std::vector<TIndex> &key)
		{
			auto p = root.get();
			auto shrd = root;
			for(size_t i = 0; i != key.size() - 1; i++)
			{
				auto value = key[i];
				auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
				{
					return obj->index == value;
				});
				if(it == p->children.end())
				{
					p->children.emplace_back(std::make_shared<TNode>(value));
					p->children.shrink_to_fit();
					p->children.back()->parent = shrd;
					shrd = p->children.back();
					p = p->children.back().get();
				}
				else
				{
					shrd = p->children[std::distance(p->children.begin(), it)];
					p = p->children[std::distance(p->children.begin(), it)].get();
				}
			}
			auto value = key.back();
			auto it = std::find_if(last_layer.begin(), last_layer.end(), [&value](const std::shared_ptr<TNode> &obj)
			{
				return obj->index == value;
			});
			size_t dist = 0;
			if(it == last_layer.end())
			{
				last_layer.emplace_back(std::make_shared<TNode>(value));
				last_layer.shrink_to_fit();
				dist = last_layer.size() - 1;
			}
			else
			{
				dist = std::distance(last_layer.begin(), it);
			}

			auto iter = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
			{
				return obj->index == value;
			});
			if(iter == p->children.end())
			{
				std::shared_ptr<TNode> ptr(last_layer[dist]);
				p->children.emplace_back(ptr);
				p->children.shrink_to_fit();
				p->children.back()->parent = shrd;
			}
		}

		template <typename TNode, typename TIndex>
		void TrieBasedInverse<TNode,TIndex>::insert(const std::vector<TIndex> &key, size_t count)
		{
			auto p = root.get();
			auto shrd = root;
			for(size_t i = 0; i != key.size() - 1; i++)
			{
				p->count += count;
				auto value = key[i];
				auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
				{
					return obj->index == value;
				});
				if(it == p->children.end())
				{
					p->children.emplace_back(std::make_shared<TNode>(value));
					p->children.shrink_to_fit();
					p->children.back()->parent = shrd;
					shrd = p->children.back();
					p = p->children.back().get();
				}
				else
				{
					shrd = p->children[std::distance(p->children.begin(), it)];
					p = p->children[std::distance(p->children.begin(), it)].get();
				}
			}
			auto value = key.back();
			auto it = std::find_if(last_layer.begin(), last_layer.end(), [&value](const std::shared_ptr<TNode> &obj)
			{
				return obj->index == value;
			});
			size_t dist = 0;
			if(it == last_layer.end())
			{
				last_layer.emplace_back(std::make_shared<TNode>(value));
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
			auto iter = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<TNode> &obj)
			{
				return obj->index == value;
			});
			if(iter == p->children.end())
			{
				std::shared_ptr<TNode> ptr(last_layer[dist]);
				p->children.emplace_back(ptr);
				p->children.shrink_to_fit();
				p->children.back()->parent = shrd;
			}
		}
	}
}
#endif
