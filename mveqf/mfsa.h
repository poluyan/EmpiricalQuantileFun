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
#ifndef MFSA_H
#define MFSA_H

#include <memory>
#include <set>
#include <vector>
#include <algorithm>
#include <mveqf/cstvect.h>

namespace mveqf
{
	namespace mfsa
	{
		template <typename TIndex>
		struct Node
		{
			bool accepting_state;
			size_t count;
			size_t in_count;
			cst::vector<std::pair<TIndex, std::shared_ptr<Node<TIndex>>>> children;

			Node<TIndex>(bool accept);
			Node<TIndex>(const Node<TIndex> *p);
			Node<TIndex> *transition(TIndex label) const;
			std::shared_ptr<Node<TIndex>> transition_shared(TIndex label) const;
			Node<TIndex> *transition(const std::vector<TIndex> &str);
			Node<TIndex> *add_node(TIndex label);
			bool same_path(Node<TIndex> *node) const;
			bool is_equal(Node<TIndex> *obj) const;
		};

		template <typename TIndex>
		Node<TIndex>::Node(bool accept): accepting_state(accept), count(0), in_count(0) {}

		template <typename TIndex>
		Node<TIndex>::Node(const Node<TIndex> *p): accepting_state(p->accepting_state), count(0), in_count(0), children(p->children)
		{
			for(const auto &i : children)
				i.second->in_count++;
		}

		template <typename TIndex>
		Node<TIndex> *Node<TIndex>::transition(TIndex label) const
		{
			auto it = std::find_if(children.begin(), children.end(),
			                       [&label](const std::pair<TIndex, std::shared_ptr<Node<TIndex>>> &obj)
			{
				return obj.first == label;
			});
			return it != children.end() ? children[std::distance(children.begin(), it)].second.get() : nullptr;
		}

		template <typename TIndex>
		std::shared_ptr<Node<TIndex>> Node<TIndex>::transition_shared(TIndex label) const
		{
			auto it = std::find_if(children.begin(), children.end(),
			                       [&label](const std::pair<TIndex, std::shared_ptr<Node<TIndex>>> &obj)
			{
				return obj.first == label;
			});
			return it != children.end() ? children[std::distance(children.begin(), it)].second : nullptr;
		}

		template <typename TIndex>
		Node<TIndex> *Node<TIndex>::transition(const std::vector<TIndex> &key)
		{
			Node<TIndex> *current = this;
			for(const auto &i : key)
			{
				current = current->transition(i);
				if(current == nullptr)
				{
					break;
				}
			}
			return current;
		}

		template <typename TIndex>
		Node<TIndex> *Node<TIndex>::add_node(TIndex label)
		{
			std::shared_ptr<Node<TIndex>> p = std::make_shared<Node<TIndex>>(false);
			p->in_count++;

			auto it = std::find_if(children.begin(), children.end(),
			                       [&label](const std::pair<TIndex, std::shared_ptr<Node<TIndex>>> &obj)
			{
				return obj.first == label;
			});
			if(it != children.end())
				children[std::distance(children.begin(), it)].second = p;
			else
			{
				children.emplace_back(std::make_pair(label, p));
				children.shrink_to_fit();
			}
			return p.get();
		}

		template <typename TIndex>
		bool Node<TIndex>::same_path(Node<TIndex> *p) const
		{
			if(this->children.size() == p->children.size())
			{
				for(const auto &i : this->children)
				{
					TIndex label = i.first;
					if(p->children.size() == 1 && p->children.front().first != label)
						return false;
					auto it = std::find_if(p->children.begin(), p->children.end(),
					                       [&label](const std::pair<TIndex, std::shared_ptr<Node<TIndex>>> &obj)
					{
						return obj.first == label;
					});
					if(it == p->children.end())
						return false;
					if(!it->second->is_equal(i.second.get()))
					{
						return false;
					}
				}
			}
			else
			{
				return false;
			}
			return true;
		}

		template <typename TIndex>
		bool Node<TIndex>::is_equal(Node<TIndex> *obj) const
		{
			bool equal = this == obj;
			if(!equal && obj != nullptr)
				equal = accepting_state == obj->accepting_state && same_path(obj);
			return equal;
		}

		template <typename TIndex>
		class MFSA
		{
		public:
			std::shared_ptr<Node<TIndex>> root;
			std::shared_ptr<Node<TIndex>> final_state;

			MFSA();
			void set_dimension(size_t dim);
			void insert(const std::vector<TIndex> &key);
			bool contains(const std::vector<TIndex> &key) const;

			void fill_tree_count(Node<TIndex> *p);
			void fill_tree_count();

			void get_number(Node<TIndex> *p, size_t &count) const;
			void get_link(Node<TIndex>* p, size_t &count) const;
			size_t count_nodes(Node<TIndex> *current, std::set<std::shared_ptr<Node<TIndex>>> &data) const;
			std::pair<size_t, size_t> get_node_link_count() const;

			size_t get_node_count() const;
			size_t get_link_count() const;
		protected:
			size_t dimension;
			std::vector<std::pair<std::shared_ptr<Node<TIndex>>, std::shared_ptr<Node<TIndex>>>> eq;
			
			std::vector<TIndex> longest_prefix(const std::vector<TIndex> &key) const;
			void replace_or_register(const std::shared_ptr<Node<TIndex>> &p, const std::vector<TIndex> &key);
			void add_path(Node<TIndex> *p, const std::vector<TIndex> &key);
			void remove_path(const std::vector<TIndex> &key);
			void clone_path(Node<TIndex> *pivot, const std::vector<TIndex> &to_pivot, const std::vector<TIndex> &key) const;
			void add(const std::vector<TIndex> &key);
		};

		template <typename TIndex>
		void MFSA<TIndex>::set_dimension(size_t dim)
		{
			dimension = dim;
		}

		template <typename TIndex>
		void MFSA<TIndex>::fill_tree_count()
		{
			fill_tree_count(root.get());
			size_t count = 0;
			for(const auto &i : root->children)
				count += i.second->count;
			root->count = count;
		}

		template <typename TIndex>
		void MFSA<TIndex>::fill_tree_count(Node<TIndex> *p)
		{
			for(const auto &i : p->children)
			{
				size_t count = 0;
				get_number(i.second.get(), count);
				i.second->count = count > 0 ? count : 1;
				fill_tree_count(i.second.get());
			}
		}

		template <typename TIndex>
		void MFSA<TIndex>::get_number(Node<TIndex> *p, size_t &count) const
		{
			for(const auto &i : p->children)
			{
				if(i.second->children.empty())
					++count;
				get_number(i.second.get(), count);
			}
		}

		template <typename TIndex>
		MFSA<TIndex>::MFSA() : root(std::make_shared<Node<TIndex>>(false)), final_state(std::make_shared<Node<TIndex>>(true)), dimension(0) {}

		template <typename TIndex>
		void MFSA<TIndex>::insert(const std::vector<TIndex> &key)
		{
			add(key);
			replace_or_register(root, key);
		}

		template <typename TIndex>
		std::vector<TIndex> MFSA<TIndex>::longest_prefix(const std::vector<TIndex> &key) const
		{
			Node<TIndex> *current = root.get();
			size_t prefix = 0;
			for(const auto &i : key)
			{
				auto it = std::find_if(current->children.begin(), current->children.end(),
				                       [&i](const std::pair<TIndex, std::shared_ptr<Node<TIndex>>> &obj)
				{
					return obj.first == i;
				});
				if(it == current->children.end())
					break;
				current = current->transition(i);
				++prefix;
			}
			return std::vector<TIndex>(key.begin(), key.begin() + prefix);
		}

		template <typename TIndex>
		void MFSA<TIndex>::replace_or_register(const std::shared_ptr<Node<TIndex>> &p, const std::vector<TIndex> &key)
		{
			TIndex label = key.front();
			auto t = p->transition_shared(label);
			std::vector<TIndex> subvect(key.begin() + 1, key.begin() + key.size());
			if(t->children.size() > 0 && subvect.size() > 0)
				replace_or_register(t, subvect);

			std::shared_ptr<Node<TIndex>> en = nullptr;
			for(size_t j = 0; j != eq.size(); j++)
			{
				if(eq[j].first->is_equal(t.get()))
				{
					en = eq[j].second;
					break;
				}
			}
			if(en == nullptr)
			{
				eq.push_back(std::make_pair(t, t));
			}
			else if(en != t)
			{
				for(const auto &i : t->children)
					i.second->in_count--;

				t->in_count--;
				en->in_count++;
				auto it = std::find_if(p->children.begin(), p->children.end(),
				                       [&label](const std::pair<TIndex, std::shared_ptr<Node<TIndex>>> &obj)
				{
					return obj.first == label;
				});
				if(it != p->children.end())
					p->children[std::distance(p->children.begin(), it)].second = en;
				else
				{
					p->children.emplace_back(std::make_pair(label, en));
					p->children.shrink_to_fit();
				}
			}
		}

		template <typename TIndex>
		void MFSA<TIndex>::add_path(Node<TIndex> *p, const std::vector<TIndex> &key)
		{
			if(key.size() > 0)
			{
				Node<TIndex> *current = p;
				for(size_t i = 0; i < key.size() - 1; i++)
				{
					current = current->add_node(key[i]);
				}

				std::shared_ptr<Node<TIndex>> p = final_state;
				p->in_count++;
				TIndex label = key.back();
				auto it = std::find_if(current->children.begin(), current->children.end(),
				                       [&label](const std::pair<TIndex, std::shared_ptr<Node<TIndex>>> &obj)
				{
					return obj.first == label;
				});
				if(it != current->children.end())
					current->children[std::distance(current->children.begin(), it)].second = p;
				else
				{
					current->children.emplace_back(std::make_pair(label, p));
					current->children.shrink_to_fit();
				}
			}
			else
			{
				p->accepting_state = true;
			}
		}

		template <typename TIndex>
		void MFSA<TIndex>::remove_path(const std::vector<TIndex> &key)
		{
			Node<TIndex> *current = root.get();
			for(const auto &i : key)
			{
				current = current->transition(i);
				for(size_t j = 0; j != eq.size(); j++)
				{
					if(eq[j].first->is_equal(current) && eq[j].second->is_equal(current))
					{
						eq.erase(eq.begin() + j);
						break;
					}
				}
			}
		}

		template <typename TIndex>
		void MFSA<TIndex>::clone_path(Node<TIndex> *pivot, const std::vector<TIndex> &to_pivot, const std::vector<TIndex> &key) const
		{
			Node<TIndex> *last_target = pivot->transition(key);
			std::shared_ptr<Node<TIndex>> last = nullptr;
			for(int i = key.size(); i >= 0; i--)
			{
				std::vector<TIndex> current_key;
				if(i > 0)
					std::copy_n(key.begin(), i, std::back_inserter(current_key));

				Node<TIndex> *target = i > 0 ? pivot->transition(current_key) : pivot;
				std::shared_ptr<Node<TIndex>> cloned;

				if(i == 0)
				{
					std::vector<TIndex> key_to_pivot(to_pivot.begin(), to_pivot.begin() + to_pivot.size() - 1);
					TIndex parent_label = to_pivot.back();

					Node<TIndex> *parent = root->transition(key_to_pivot);
					cloned = std::make_shared<Node<TIndex>>(pivot);
					pivot->in_count--;
					cloned->in_count++;
					auto it = std::find_if(parent->children.begin(), parent->children.end(),
					                       [&parent_label](const std::pair<TIndex, std::shared_ptr<Node<TIndex>>> &obj)
					{
						return obj.first == parent_label;
					});
					if(it != parent->children.end())
						parent->children[std::distance(parent->children.begin(), it)].second = cloned;
					else
					{
						parent->children.emplace_back(std::make_pair(parent_label, cloned));
						parent->children.shrink_to_fit();
					}
				}
				else
				{
					cloned = std::make_shared<Node<TIndex>>(target);
				}

				if(last != nullptr)
				{
					last_target->in_count--;
					last->in_count++;
					TIndex last_label = key[i]; // last = nullptr at first iteration
					auto it = std::find_if(cloned->children.begin(), cloned->children.end(),
					                       [&last_label](const std::pair<TIndex, std::shared_ptr<Node<TIndex>>> &obj)
					{
						return obj.first == last_label;
					});

					if(it != cloned->children.end())
						cloned->children[std::distance(cloned->children.begin(), it)].second = last;
					else
					{
						cloned->children.emplace_back(std::make_pair(last_label, last));
						cloned->children.shrink_to_fit();
					}
					last_target = target;
				}
				last = cloned;
			}
		}

		template <typename TIndex>
		void MFSA<TIndex>::add(const std::vector<TIndex> &key)
		{
			std::vector<TIndex> prefix = longest_prefix(key);
			std::vector<TIndex> suffix(key.begin() + prefix.size(), key.begin() + key.size());

			size_t i = 0;
			Node<TIndex> *current = root.get();
			for(; i < prefix.size(); i++)
			{
				auto t = prefix[i];
				auto it = std::find_if(current->children.begin(), current->children.end(),
				                       [&t](const std::pair<TIndex, std::shared_ptr<Node<TIndex>>> &obj)
				{
					return obj.first == t;
				});
				current = it != current->children.end() ? current->transition(prefix[i]) : nullptr;
				if(current == nullptr || current->in_count > 1)
				{
					break;
				}
			}

			if(current == root.get() || i == prefix.size())
			{
				remove_path(prefix);
			}
			else
			{
				std::vector<TIndex> subvect(prefix.begin(), prefix.begin() + i);
				remove_path(subvect);
				if(current != nullptr)
				{
					std::vector<TIndex> key_to_first(prefix.begin(), prefix.begin() + i + 1);
					std::vector<TIndex> duplicate(prefix.begin() + i + 1, prefix.begin() + prefix.size());
					clone_path(current, key_to_first, duplicate);
				}
			}
			add_path(root->transition(prefix), suffix);
		}

		template <typename TIndex>
		bool MFSA<TIndex>::contains(const std::vector<TIndex> &key) const
		{
			Node<TIndex> *p = root->transition(key);
			return p != nullptr && p->accepting_state;
		}

		template <typename TIndex>
		size_t MFSA<TIndex>::count_nodes(Node<TIndex> *current, std::set<std::shared_ptr<Node<TIndex>>> &data) const
		{
			for(const auto &i : current->children)
			{
				if(data.find(i.second) == data.end())
					data.insert(i.second);
				count_nodes(i.second.get(), data);
			}
			return data.size();
		}

		template <typename TIndex>
		std::pair<size_t, size_t> MFSA<TIndex>::get_node_link_count() const
		{
			std::set<std::shared_ptr<Node<TIndex>>> temp;
			size_t t = count_nodes(root.get(), temp);
			size_t c = root->children.size();
			for(const auto &i : temp)
				c += i->children.size();
			return std::pair<size_t, size_t>(t + 1, c);
		}

		template <typename TIndex>
		size_t MFSA<TIndex>::get_node_count() const
		{
			std::set<std::shared_ptr<Node<TIndex>>> temp;
			size_t t = count_nodes(root.get(), temp);
			return t + 1;
		}

		template <typename TIndex>
		size_t MFSA<TIndex>::get_link_count() const
		{
			std::set<std::shared_ptr<Node<TIndex>>> temp;
			count_nodes(root.get(), temp);
			size_t c = root->children.size();
			for(const auto &i : temp)
				c += i->children.size();
			return c;
		}

		template <typename TIndex>
		void MFSA<TIndex>::get_link(Node<TIndex>* p, size_t &count) const
		{
			count += p->children.size();
			for(const auto &i: p->children)
				get_link(i.second.get(), count);
		}
	}
}

#endif
