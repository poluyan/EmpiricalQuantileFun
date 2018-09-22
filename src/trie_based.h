#ifndef TRIE_H
#define TRIE_H

#include <vector>
#include <algorithm>
#include <memory>

namespace trie_based
{
    struct Node
    {
        int index;
        std::vector<std::shared_ptr<Node>> children;
        size_t count;
        Node() : index(0), count(0) { }
        Node(int ind) : index(ind), count(0) { }
    };

    class TrieBased
    {
    public:
        size_t dimension;
        std::vector<std::shared_ptr<Node>> last_layer;
    public:
        std::shared_ptr<Node> root;
        TrieBased(size_t dim) : dimension(dim)
        {
            root = std::make_shared<Node>();
        };
        ~TrieBased() {};
        void insert(const std::vector<int> &key);
        bool search(const std::vector<int> &key) const;
        void fill_tree_count();
        void fill_tree_count(Node *p);
        void get_number(Node *p, size_t &count) const;
        bool empty() const;
        void remove_tree();
        void is_all_empty(Node *p) const;
        size_t get_total_count() const;
        void delete_last(int dim);
        std::vector<int> get_and_remove_last();
    };
    bool TrieBased::empty() const
    {
        return root->children.empty();
    }
    size_t TrieBased::get_total_count() const
    {
        size_t count = 0;
        for(auto &i : last_layer)
        {
            count += i.use_count();
        }
        return count - last_layer.size();
    }
    void TrieBased::insert(const std::vector<int> &key)
    {
        auto p = root.get();
        for(size_t i = 0; i != key.size() - 1; i++)
        {
            auto value = key[i];
            auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<Node> &obj)
            {
                return obj->index == value;
            });
            if(it == p->children.end())
            {
                p->children.emplace_back(std::make_shared<Node>(value));
                p = p->children.back().get();
            }
            else
            {
                p = p->children[std::distance(p->children.begin(), it)].get();
            }
        }
        auto value = key.back();
        auto it = std::find_if(last_layer.begin(), last_layer.end(), [&value](const std::shared_ptr<Node> &obj)
        {
            return obj->index == value;
        });
        size_t dist = 0;
        if(it == last_layer.end())
        {
            last_layer.push_back(std::make_shared<Node>(value));
            dist = last_layer.size() - 1;
        }
        else
        {
            dist = std::distance(last_layer.begin(), it);
        }

        it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<Node> &obj)
        {
            return obj->index == value;
        });
        if(it == p->children.end())
        {
            std::shared_ptr<trie_based::Node> ptr(last_layer[dist]);
            p->children.push_back(ptr);
        }
    }
    bool TrieBased::search(const std::vector<int> &key) const
    {
        auto p = root.get();
        for(size_t i = 0; i != key.size(); i++)
        {
            auto value = key[i];
            auto it = std::find_if(p->children.begin(), p->children.end(), [&value](const std::shared_ptr<Node> &obj)
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

    std::vector<int> TrieBased::get_and_remove_last()
    {
        std::vector<int> sample;

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
                           [&sample](const std::shared_ptr<Node> &obj)
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

    void TrieBased::remove_tree()
    {
        auto p = root.get();
        while(!p->children.empty())
        {
            auto t = get_and_remove_last();
        }
    }

    void TrieBased::delete_last(int dim)
    {
        if(dim < 0)
            return;

        auto p = root.get();
        if(p->children.empty())
            return;

        for(size_t i = 0; i != dim; ++i)
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
    void TrieBased::fill_tree_count()
    {
        fill_tree_count(root.get());
    }
    void TrieBased::fill_tree_count(Node *p)
    {
        for(auto &i : p->children)
        {
            size_t count = 0;
            get_number(i.get(), count);
            i->count = count > 0 ? count : 1;
            fill_tree_count(i.get());
        }
    }
    void TrieBased::get_number(Node *p, size_t &count) const
    {
        for(auto &i : p->children)
        {
            if(i->children.empty())
                ++count;
            get_number(i.get(), count);
        }
    }
}

#endif
