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

#ifndef QUANTILE_H
#define QUANTILE_H

#include "trie_based.h"
#include "print2file.h"

namespace empirical_quantile
{
    template <typename T, typename U>
    class ImplicitQuantile
    {
    private:
        typedef trie_based::TrieBased<trie_based::NodeCount<T>,T> sample_type;
        sample_type sample;

        std::vector<U> lb;
        std::vector<U> ub;
        std::vector<std::vector<U>> grids;
        std::vector<size_t> grid_number;

        std::pair<size_t, U> quantile_transform(const std::vector<std::pair<T, T>> &row,
                size_t sample_size, size_t ind, U val01) const;
    public:
        ImplicitQuantile(std::vector<U> in_lb,
                         std::vector<U> in_ub,
                         std::vector<size_t> in_gridn,
                         std::vector<std::vector<T>> in_sample);
        void transform(const std::vector<U>& in01, std::vector<U>& out) const;
    };

    template <typename T, typename U>
    ImplicitQuantile<T, U>::ImplicitQuantile(std::vector<U> in_lb,
            std::vector<U> in_ub,
            std::vector<size_t> in_gridn,
            std::vector<std::vector<T>> in_sample)
    {
        lb = in_lb;
        ub = in_ub;
        grid_number = in_gridn;

        grids.resize(grid_number.size());
        for(size_t i = 0; i != grids.size(); i++)
        {
            std::vector<U> grid(grid_number[i] + 1);
            U startp = lb[i];
            U endp = ub[i];
            U es = endp - startp;
            for(size_t j = 0; j != grid.size(); j++)
            {
                grid[j] = startp + j*es/U(grid_number[i]);
            }
            grids[i] = grid;
            //dx[i] = es/(float(grid_number[i])*2);
        }

        sample.set_dimension(grids.size());
        for(const auto & i : in_sample)
            sample.insert(i);
        sample.fill_tree_count();
    }

    template <typename T, typename U>
    void ImplicitQuantile<T, U>::transform(const std::vector<U>& in01, std::vector<U>& out) const
    {
        auto p = sample.root.get();
        for(size_t i = 0; i != in01.size(); i++)
        {
            std::vector<std::pair<T, T>> row;
            size_t cc = 0;
            for(size_t j = 0; j != p->children.size(); j++)
            {
                row.push_back(std::make_pair(p->children[j]->index, p->children[j]->count));
                cc += p->children[j]->count;
            }

            auto rez2 = quantile_transform(row, cc, i, in01[i]);
            out[i] = rez2.second;

            T index = 0;
            for(size_t j = 1; j < p->children.size(); j++)
            {
                if(p->children[j]->index == T(rez2.first))
                    index = j;
            }
            p = p->children[index].get();
        }
    }
    template <typename T, typename U>
    std::pair<size_t, U> ImplicitQuantile<T, U>::quantile_transform(const std::vector<std::pair<T, T>> &row, size_t sample_size, size_t ind, U val01) const
    {
        size_t m = 0, count = grids[ind].size() - 1, step, index1 = 0, index2 = 0;
        U cdf1, cdf2, sample_size_u = static_cast<U>(sample_size);
        auto first = grids[ind].begin();
        auto it = grids[ind].begin();

        while(count > 0)
        {
            it = first;
            step = count / 2;
            std::advance(it, step);
            m = std::distance(grids[ind].begin(), it);
            
            index1 = 0;
            index2 = 0;
            for(const auto & i:row)
            {
                if(static_cast<size_t>(i.first) < m)
                {
                    index1 += i.second;
                }
                if(static_cast<size_t>(i.first) < m + 1)
                {
                    index2 += i.second;
                }
            }

            cdf1 = index1/sample_size_u;
            cdf2 = index2/sample_size_u;

            if(val01 > cdf1 && val01 < cdf2)
                break;

            if(cdf1 < val01)
            {
                first = ++it;
                count -= step + 1;
            }
            else
                count = step;
        }
        if(index1 == index2)
            return it == grids[ind].begin() ? std::make_pair(size_t(0), grids[ind].front()) : std::make_pair(size_t(grids[ind].size() - 1), grids[ind].back());
//        U x1 = grids[ind][m + 1], y1 = cdf2;
//        return std::make_pair(m, x0 + (val01 - y0) * (x1 - x0) / (y1 - y0));
        return std::make_pair(m, grids[ind][m] + (val01 - cdf1) * (grids[ind][m + 1] - grids[ind][m]) / (cdf2 - cdf1));
    }
}

#endif