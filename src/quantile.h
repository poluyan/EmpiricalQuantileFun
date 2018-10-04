#ifndef QUANTILE_H
#define QUANTILE_H

#include "trie_based.h"
#include "print2file.h"

template <typename T>
class Quantile
{
private:
    typedef trie_based::TrieBased<trie_based::NodeCount<T>,T> sample_type;
    sample_type sample;
    std::vector<float> lb;
    std::vector<float> ub;
    std::vector<std::vector<float>> grids;
    std::vector<size_t> grid_number;
public:
    Quantile(std::vector<float> in_lb,
             std::vector<float> in_ub,
             std::vector<size_t> in_gridn,
             std::vector<std::vector<int>> in_sample)
    {
        lb = in_lb;
        ub = in_ub;
        grid_number = in_gridn;

        grids.resize(grid_number.size());
        for(size_t i = 0; i != grids.size(); i++)
        {
            std::vector<float> grid(grid_number[i] + 1);
            float startp = lb[i];
            float endp = ub[i];
            float es = endp - startp;
            for(size_t j = 0; j != grid.size(); j++)
            {
                grid[j] = startp + j*es/float(grid_number[i]);
            }
            grids[i] = grid;
            //dx[i] = es/(float(grid_number[i])*2);
        }

        sample.set_dimension(grids.size());
        for(const auto & i : in_sample)
            sample.insert(i);
        sample.fill_tree_count();
    }
    void tranform(const std::vector<float>& in01, std::vector<float>& out)
    {
        auto p = sample.root.get();
        for(size_t i = 0; i != in01.size(); i++)
        {
            std::vector<std::pair<int,int>> row;
            int cc = 0;
            for(size_t j = 0; j != p->children.size(); j++)
            {
                row.push_back(std::make_pair(p->children[j]->index, p->children[j]->count));
                cc += p->children[j]->count;
            }

            auto rez2 = ecdf1d_pair_fromgrid_trie(row, cc, i, in01[i]);
            out[i] = rez2.second;

            int index = 0;
            for(size_t j = 1; j < p->children.size(); j++)
            {
                if(p->children[j]->index == rez2.first)
                    index = j;
            }
            p = p->children[index].get();
        }
    }
    std::pair<size_t, float> ecdf1d_pair_fromgrid_trie(const std::vector<std::pair<int,int>> &row, size_t sample_size, size_t ind, float val01) const
    {
        size_t l = 0, r = grids[ind].size() - 1;

        size_t m = 0, index1 = 0, index2 = 0;
        float cdf1, cdf2;

        while(l <= r)
        {
            m = l + (r - l) / 2;

            index1 = 0;
            index2 = 0;
            for(size_t i = 0, n = row.size(); i != n; ++i)
            {
                if(row[i].first < m)
                {
                    index1 += row[i].second;
                }
                if(row[i].first < m + 1)
                {
                    index2 += row[i].second;
                }
            }
            cdf1 = index1/float(sample_size);
            cdf2 = index2/float(sample_size);

            if((val01 > cdf1) && (val01 < cdf2))
                break;

            if(val01 > cdf1)
                l = m + 1;
            else
                r = m - 1;
        }

        float x0 = grids[ind][m], y0 = cdf1, x1 = grids[ind][m + 1], y1 = cdf2;
        return std::make_pair(m, x0 + (val01 - y0) * (x1 - x0) / (y1 - y0));
    }
};

#endif
