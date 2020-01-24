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
#ifndef MVEQF_H
#define MVEQF_H

#include "kquantile.h"
#include "mvff.h"

namespace mveqf
{

template <typename T, typename U>
class MVEQF
{
protected:
    size_t kernel_type;
    std::shared_ptr<trie_based::Trie<trie_based::NodeCount<T>,T>> sample;
    std::shared_ptr<ImplicitTrieKQuantile<T, U>> qf;
public:
    MVEQF(): kernel_type(0) {}
    void set_kernel_type(size_t kt)
    {
        kernel_type = kt;
    }
    void set_sample_and_bounds(const std::vector<std::vector<U>> &init_points, const std::vector<U> &lb, const std::vector<U> &ub)
    {
        const size_t dimension = lb.size();

        const size_t N = 500;
        const U threshold = 1e-8;
        const size_t multi = 100000;

        auto init_sample_points = std::make_shared<std::vector<std::vector<U>>>(init_points.size(), std::vector<U>(dimension));
        for(auto it = init_sample_points->begin(); it != init_sample_points->end(); ++it)
        {
            size_t i = std::distance(init_sample_points->begin(), it);
            for(size_t j = 0; j != dimension; j++)
            {
                (*it)[j] = init_points[i][j];
            }
        }

        mveqf::kde::KDE<U> init_pdf;
        init_pdf.set_dimension(dimension);
        init_pdf.set_kernel_type(kernel_type);
        init_pdf.set_sample_shared(init_sample_points);

        std::vector<size_t> gridn(lb.size(), N);
        qf = std::make_shared<ImplicitTrieKQuantile<T, U>>(lb, ub, gridn, kernel_type);

        const auto grid = qf->get_grid();
        const auto dx = qf->get_dx();

        // finding start dot over created grid
        std::vector<std::vector<T>> points;
        for(const auto &k : init_points)
        {
            std::vector<int> startdot(dimension);
            for(size_t i = 0; i != startdot.size(); i++)
            {
                std::vector<U> val(grid[i].size());
                for(size_t j = 0; j != val.size(); j++)
                {
                    val[j] = grid[i][j] + dx[i];
                }
                auto pos1 = std::lower_bound(val.begin(), val.end(), k[i]);
                startdot[i] = std::distance(val.begin(), pos1) - 1;
            }
            points.push_back(startdot);
        }

        std::set<std::vector<T>> visited;

        std::function<U(const std::vector<U> &)> pdffunc = std::bind(&mveqf::kde::KDE<U>::pdf, std::ref(init_pdf), std::placeholders::_1);

        sample = std::make_shared<mveqf::trie_based::Trie<mveqf::trie_based::NodeCount<T>,T>>();
        sample->set_dimension(dimension);
//
        size_t counter = 0, fe_count = 0;
        std::cout << "FloodFill" << std::endl;
        mveqf::mvff::FloodFill_MultipleGrids_VonNeumann_trie<U>(grid, points, sample, dx, counter, fe_count, pdffunc, threshold, multi);
        std::cout << "FloodFill" << std::endl;

        qf->set_sample_shared(sample);
    }
    void transform(const std::vector<U>& in01, std::vector<U>& out) const
    {
        qf->transform(in01, out);
    }
    std::vector<U> transform(const std::vector<U>& in01) const
    {
        return qf->transform(in01);
    }
};

}

#endif
