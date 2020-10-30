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
#ifndef MVFF_H
#define MVFF_H

#include <trie.h>
#include <vector>
#include <set>

namespace mveqf
{
	namespace mvff
	{
		template <typename T, typename V>
		struct cell
		{
			std::vector<T> dot;
			V value;

			cell() {}
			cell(std::vector<T> _dot, V _value):dot(_dot),value(_value) {}
			cell(const cell &a):dot(a.dot), value(a.value) {}
			bool operator==(const cell& a) const
			{
				return dot == a.dot;
			}
			bool operator<(const cell& a) const
			{
				return dot < a.dot;
			}
		};
		template <typename T>
		void FloodFill_MultipleGrids_VonNeumann_trie(const std::vector<std::vector<T>> &grids,
		    std::vector<std::vector<int>> &pp,
		    std::shared_ptr<mveqf::trie::Trie<mveqf::trie_based::NodeCount<int>,int>> samples,
		    const std::vector<T> &dx,
		    size_t &counter,
		    size_t &fe_count,
		    std::function<T(const std::vector<T> &)> f,
		    const T threshold,
		    const size_t multi)

		{

			bool multithread = true;
			const auto nthreads = std::thread::hardware_concurrency();
//    std::cout << "std::thread::hardware_concurrency " << nthreads << std::endl;

			std::vector<T> fvalues;
			std::vector<T> dot(grids.size());
			std::map<std::vector<int>, T> fmap;

			mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<int>,int> visited;
			mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<int>,int> points;
			mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<int>,int> not_coumputed;
			mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<int>,int> ss;

			visited.set_dimension(grids.size());
			points.set_dimension(grids.size());
			not_coumputed.set_dimension(grids.size());
			ss.set_dimension(grids.size());

			for(const auto &start : pp)
			{
				points.insert(start);
//        std::cout << "init points " << '\t' << start.front() << '\t' << start.back() << std::endl;
			}
//    std::cout << std::endl;
//    std::cin.get();
			++fe_count;

			T min = std::numeric_limits<T>::max();
			T max = std::numeric_limits<T>::min();

			while(!points.empty() || !not_coumputed.empty())
			{
				while(!points.empty())
				{
					auto point = points.get_and_remove_last();
					if(visited.search(point) || ss.search(point))
					{
						continue;
					}
					ss.insert(point);

					auto init_point = point;
					for(size_t i = 0; i != point.size(); i++)
					{
						point = init_point;
						point[i] = point[i] + 1;

						if(point[i] < 0 || point[i] > static_cast<int>(grids[i].size() - 2))
						{
//                    std::cout << "here " << point[i] << std::endl;
							continue;
						}

						if(!visited.search(point) && !samples->search(point))
						{
//                    std::cout << "insert " << point.front() << '\t' << point.back() << std::endl;
							not_coumputed.insert(point);
						}
					}
					for(size_t i = 0; i != point.size(); i++)
					{
						point = init_point;
						point[i] = point[i] - 1;

						if(point[i] < 0 || point[i] > static_cast<int>(grids[i].size() - 2))
						{
//                    std::cout << "here " << point[i] << std::endl;
							continue;
						}

						if(!visited.search(point) && !ss.search(point))
						{
//                    std::cout << "insert " << point.front() << '\t' << point.back() << std::endl;
							not_coumputed.insert(point);
						}
					}
				}

				size_t number_to_points = 0;
				while(!not_coumputed.empty())
				{
					if(!multithread)
					{
						auto point = not_coumputed.get_and_remove_last();
						std::vector<T> values(point.size());
						for(size_t j = 0; j != values.size(); j++)
						{
							values[j] = grids[j][point[j]] + dx[j];
						}
						T value = f(values);

						++fe_count;

						if(value > threshold)
						{
							++number_to_points;
							points.insert(point);
							fvalues.push_back(value);
							fmap.insert(std::make_pair(point, value));
							if(min > value)
								min = value;
							if(max < value)
								max = value;
						}
						else
						{
							visited.insert(point);
						}
					}
					else
					{

						std::vector<cell<int,T>> to_compute;
						for(size_t i = 0; i != nthreads*10000000; i++)
						{
							if(not_coumputed.empty())
								break;

							auto point = not_coumputed.get_and_remove_last();
							to_compute.push_back(cell<int,T>(point, -1));
						}


						auto first = to_compute.begin();
						auto last = to_compute.end();
						const auto size = last - first;
						const auto size_per_thread = size / nthreads;

						std::vector<std::future<void>> futures;
						for(unsigned int i = 0; i < nthreads - 1; i++)
						{
							futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, &to_compute, &grids, &dx, f]()
							{
								for(auto it = start; it != start + size_per_thread; ++it)
								{
									size_t i = std::distance(to_compute.begin(), it);
									std::vector<T> values(to_compute[i].dot.size());
									for(size_t j = 0; j != values.size(); j++)
									{
										values[j] = grids[j][to_compute[i].dot[j]] + dx[j];
									}
									to_compute[i].value = f(values);
								}
							}));
						}
						futures.emplace_back(
						  std::async([start = first + (nthreads - 1) * size_per_thread, last, &to_compute, &grids, &dx, f]()
						{
							for(auto it = start; it != last; ++it)
							{
								size_t i = std::distance(to_compute.begin(), it);
								std::vector<T> values(to_compute[i].dot.size());
								for(size_t j = 0; j != values.size(); j++)
								{
									values[j] = grids[j][to_compute[i].dot[j]] + dx[j];
								}
								to_compute[i].value = f(values);
							}

						}));

						for(auto &&future : futures)
						{
							if(future.valid())
							{
								future.get();
							}
							else
							{
								throw std::runtime_error("Something going wrong.");
							}
						}
						for(size_t i = 0; i != to_compute.size(); i++)
						{
							if(to_compute[i].value > threshold)
							{
								points.insert(to_compute[i].dot);
								fvalues.push_back(to_compute[i].value);
								fmap.insert(std::make_pair(to_compute[i].dot, to_compute[i].value));
								++number_to_points;

								if(min > to_compute[i].value)
									min = to_compute[i].value;
								if(max < to_compute[i].value)
									max = to_compute[i].value;
							}
							else
							{
								visited.insert(to_compute[i].dot);
							}
						}
					}
				}
				//std::cout << number_to_points << std::endl;
			}
			counter = 0;
			visited.remove_tree();
			points.remove_tree();
			not_coumputed.remove_tree();

//    T delta = std::numeric_limits<T>::max();
//    std::sort(fvalues.begin(), fvalues.end());
//    for(size_t j = 1; j != fvalues.size(); j++)
//    {
//        T temp = fvalues[j] - fvalues[j - 1];
//        if(temp > 0 && delta > temp)
//            delta = temp;
//    }
//    std::cout << std::scientific << delta << std::endl;
//
//    delta = 1.0/delta;

			/*while(!ss.empty())
			{
			    std::vector<cell<int,T>> to_compute;
			    for(size_t i = 0; i != nthreads*10000000; i++)
			    {
			        if(ss.empty())
			            break;

			        auto point = ss.get_and_remove_last();
			        to_compute.push_back(cell<int,T>(point, -1));
			    }


			    auto first = to_compute.begin();
			    auto last = to_compute.end();
			    const auto size = last - first;
			    const auto size_per_thread = size / nthreads;

			    std::vector<std::future<void>> futures;
			    for(unsigned int i = 0; i < nthreads - 1; i++)
			    {
			        futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, &to_compute, &grids, &dx, f]()
			        {
			            for(auto it = start; it != start + size_per_thread; ++it)
			            {
			                size_t i = std::distance(to_compute.begin(), it);
			                std::vector<T> values(to_compute[i].dot.size());
			                for(size_t j = 0; j != values.size(); j++)
			                {
			                    values[j] = grids[j][to_compute[i].dot[j]] + dx[j];
			                }
			                to_compute[i].value = f(values);
			            }
			        }));
			    }
			    futures.emplace_back(
			        std::async([start = first + (nthreads - 1) * size_per_thread, last, &to_compute, &grids, &dx, f]()
			    {
			        for(auto it = start; it != last; ++it)
			        {
			            size_t i = std::distance(to_compute.begin(), it);
			            std::vector<T> values(to_compute[i].dot.size());
			            for(size_t j = 0; j != values.size(); j++)
			            {
			                values[j] = grids[j][to_compute[i].dot[j]] + dx[j];
			            }
			            to_compute[i].value = f(values);
			        }
			    }));

			    for(auto &&future : futures)
			    {
			        if(future.valid())
			        {
			            future.get();
			        }
			        else
			        {
			            throw std::runtime_error("Something going wrong.");
			        }
			    }
			    for(size_t i = 0; i != to_compute.size(); i++)
			    {
			        //samples->insert(to_compute[i].dot, std::round(delta*to_compute[i].value) + 1);
			        to_compute[i].value = (to_compute[i].value - min)/(max - min);
			        samples->insert(to_compute[i].dot, 1000*to_compute[i].value + 1);
			    }
			}*/
//    std::cout << ss.get_total_count() << std::endl;
//    std::cout << fmap.size() << std::endl;
			size_t count = 0;
//    for(const auto & i : fmap)
//    {
//        if(ss.search(i.first))
//        {
//            T val = (i.second - min)/(max - min);
//            samples->insert(i.first, 1000*val + 1);
//            ++count;
//        }
//    }
			while(!ss.empty())
			{
				auto point = ss.get_and_remove_last();
				auto it = fmap.find(point);

//        for(size_t i = 0; i != point.size(); i++)
//        {
//            if(point[i] < 0 || point[i] > static_cast<int>(grids[i].size() - 2))
//            {
//                std::cout << "what??" << std::endl;
//            }
//        }
//        std::vector<T> values(point.size());
//        for(size_t j = 0; j != values.size(); j++)
//        {
//            values[j] = grids[j][point[j]] + dx[j];
//            if(values[j] > grids[j].back())
//            {
//                std::cout << "wow!! " << values[j] << std::endl;
//            }
//        }



				if(it != fmap.end())
				{
					T val = (it->second - min)/(max - min);
					if(!samples->search(point))
					{
						samples->insert(it->first, multi*val + 1);
						++count;
					}

//            for(size_t i = 0; i != point.size(); i++)
//            {
//                std::cout << point[i] << '\t';
//            }
//            std::cout << std::endl;
				}
			}
//    std::cout << count << '\t' << grids.front().size() << '\t' << grids.back().size() << std::endl;


			if(!count)
			{
				for(size_t i = 0; i != pp.size(); i++)
					if(!samples->search(pp[i]))
						samples->insert(pp[i], 1);
			}

		}

	}

}

#endif
