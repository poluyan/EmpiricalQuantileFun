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
#ifndef TEST2D_H
#define TEST2D_H

#include <kquantile.h>
#include <trie.h>
#include <trie_based.h>

#include <utility/data_io.h>
#include <utility/timer.h>

void test_2d1();
void test_2d2();

void grid_test_2d();
void test_2d_func();

void worst_space_test_grid_2d();

void test_2d_uniform_vs_nonuniform();

void test_2d_uniform_vs_nonuniform_trie();


void test_2d_discrete();


template<class InputIt, class T>
void parallel_step(InputIt first,
                   InputIt last,
                   std::vector<std::vector<T>> &pdf,
                   const T lb,
                   const T ub,
                   std::function<T(const std::vector<T> &)> f)
{
	T es = ub - lb;
	for(auto it = first; it != last; ++it)
	{
		size_t i = std::distance(pdf.begin(), it);
		T a = lb + i*es/(pdf.size() - 1);
		for(size_t j = 0; j != it->size(); j++)
		{
			std::vector<T> temp = {a, lb + j*es/(pdf.size() - 1)};
			pdf[i][j] = f(temp);
		}
	}
}

template<class InputIt, class T>
void parallel_pdf(InputIt first,
                  InputIt last,
                  std::vector<std::vector<T>> &pdf,
                  const T lb,
                  const T ub,
                  std::function<T(const std::vector<T> &)> f)
{
	const auto size = last - first;
	const auto nthreads = std::thread::hardware_concurrency();
	const auto size_per_thread = size / nthreads;

	std::vector<std::future<void>> futures;
	for(unsigned int i = 0; i < nthreads - 1; i++)
	{
		futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, &pdf, lb, ub, f]()
		{
			T es = ub - lb;
			for(auto it = start; it != start + size_per_thread; ++it)
			{
				size_t i = std::distance(pdf.begin(), it);
				T a = lb + i*es/(pdf.size() - 1);
				for(size_t j = 0; j != it->size(); j++)
				{
					std::vector<T> temp = {a, lb + j*es/(pdf.size() - 1)};
					pdf[i][j] = f(temp);
				}
			}
		}));
	}
	futures.emplace_back(
	  std::async([start = first + (nthreads - 1) * size_per_thread, last, &pdf, lb, ub, f]()
	{
		T es = ub - lb;
		for(auto it = start; it != last; ++it)
		{
			size_t i = std::distance(pdf.begin(), it);
			T a = lb + i*es/(pdf.size() - 1);
			for(size_t j = 0; j != it->size(); j++)
			{
				std::vector<T> temp = {a, lb + j*es/(pdf.size() - 1)};
				pdf[i][j] = f(temp);
			}
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
}

template <typename T>
void print_pdf(std::string fname, std::function<T(const std::vector<T> &)> f, const T lb, const T ub, const size_t n)
{
	std::vector<std::vector<T> > pdf(n, std::vector<T>(n));
	//parallel_step(pdf.begin(), pdf.end(), pdf, lb, ub, f);
	parallel_pdf(pdf.begin(), pdf.end(), pdf, lb, ub, f);
	data_io::write_default2d(fname, pdf, 5);
}


template <typename T>
void test_2d_kquantile()
{
	const size_t kt = 0;
	const size_t dim = 2;
	const size_t N = 1000;
	const size_t nrolls = 10000;
	const T threshold = 1e-2;
	const size_t multi = 1000000;
	std::vector<size_t> gridn(dim, N);

//    auto sample = std::make_shared<trie_based::TrieBased<trie_based::NodeCount<int>,int>>();
//    sample->set_dimension(gridn.size());
//    sample->insert(std::vector<int>{5,5});
//    sample->insert(std::vector<int>{4,5});
//    sample->insert(std::vector<int>{5,4});
//    sample->insert(std::vector<int>{1,3});



	auto sample = std::make_shared<mveqf::trie_based::Trie<mveqf::trie_based::NodeCount<int>,int>>();
	sample->set_dimension(gridn.size());
	// 10x10
	/*sample->insert(std::vector<int>{3,2}, 1);
	sample->insert(std::vector<int>{2,3}, 1);
	sample->insert(std::vector<int>{3,4}, 1);
	sample->insert(std::vector<int>{4,3}, 1);
	sample->insert(std::vector<int>{3,3}, 2);*/

	///


	auto kdesample = std::make_shared<std::vector<std::vector<T>>>();

	///first
//    sample->push_back(std::vector<T>{0.0, 0.0});

	/// second
	// 1125, 1126
	/*std::vector<std::vector<size_t> > dots =
	{
	    {484,197},
	    {731,256},
	    {582,288},
	    {505,291},
	    {505,378},
	    {699,387},
	    {799,407},
	    {793,439},
	    {591,444},
	    {736,452},
	    {751,460},
	    {308,463},
	    {170,472},
	    {243,483},
	    {406,496},
	    {388,498},
	    {442,498},
	    {691,511},
	    {485,516},
	    {598,518},
	    {522,529},
	    {32,532},
	    {472,534},
	    {578,535},
	    {596,549},
	    {566,558},
	    {479,579},
	    {209,587},
	    {858,594},
	    {601,621},
	    {787,623},
	    {749,635},
	    {503,637},
	    {878,644},
	    {428,651},
	    {585,652},
	    {494,654},
	    {687,654},
	    {803,667},
	    {650,684},
	    {595,700},
	    {315,702},
	    {549,702},
	    {561,723},
	    {596,731},
	    {625,735},
	    {664,742},
	    {731,744},
	    {732,788},
	    {736,863}
	};
	for(const auto & i : dots)
	{
	    std::vector<T> t = {-1.75 + 3.5*(i.front())/T(1125), -1.75 + 3.5*(T(1126) - i.back())/T(1126)};
	    std::cout << std::scientific << t.front() << '\t' << t.back() << std::endl;
	    kdesample->push_back(t);
	}*/

	std::ifstream gridIn;
	gridIn.open("maps/Rpics/sampled_implicit.dat");
	if(!gridIn.is_open())
	{
		std::cout << "Error opening file." << std::endl;
		return;
	}
	std::vector<T> lbv(dim, std::numeric_limits<T>::max());
	std::vector<T> ubv(dim, std::numeric_limits<T>::min());
	std::vector<T> t(dim);
	while(!gridIn.eof())
	{
		gridIn >> t[0];
		gridIn >> t[1];

		for(size_t i = 0; i != t.size(); i++)
		{
			if(t[i] < lbv[i])
				lbv[i] = t[i];
			if(t[i] > ubv[i])
				ubv[i] = t[i];
		}

		kdesample->push_back(t);
		if(kdesample->size() > 10000)
			break;
	}
	gridIn.close();


	const T lb = -2, ub = 2;

	mveqf::kde::KDE<T> obj;
	obj.set_dimension(dim);
	obj.set_kernel_type(kt);
	obj.set_sample_shared(kdesample);

	std::function<T(const std::vector<T> &)> fff = std::bind(&mveqf::kde::KDE<T>::pdf, std::ref(obj), std::placeholders::_1);

	print_pdf<T>("maps/Trie/kpdf2d_init.dat", fff, lb, ub, 1000);

	///

	std::vector<std::vector<T>> grids(gridn.size());
	std::vector<T> dx(gridn.size());

	for(size_t i = 0; i != grids.size(); i++)
	{
		size_t num_points = gridn[i];
		std::vector<T> onegrid(num_points);
		T es = ub - lb;
		for(size_t j = 0; j != onegrid.size(); j++)
		{
			onegrid[j] = lb + j*es/(num_points);
		}
		grids[i] = onegrid;
		dx[i] = es/(num_points*2);
	}

	// finding start dot over created grid
	std::vector<std::vector<int>> points;
	for(auto it = kdesample->begin(); it != kdesample->end(); ++it)
	{
		std::vector<int> startdot(dim);
		for(size_t i = 0; i != startdot.size(); i++)
		{
//            std::vector<T> val(grids[i].size());
//            for(size_t j = 0; j != val.size(); j++)
//            {
//                val[j] = grids[i][j] + dx[i];
//            }
//            auto pos1 = std::lower_bound(val.begin(), val.end(), (*it)[i]);
//            startdot[i] = std::distance(val.begin(), pos1) - 1;
			auto pos1 = std::lower_bound(grids[i].begin(), grids[i].end(), (*it)[i]);
			if(pos1 == grids[i].end())
				startdot[i] = grids[i].size() - 2;
			else if(pos1 == grids[i].begin())
				startdot[i] = 0;
			else
				startdot[i] = std::distance(grids[i].begin(), pos1) - 1;
		}
//        std::cout << startdot.front() << '\t' << startdot.back() << std::endl;
//        std::cin.get();
		points.push_back(startdot);
	}

	std::set<std::vector<int>> visited;

	//std::vector<std::vector<T> > samples;

	std::function<T(const std::vector<T> &)> pdffunc = std::bind(&mveqf::kde::KDE<T>::pdf, std::ref(obj), std::placeholders::_1);
//    std::cout << pdffunc(std::vector<float>{-1.0,1.0}) << std::endl;
//    std::cout << obj.pdf(std::vector<float>{-1.0,1.0}) << std::endl;
//
//    std::cin.get();

	size_t counter = 0, fe_count = 0;
	std::cout << "FloodFill start" << std::endl;
	mveqf::mvff::FloodFill_MultipleGrids_VonNeumann_trie<T>(grids, points, sample, dx, counter, fe_count, pdffunc, threshold, multi);
	std::cout << "FloodFill end" << std::endl;

	//std::cout << samples.size() << std::endl;
	//std::cin.get();
	///

	std::mt19937_64 generator;
	generator.seed(1);
	std::uniform_real_distribution<T> ureal01(0.0,1.0);

	std::vector<std::vector<T> > sampled;
	std::vector<std::vector<T> > values01;

	mveqf::ImplicitTrieKQuantile<int, T> quant(std::vector<T>(gridn.size(), lb), std::vector<T>(gridn.size(), ub), gridn, kt);
	quant.set_sample_shared(sample);

	timer::Timer time_cpp11;
	time_cpp11.reset();

	std::vector<T> temp1(gridn.size());
	std::vector<T> temp2(temp1.size());

	std::cout << "Transform start" << std::endl;
	for(size_t i = 0; i != nrolls; ++i)
	{
//        timer::Timer t;
		for(size_t j = 0; j != temp1.size(); j++)
		{
			temp1[j] = ureal01(generator);
		}
		quant.transform(temp1,temp2);
		values01.push_back(temp1);
		sampled.push_back(temp2);
//        std::cout << std::fixed << t.elapsed_seconds() << std::endl;
	}
	std::cout << "Transform end" << std::endl;
	/*for(size_t i = 0; i != temp1.size(); ++i)
	{
	    for(size_t j = 0; j != temp1.size(); j++)
	    {
	        temp1[j] = ureal01(generator);
	    }
	    temp1[i] = 0.0;
	    quant.ktransform(temp1,temp2);
	    values01.push_back(temp1);
	    sampled.push_back(temp2);
	}
	for(size_t i = 0; i != temp1.size(); ++i)
	{
	    for(size_t j = 0; j != temp1.size(); j++)
	    {
	        temp1[j] = ureal01(generator);
	    }
	    temp1[i] = 1.0;
	    quant.ktransform(temp1,temp2);
	    values01.push_back(temp1);
	    sampled.push_back(temp2);
	}

	for(size_t i = 0; i != temp1.size(); ++i)
	{
	    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
	    {
	        temp1[j] = 0.0;
	    }
	    for(size_t j = i + 1; j < temp1.size(); j++)
	    {
	        temp1[j] = ureal01(generator);
	    }
	    quant.ktransform(temp1,temp2);
	    values01.push_back(temp1);
	    sampled.push_back(temp2);
	}
	for(size_t i = 0; i != temp1.size(); ++i)
	{
	    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
	    {
	        temp1[j] = 1.0;
	    }
	    for(size_t j = i + 1; j < temp1.size(); j++)
	    {
	        temp1[j] = ureal01(generator);
	    }
	    quant.ktransform(temp1,temp2);
	    values01.push_back(temp1);
	    sampled.push_back(temp2);
	}

	for(size_t i = 0; i != temp1.size() - 1; ++i)
	{
	    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
	    {
	        temp1[j] = ureal01(generator);
	    }
	    for(size_t j = i + 1; j < temp1.size(); j++)
	    {
	        temp1[j] = 0.0;
	    }
	    quant.ktransform(temp1,temp2);
	    values01.push_back(temp1);
	    sampled.push_back(temp2);
	}
	for(size_t i = 0; i != temp1.size() - 1; ++i)
	{
	    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
	    {
	        temp1[j] = ureal01(generator);
	    }
	    for(size_t j = i + 1; j < temp1.size(); j++)
	    {
	        temp1[j] = 1.0;
	    }
	    quant.ktransform(temp1,temp2);
	    values01.push_back(temp1);
	    sampled.push_back(temp2);
	}

	for(size_t i = 0; i != sampled.size(); i++)
	{
	    for(size_t j = 0; j != sampled[i].size(); j++)
	    {
	        if(sampled[i][j] < lb || sampled[i][j] > ub)
	        {
	            std::cout << "out of bounds!" << std::endl;
	        }
	    }
	}*/

	data_io::write_default2d("maps/sampled_implicit.dat", sampled, 5);

	std::shared_ptr<std::vector<std::vector<T> >> ss = std::make_shared<std::vector<std::vector<T>>>();
	for(const auto & i : sampled)
		ss->push_back(i);

	mveqf::kde::KDE<T> obj1;
	obj1.set_dimension(dim);
	obj1.set_kernel_type(kt);
	obj1.set_sample_shared(ss);
	std::function<T(const std::vector<T> &)> ff = std::bind(&mveqf::kde::KDE<T>::pdf, std::ref(obj1), std::placeholders::_1);

	print_pdf<T>("maps/Trie/kpdf2d.dat", ff, lb, ub, 1000);
}

#endif
