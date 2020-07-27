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

#ifndef TEST_H
#define TEST_H

#include <trie_based.h>
#include <explicit.h>
#include <implicit.h>

#include <sdot.h>

#include <utility/data_io.h>
#include <utility/timer.h>

#include <random>

void explicit_quantile(float lb, float ub, std::vector<size_t> gridn, std::vector<std::vector<int> > &sample, size_t nrolls);
void implicit_quantile_class(float lb, float ub, std::vector<size_t> gridn, std::vector<std::vector<int> > &sample, size_t nrolls);
void implicit_quantile_class_sorted(float lb, float ub, std::vector<size_t> gridn, std::vector<std::vector<int> > &sample, size_t nrolls);
void implicit_quantile_class_sorted_interp(float lb, float ub, std::vector<size_t> gridn, std::vector<std::vector<int> > &sample, size_t nrolls);
void implicit_quantile_graph_sorted(float lb, float ub, std::vector<size_t> gridn, size_t nrolls);

template <typename TFloat, typename TIndex>
void check_transform(const mveqf::Quantile<TIndex, TFloat> &quant, const std::vector<size_t> &gridn, std::string fname, size_t nrolls)
{
	std::mt19937_64 generator;
	generator.seed(1);
	std::uniform_real_distribution<TFloat> ureal01(0.0,1.0);

	std::vector<std::vector<TFloat> > sampled;
	std::vector<std::vector<TFloat> > values01;

	timer::Timer time_cpp11;
	time_cpp11.reset();

	std::vector<TFloat> temp1(gridn.size());
	std::vector<TFloat> temp2(temp1.size());

	for(size_t i = 0; i != nrolls; ++i)
	{
		for(size_t j = 0; j != temp1.size(); j++)
		{
			temp1[j] = ureal01(generator);
		}
		quant.transform(temp1,temp2);
		values01.push_back(temp1);
		sampled.push_back(temp2);
	}

	for(size_t i = 0; i != temp1.size(); ++i)
	{
		for(size_t j = 0; j != temp1.size(); j++)
		{
			temp1[j] = ureal01(generator);
		}
		temp1[i] = 0.0;
		quant.transform(temp1,temp2);
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
		quant.transform(temp1,temp2);
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
		quant.transform(temp1,temp2);
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
		quant.transform(temp1,temp2);
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
		quant.transform(temp1,temp2);
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
		quant.transform(temp1,temp2);
		values01.push_back(temp1);
		sampled.push_back(temp2);
	}

	std::cout << "total time: " << time_cpp11.elapsed_seconds() << std::endl;
	std::cout << "time per transform: " << std::scientific << time_cpp11.elapsed_seconds()/double(sampled.size()) << std::endl;
	data_io::write_default2d(fname, sampled, 10);
	//data_io::write_default2d("maps/z.dat", values01, 15);
}


template <typename TFloat, typename TIndex>
void perform_tests(TFloat lb, TFloat ub, const std::vector<size_t> &gridn, const std::vector<std::vector<TIndex> > &sample, size_t nrolls)
{
	mveqf::ExplicitQuantile<TIndex, TFloat> qex(
	  std::vector<TFloat>(gridn.size(), lb),
	  std::vector<TFloat>(gridn.size(), ub),
	  gridn);
	qex.set_sample(sample);
	check_transform<TFloat, TIndex>(qex, gridn, "maps/sampled_explicit.dat", nrolls);

	mveqf::ImplicitQuantile<TIndex, TFloat> qim(
	  std::vector<TFloat>(gridn.size(), lb),
	  std::vector<TFloat>(gridn.size(), ub),
	  gridn);
	qim.set_sample(sample);
	check_transform<TFloat, TIndex>(qim, gridn, "maps/sampled_implicit.dat", nrolls);

	mveqf::ImplicitQuantileSorted<TIndex, TFloat> qims(
	  std::vector<TFloat>(gridn.size(), lb),
	  std::vector<TFloat>(gridn.size(), ub),
	  gridn);
	qims.set_sample(sample);
	check_transform<TFloat, TIndex>(qims, gridn, "maps/sampled_implicit_sorted.dat", nrolls);

	mveqf::ImplicitQuantileSortedInterp<TIndex, TFloat> qimsi(
	  std::vector<TFloat>(gridn.size(), lb),
	  std::vector<TFloat>(gridn.size(), ub),
	  gridn);
	qimsi.set_sample(sample);
	check_transform<TFloat, TIndex>(qimsi, gridn, "maps/sampled_implicit_sorted_interp.dat", nrolls);
}

template<typename TFloat>
void ot_1d()
{
	const int n = 100; // number of init dots
	const int dimension = 2;
	std::vector<std::vector<TFloat>> init_dots(n, std::vector<TFloat>(dimension));
	std::vector<size_t> weights(n);

	std::mt19937_64 generator;
	generator.seed(1);
	std::uniform_real_distribution<TFloat> ureal(-100.0,100.0);
	std::uniform_real_distribution<TFloat> ureal01(0.0,1.0);
	for(size_t i = 0; i != n; i++)
	{
		for(size_t j = 0; j != dimension; j++)
		{
			init_dots[i][j] = ureal(generator);
			weights[i] = std::abs(ureal(generator))*100*ureal01(generator) + 1;
		}
	}
	mveqf::ot::SDOT<TFloat> obj(init_dots, weights);
	std::vector<size_t> counter(init_dots.size(), 0);
	std::vector<size_t> times(init_dots.size());
	std::vector<std::vector<size_t>> unique_values;
//	std::vector<std::vector<TFloat>> sampled;
	for(size_t i = 0; i != 5000; i++)
	{
		std::vector<TFloat> temp(dimension);
		for(size_t j = 0; j != dimension; j++)
		{
			temp[j] = ureal01(generator);
		}


		++counter[obj.transform_index(temp)];
		//std::cout << obj.transform_index(temp) << std::endl;
		auto t = obj.transform(temp);//obj.transform_float(temp);

		auto it = find_if(unique_values.begin(), unique_values.end(), [&t](const std::vector<size_t>& s)
		{
			size_t c = 0;
			for(size_t j = 0; j != s.size(); j++)
			{
				if(t[j] == s[j])
					++c;
			}
			return s.size() == c;
		});
		if(it == unique_values.end())
		{
			unique_values.push_back(t);
			times[unique_values.size() - 1] = 1;
		}
		else
			++times[std::distance(unique_values.begin(), it)];

//		for(size_t j = 0; j != dimension; j++)
//			std::cout << t[j] << '\t';
//		std::cout << std::endl;


//		sampled.push_back(t);
	}
	std::cout << std::endl;
	for(size_t i = 0; i != counter.size(); i++)
		std::cout << counter[i] << '\t' << times[i] << std::endl;
	std::cout << std::endl;
	std::cout << unique_values.size() << std::endl;

	data_io::write_default2d("maps/ot/sample.dat", init_dots, 10);
//	data_io::write_default2d("maps/ot/sampled.dat", sampled, 10);
}

#endif
