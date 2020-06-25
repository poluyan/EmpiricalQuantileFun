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
#ifndef TESTND_H
#define TESTND_H

#include <kquantile.h>
#include <trie.h>
#include <trie_based.h>

#include <utility/data_io.h>
#include <utility/timer.h>

namespace testNd
{

	void test_grid_10d();
	void test_Nd(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nsamples, size_t Nrolls);

	std::pair<double, double> test_Nd_time(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nsamples, size_t Nrolls);
	void grid_test_Nd();
	void dim_test_Nd();



	template <typename T>
	bool increase(const std::vector<std::vector<T>> &v, std::vector<size_t> &it)
	{
		for(size_t i = 0, size = it.size(); i != size; i++)
		{
			const size_t index = size - 1 - i;
			++it[index];
			if(it[index] == v[index].size())
			{
				it[index] = 0;
			}
			else
			{
				return true;
			}
		}
		return false;
	}

	template <typename T>
	std::vector<T> get_line(const std::vector<std::vector<T>> &v, std::vector<size_t> &it)
	{
		std::vector<T> rez(v.size());
		for(size_t i = 0, size = v.size(); i != size; i++)
		{
			rez[i] = v[i][it[i]];
		}
		return rez;
	}

	template <typename T>
	std::vector<std::vector<T>> iterate(const std::vector<std::vector<T>> &v)
	{
		std::vector<size_t> it(v.size(), 0);
		std::vector<std::vector<T>> values;
		do
		{
			values.push_back(get_line(v, it));
		}
		while(increase(v, it));
		return values;
	}

	template <typename T, typename I>
	void iterate_trie(const std::vector<std::vector<T>> &v, std::shared_ptr<mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<I>,I>> sample)
	{
		std::vector<size_t> it(v.size(), 0);
		do
		{
			sample->insert(get_line(v, it));
		}
		while(increase(v, it));
	}

	std::vector<double> worst_space(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nrolls, size_t seed_append);

	void worst_space_test_dim();
	void worst_space_test_grid();

	void worst_space_check();

	void sample_size_test(size_t dim);











	template <typename TIndex, typename TFloat>
	void explicit_quantile_descrete(TFloat lb, TFloat ub, std::vector<std::vector<TFloat>> &sample, size_t nrolls)
	{
		std::mt19937_64 generator;
		generator.seed(1);
		std::uniform_real_distribution<TFloat> ureal01(0.0,1.0);

		std::vector<std::vector<TIndex> > sampled;
		std::vector<std::vector<TFloat> > values01;

		mveqf::ExplicitQuantile<TIndex, TFloat> quant;
		quant.set_grid_from_sample(
		  std::vector<TFloat>(sample.front().size(), lb),
		  std::vector<TFloat>(sample.front().size(), ub),
		  sample,
		  1.0e-6);
		quant.set_sample(sample);
		
		auto gridn = quant.get_grid_number();
		for(const auto & i : gridn)
			std::cout << i << std::endl;

		timer::Timer time_cpp11;
		time_cpp11.reset();

		std::vector<TFloat> temp1(sample.front().size());
		std::vector<TIndex> temp2(temp1.size());

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

//    std::vector<float> a = {7.546996e-01,
//                            2.108486e-01,
//                            9.091859e-01,
//                            9.425417e-01,
//                            6.641316e-01,
//                            8.285044e-01,
//                            6.402776e-01,
//                            1.946763e-01,
//                            2.501838e-01,
//                            8.654959e-01
//                           };
//    std::vector<float> b(a.size(), 0);
//    quant.transform(a, b);
//    sampled.push_back(b);

//    a = {0.68};
//    b = {0.0};
//    quant.transform(a, b);
//    sampled.push_back(b);

		std::cout << "total time: " << time_cpp11.elapsed_seconds() << std::endl;
		std::cout << "time per transform: " << std::scientific << time_cpp11.elapsed_seconds()/double(sampled.size()) << std::endl;
		data_io::write_default2d("maps/sampled_explicit_discrete.dat", sampled, 10);
		//data_io::write_default2d("maps/z.dat", values01, 15);
	}


	void test_Nd_discrete();
}
#endif
