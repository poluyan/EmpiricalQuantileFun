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
#ifndef RUN_TEST_H
#define RUN_TEST_H


#include <utility/timer.h>
#include <trie_based.h>
#include <explicit.h>
#include <implicit.h>
#include <test/test.h>
#include <test/test1d.h>
#include <test/test2d.h>
#include <test/test3d.h>
#include <test/testNd.h>
#include <utility/data_io.h>
#include <test/test_kde.h>
#include <kquantile.h>
#include <mveqf.h>

#include <random>
#include <vector>

void test_mveqf()
{
	std::vector<std::vector<double>> sample =
	{
		{-2.442222e-01, 1.137655e+00},
		{5.242222e-01, 9.542629e-01},
		{6.066667e-02, 8.547957e-01},
		{-1.788889e-01, 8.454707e-01},
		{-1.788889e-01, 5.750444e-01},
		{4.246667e-01, 5.470693e-01},
		{7.357778e-01, 4.849023e-01},
		{7.171111e-01, 3.854352e-01},
		{8.866667e-02, 3.698934e-01},
		{5.397778e-01, 3.450266e-01},
		{5.864444e-01, 3.201599e-01},
		{-7.917778e-01, 3.108348e-01},
		{-1.221111e+00, 2.828597e-01},
		{-9.940000e-01, 2.486679e-01},
		{-4.868889e-01, 2.082593e-01},
		{-5.428889e-01, 2.020426e-01},
		{-3.748889e-01, 2.020426e-01},
		{3.997778e-01, 1.616341e-01},
		{-2.411111e-01, 1.460924e-01},
		{1.104444e-01, 1.398757e-01},
		{-1.260000e-01, 1.056838e-01},
		{-1.650444e+00, 9.635879e-02},
		{-2.815556e-01, 9.014210e-02},
		{4.822222e-02, 8.703375e-02},
		{1.042222e-01, 4.351687e-02},
		{1.088889e-02, 1.554174e-02},
		{-2.597778e-01, -4.973357e-02},
		{-1.099778e+00, -7.460036e-02},
		{9.193333e-01, -9.635879e-02},
		{1.197778e-01, -1.802842e-01},
		{6.984444e-01, -1.865009e-01},
		{5.802222e-01, -2.238011e-01},
		{-1.851111e-01, -2.300178e-01},
		{9.815556e-01, -2.517762e-01},
		{-4.184444e-01, -2.735346e-01},
		{7.000000e-02, -2.766430e-01},
		{-2.131111e-01, -2.828597e-01},
		{3.873333e-01, -2.828597e-01},
		{7.482222e-01, -3.232682e-01},
		{2.722222e-01, -3.761101e-01},
		{1.011111e-01, -4.258437e-01},
		{-7.700000e-01, -4.320604e-01},
		{-4.200000e-02, -4.320604e-01},
		{-4.666667e-03, -4.973357e-01},
		{1.042222e-01, -5.222025e-01},
		{1.944444e-01, -5.346359e-01},
		{3.157778e-01, -5.563943e-01},
		{5.242222e-01, -5.626110e-01},
		{5.273333e-01, -6.993783e-01},
		{5.397778e-01, -9.325044e-01}
	};
	data_io::write_default2d("maps/sample2d.dat", sample, 5);

	mveqf::MVEQF<int, double> qf;
	qf.set_kernel_type(0);
	std::vector<size_t> empty;
	qf.set_sample_and_bounds(sample, std::vector<double>(2, -2), std::vector<double>(2, 2), empty, 500, 1e-8, 100000);

	std::mt19937_64 generator;
	generator.seed(1);
	std::uniform_real_distribution<double> ureal01(0.0,1.0);
	size_t nrolls = 500;
	std::vector<std::vector<double>> sampled(nrolls, std::vector<double>(2));
	for(size_t i = 0; i != nrolls; i++)
	{
		std::vector<double> t = {ureal01(generator), ureal01(generator)};
		sampled[i] = qf.transform(t);
	}
	data_io::write_default2d("maps/sampled2d.dat", sampled, 5);


	nrolls = 50000;
	sampled = std::vector<std::vector<double>>(nrolls, std::vector<double>(2));
	auto vals01 = sampled;
	for(auto & i : vals01)
	{
		for(auto & j : i)
		{
			j = ureal01(generator);
		}
	}

	const auto nthreads = std::thread::hardware_concurrency();
	auto first = sampled.begin();
	auto last = sampled.end();
	const auto size = last - first;
	const auto size_per_thread = size / nthreads;

	std::vector<std::future<void>> futures;
	for(unsigned int i = 0; i < nthreads - 1; i++)
	{
		futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, &vals01, &sampled, &qf]()
		{
			for(auto it = start; it != start + size_per_thread; ++it)
			{
				qf.transform(vals01[std::distance(sampled.begin(), it)], *it);
			}
		}));
	}
	futures.emplace_back(
	  std::async([start = first + (nthreads - 1) * size_per_thread, last, &vals01, &sampled, &qf]()
	{
		for(auto it = start; it != last; ++it)
		{
			qf.transform(vals01[std::distance(sampled.begin(), it)], *it);
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
	data_io::write_default2d("maps/sampled2d_1m.dat", sampled, 5);
}


void test_mveqf2()
{
	std::vector<std::vector<double>> sample =
	{
		{-8.4772980213e-01,	-1.3907787800e+00},
		{-6.2346875668e-02,	3.7843622267e-02},
		{-3.8475060463e-01,	1.5340744257e+00},
		{7.1187913418e-03,	1.3396510482e-01},
		{3.4900745749e-01,	7.0569372177e-01},
		{-9.4249987602e-01,	4.6853682399e-01},
		{8.1792432070e-01,	-5.3509891033e-01},
		{-1.7806731164e-01,	4.4960036874e-01},
		{-5.1068866253e-01,	1.2097091675e+00},
		{2.0778149366e-02,	4.8589122295e-01},
		{-5.2311074734e-01,	1.0469725132e+00},
		{-3.7779271603e-02,	5.5113613605e-01},
		{-4.4691383839e-01,	-1.4604777098e+00},
		{-8.7871319056e-01,	-1.5926430225e+00},
		{6.1549007893e-01,	7.4339038134e-01},
		{8.1910520792e-01,	-2.2428214550e-02},
		{2.1755510569e-01,	7.1706706285e-01},
		{-7.2723817825e-01,	5.9097236395e-01},
		{1.0715749264e+00,	-7.5527334213e-01},
		{-1.2566597462e+00,	7.6995909214e-01},
		{1.0549577475e+00,	-8.8587015867e-01},
		{1.8394696712e-01,	3.5138349980e-02},
		{1.3489142060e-01,	1.7980568409e+00},
		{4.3190231919e-01,	1.5910146236e+00},
		{9.9428802729e-01,	-7.1379506588e-01},
		{8.3084392548e-01,	2.8293383121e-01},
		{2.9168587923e-01,	1.3119093180e+00},
		{-1.5170594454e+00,	8.6988085508e-01},
		{-7.4240112305e-01,	-1.6298339367e+00},
		{-7.2394043207e-01,	-1.7632457018e+00},
		{-9.5290309191e-01,	-1.7981760502e+00},
		{-8.4711390734e-01,	-1.1068773270e+00},
		{1.1934005022e+00,	-5.0828325748e-01},
		{1.2606358528e+00,	-1.1707264185e+00},
		{7.6515805721e-01,	1.0476279259e+00},
		{7.6225411892e-01,	1.7606880665e+00},
		{3.3772689104e-01,	2.5667250156e-01},
		{-6.0842096806e-01,	1.6842957735e+00},
		{5.4336506128e-01,	4.2362573743e-01},
		{4.0705797076e-01,	5.1623529196e-01},
		{6.0142117739e-01,	9.6227097511e-01},
		{-3.2689183950e-02,	1.1448560953e+00},
		{7.8598010540e-01,	1.4907100201e+00},
		{5.6932246685e-01,	-1.1631521583e-01},
		{-1.5453243256e+00,	6.8581354618e-01},
		{9.1223746538e-02,	1.3116534948e+00},
		{7.9869186878e-01,	-3.8744705915e-01},
		{-8.3318686485e-01,	1.1196424961e+00},
		{1.8072992563e-02,	1.7401126623e+00},
		{1.0915973186e+00,	-3.8188034296e-01}
	};
	data_io::write_default2d("maps/sample2d.dat", sample, 5);

	mveqf::MVEQF<int, double> qf;
	qf.set_kernel_type(4);
	std::vector<size_t> empty;
	qf.set_sample_and_bounds(sample, std::vector<double>(2, -2), std::vector<double>(2, 2), empty, 500, 1e-8, 100000);

	std::mt19937_64 generator;
	generator.seed(1);
	std::uniform_real_distribution<double> ureal01(0.0,1.0);
	size_t nrolls = 500;
	std::vector<std::vector<double>> sampled(nrolls, std::vector<double>(2));
	for(size_t i = 0; i != nrolls; i++)
	{
		std::vector<double> t = {ureal01(generator), ureal01(generator)};
		sampled[i] = qf.transform(t);
	}
	data_io::write_default2d("maps/sampled2d.dat", sampled, 5);


	nrolls = 50000;
	sampled = std::vector<std::vector<double>>(nrolls, std::vector<double>(2));
	auto vals01 = sampled;
	for(auto & i : vals01)
	{
		for(auto & j : i)
		{
			j = ureal01(generator);
		}
	}

	const auto nthreads = std::thread::hardware_concurrency();
	auto first = sampled.begin();
	auto last = sampled.end();
	const auto size = last - first;
	const auto size_per_thread = size / nthreads;

	std::vector<std::future<void>> futures;
	for(unsigned int i = 0; i < nthreads - 1; i++)
	{
		futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, &vals01, &sampled, &qf]()
		{
			for(auto it = start; it != start + size_per_thread; ++it)
			{
				qf.transform(vals01[std::distance(sampled.begin(), it)], *it);
			}
		}));
	}
	futures.emplace_back(
	  std::async([start = first + (nthreads - 1) * size_per_thread, last, &vals01, &sampled, &qf]()
	{
		for(auto it = start; it != last; ++it)
		{
			qf.transform(vals01[std::distance(sampled.begin(), it)], *it);
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
	data_io::write_default2d("maps/sampled2d_1m.dat", sampled, 5);
}


void test_mveqf3()
{
	std::vector<std::vector<double>> sample =
	{
		{5.000000e+01,	5.000000e+01},
		{-5.000000e+01,	5.000000e+01},
		{5.000000e+01,	-5.000000e+01},
		{-5.000000e+01,	-5.000000e+01},
		{-3.603057e+01,	7.704624e+01},
		{-3.210936e+01,	-3.876422e+01},
		{-5.673265e+00,	-6.626629e+01},
		{-9.162707e+01,	1.064726e+01},
		{-5.508124e+01,	-9.601210e+01},
		{-9.279164e+01,	-2.018437e-01}
	};
	data_io::write_default2d("maps/sample2d.dat", sample, 5);

//    std::vector<size_t> count = {1, 25, 84, 35, 101, 36, 77, 54, 76, 62};
//    for(size_t i = 0; i != count.size(); i++)
//    {
//        for(size_t j = 0; j < count[i] - 1; j++)
//        {
//            sample.push_back(sample[i]);
//        }
//    }

	mveqf::MVEQF<int, double> qf;
	qf.set_kernel_type(0);
	std::vector<size_t> empty;// = {1, 25, 84, 35, 101, 36, 77, 54, 76, 62};
	qf.set_sample_and_bounds(sample, std::vector<double>(2, -100), std::vector<double>(2, 100), empty, 500, 1e-8, 100000);

	std::mt19937_64 generator;
	generator.seed(1);
	std::uniform_real_distribution<double> ureal01(0.0,1.0);
	size_t nrolls = 500;
	std::vector<std::vector<double>> sampled(nrolls, std::vector<double>(2));
	for(size_t i = 0; i != nrolls; i++)
	{
		std::vector<double> t = {ureal01(generator), ureal01(generator)};
		sampled[i] = qf.transform(t);
	}
	data_io::write_default2d("maps/sampled2d.dat", sampled, 5);


	nrolls = 50000;
	sampled = std::vector<std::vector<double>>(nrolls, std::vector<double>(2));
	auto vals01 = sampled;
	for(auto & i : vals01)
	{
		for(auto & j : i)
		{
			j = ureal01(generator);
		}
	}

	const auto nthreads = std::thread::hardware_concurrency();
	auto first = sampled.begin();
	auto last = sampled.end();
	const auto size = last - first;
	const auto size_per_thread = size / nthreads;

	std::vector<std::future<void>> futures;
	for(unsigned int i = 0; i < nthreads - 1; i++)
	{
		futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, &vals01, &sampled, &qf]()
		{
			for(auto it = start; it != start + size_per_thread; ++it)
			{
				qf.transform(vals01[std::distance(sampled.begin(), it)], *it);
			}
		}));
	}
	futures.emplace_back(
	  std::async([start = first + (nthreads - 1) * size_per_thread, last, &vals01, &sampled, &qf]()
	{
		for(auto it = start; it != last; ++it)
		{
			qf.transform(vals01[std::distance(sampled.begin(), it)], *it);
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
	data_io::write_default2d("maps/sampled2d_1m.dat", sampled, 5);
}


void test_mveqf4()
{

	std::vector<std::vector<double>> sample;

	std::mt19937_64 generator;
	generator.seed(1);
	std::uniform_real_distribution<double> ureal01(0.0,1.0);
	std::uniform_real_distribution<double> ureal11(-0.0,1.0);
	std::uniform_real_distribution<double> ureal360(0.0,360.0);
	std::lognormal_distribution<> lognormal(1.6, 0.5);
	std::normal_distribution<> norm(0, 1.0);

	for(size_t i = 0; i != 5000; i++)
	{
		std::vector<double> shift = {5, 5, 5};
		double tt = lognormal(generator);
		while(tt > 25)
			tt = lognormal(generator);
		std::vector<double> t = {-5.0 + tt, norm(generator), norm(generator)};
		for(size_t j = 0; j != t.size(); j++)
			t[j] += shift[j];
		double x = t[0];
		double y = t[1];
		double tau = 3.1415/4;//acos(-1.0)/2.0;
		t[0] = x*std::cos(tau) - y*std::sin(tau) + 11;
		t[1] = y*std::cos(tau) + x*std::sin(tau);


		x = t[0];
		y = t[2];
		tau = 3.1415/4;//acos(-1.0)/2.0;
		t[0] = x*std::cos(tau) - y*std::sin(tau);
		t[2] = y*std::cos(tau) + x*std::sin(tau);

		for(size_t j = 0; j != t.size(); j++)
		{
			if(t[j] < 0)
				std::cout << "!!!" << '\t' << j << std::endl;
		}

		sample.push_back(t);
	}
	data_io::write_default2d("maps/Trie/3d2/sample3d.dat", sample, 5);

	mveqf::MVEQF<int, double> qf;
	qf.set_kernel_type(0);
	std::vector<size_t> empty;
	qf.set_sample_and_bounds(sample, std::vector<double>(3, -10), std::vector<double>(3, 10), empty, 500, 1e-8, 100000);

	size_t nrolls = 500;
	std::vector<std::vector<double>> sampled(nrolls, std::vector<double>(3));
	for(size_t i = 0; i != nrolls; i++)
	{
		std::vector<double> t = {ureal01(generator), ureal01(generator)};
		sampled[i] = qf.transform(t);
	}
	data_io::write_default2d("maps/Trie/3d2/sampled3d.dat", sampled, 5);


	/*nrolls = 50000;
	sampled = std::vector<std::vector<double>>(nrolls, std::vector<double>(3));
	auto vals01 = sampled;
	for(auto & i : vals01)
	{
	    for(auto & j : i)
	    {
	        j = ureal01(generator);
	    }
	}

	const auto nthreads = std::thread::hardware_concurrency();
	auto first = sampled.begin();
	auto last = sampled.end();
	const auto size = last - first;
	const auto size_per_thread = size / nthreads;

	std::vector<std::future<void>> futures;
	for(unsigned int i = 0; i < nthreads - 1; i++)
	{
	    futures.emplace_back(std::async([start = first + i * size_per_thread, size_per_thread, &vals01, &sampled, &qf]()
	    {
	        for(auto it = start; it != start + size_per_thread; ++it)
	        {
	            qf.transform(vals01[std::distance(sampled.begin(), it)], *it);
	        }
	    }));
	}
	futures.emplace_back(
	    std::async([start = first + (nthreads - 1) * size_per_thread, last, &vals01, &sampled, &qf]()
	{
	    for(auto it = start; it != last; ++it)
	    {
	        qf.transform(vals01[std::distance(sampled.begin(), it)], *it);
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
	data_io::write_default2d("maps/sampled3d_1m.dat", sampled, 5);*/
}

void run_test()
{
	timer::Timer time;
//
//    std::cout << tt.front() == tt.back() << std::endl;

//    simple_empirical_1d();
//    simple1d_example();

//    test_1d1();
//    test_1d2();
//    test_1d3();
//    test_1d4();
//    test_1d5();


	test_2d1();
//	test_2d2();
//    test_3d1();
//	test_3d2();

	//testNd::test_Nd_discrete();

//    test_1d_func();
//    test_2d_func();

//    kquantile::test_1d_kquantile();
//    test_2d_kquantile<double>();
//    test_3d_kquantile<double>();
//  test_3d_kquantile2<double>();

//    test_1d_uniform_vs_nonuniform();
//    test_2d_uniform_vs_nonuniform();

//	test_2d_uniform_vs_nonuniform_trie();
	//test_2d_random_area_trie();

//	ot_1d<double>();

//    test_grid_10d();

	/// 2-dimensional test

//    std::vector<size_t> g = {13, 20};
//    //std::vector<size_t> g(N, 10);
//    std::vector<float> lb = {-3, -2};
//    std::vector<float> ub = {1, 4};
//    test_Nd(g, lb, ub, 100, 1e5);

	/// N-dimensional test

//    std::mt19937_64 generator;
//    generator.seed(1);
//
//    std::uniform_int_distribution<size_t> dis(2, 20);
//    size_t N = 10;
//    std::vector<size_t> g(N);
//    for(size_t i = 0; i != N; i++)
//    {
//        g[i] = dis(generator);
//        std::cout << g[i] << '\t';
//    }
//    std::cout << std::endl;
//    std::vector<float> lb(N, -10);
//    std::vector<float> ub(N, 10);
//    test_Nd(g, lb, ub, 50000, 1e2);


//    size_t N = 2;
//    std::vector<size_t> g(N);
//    for(size_t i = 0; i != N; i++)
//    {
//        g[i] = 10000;
//    }
//    std::cout << std::endl;
//    std::vector<float> lb(N, -10);
//    std::vector<float> ub(N, 10);
//    test_Nd(g, lb, ub, 10000, 1e5);

	/// grid and dim test
//    grid_test_Nd();
//	testNd::dim_test_Nd();


//    size_t N = 1;
//    std::vector<size_t> g(N);
//    for(size_t i = 0; i != N; i++)
//    {
//        g[i] = 100;
//    }
//    std::cout << std::endl;
//    std::vector<float> lb(N, -10);
//    std::vector<float> ub(N, 10);
//    test_Nd(g, lb, ub, 50, 1e5);

//    grid_test_2d();
//    sample_size_test(10);

//    worst_space_test_dim();
//    worst_space_test_grid();
//    worst_space_test_grid_2d();

//    worst_space_check();
	//std::vector<std::vector<double>> x(100,std::vector<double>(10));


	/// kde
//    mveqf::kde::test1d();
//    mveqf::kde::test1d_1();
//    mveqf::kde::test2d();
//    mveqf::kde::test2d_2();

//    test_mveqf4();

	///
	std::cout << time.elapsed_seconds() << std::endl;
}

// // trie based layer test
//    trie_based::TrieBased<trie_based::NodeCount<int>,int>  sample;
//    sample.set_dimension(3);
//    std::vector<std::vector<int>> in_sample;
//    in_sample =
//    {
//        {1,0,0}, /* baa */
//        {2,0,0}, /* caa */
//        {4,0,0}, /* eaa */
//        {0,2,0}, /* aca */
//        {4,4,0}, /* eea */
//        {4,3,0}, /* eda */
//        {3,3,0}, /* dda */
//
//        {0,0,1}, /* aab */
//
//        {3,0,2}, /* dac */
//        {0,3,2}, /* adc */
//
//        {0,3,3}, /* add */
//
//        {2,0,4}, /* cae */
//        {2,1,4}, /* cbe */
//        {2,2,4} /* cce */
//    };
//    for(const auto & i : in_sample)
//        sample.insert(i);
//    sample.fill_tree_count();
//    //sort();
//    auto rez = sample.get_layer_count();
//    std::cout << rez.size() << std::endl;
//    std::cout << sample.last_layer.size() << std::endl;
//
//    std::cout << std::endl;
//
//    for(const auto &i : rez)
//    {
//        std::cout << i.first << std::endl;
//        for(const auto &j : i.second)
//            std::cout << j << ' ';
//        std::cout << std::endl;
//    }


#endif
