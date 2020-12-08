#include <iostream>
#include <vector>
#include <random>
#include <mveqf/sdot.h>

int main()
{
	const size_t n = 10; // number of init dots
	const size_t dimension = 2;
	std::vector<std::vector<float>> init_dots(n, std::vector<float>(dimension));
	std::vector<size_t> weights(n);

	std::mt19937_64 generator;
	generator.seed(1);
	std::uniform_real_distribution<float> ureal(-100.0, 100.0);
	std::uniform_real_distribution<float> ureal01(0.0,1.0);
	for(size_t i = 0; i != n; i++)
	{
		for(size_t j = 0; j != dimension; j++)
		{
			init_dots[i][j] = ureal(generator);
			weights[i] = std::abs(ureal(generator))*100*ureal01(generator) + 1;
		}
	}

	mveqf::ot::SDOT<float> obj(init_dots, weights);
	std::vector<size_t> counter(n, 0);
	std::vector<size_t> times(n);
	std::vector<std::vector<size_t>> unique_values;
	//	std::vector<std::vector<TFloat>> sampled;
	for(size_t i = 0; i != 5000; i++)
	{
		std::vector<float> temp(dimension);
		for(size_t j = 0; j != dimension; j++)
		{
			temp[j] = ureal01(generator);
		}


		++counter[obj.transform_index(temp)];
		//std::cout << obj.transform_index(temp) << std::endl;
		std::vector<size_t> t = obj.transform(temp);//obj.transform_float(temp);
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
	std::cout << "counter\t\ttimes" << std::endl;
	for(size_t i = 0; i != counter.size(); i++)
		std::cout << counter[i] << "\t\t" << times[i] << std::endl;
	std::cout << std::endl;
	std::cout << unique_values.size() << std::endl;
}
