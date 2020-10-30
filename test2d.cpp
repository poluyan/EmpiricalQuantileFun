#include <iostream>
#include <vector>
#include <random>
#include <mveqf/implicit.h>

int main()
{
	size_t dimension = 2;

	std::vector<size_t> grid = {9, 10}; // grid

	// sample grid points
	std::vector<std::vector<int>> sample_implicit =
	{
		{2,6}, {3,2}, {3,3}, {3,5}, {3,6}, {3,7}, {4,5}, {4,6},
		{4,7}, {5,3}, {5,4}, {5,5}, {5,6}, {5,7},	{6,3}, {6,4}
	};

	std::vector<float> lb(dimension, -3.0f); // lower bound
	std::vector<float> ub(dimension, 3.0f);  // upper bound

	mveqf::ImplicitQuantile<int, float> mveqfunc(lb, ub, grid);
	mveqfunc.set_sample(sample_implicit);

	std::mt19937_64 generator;
	generator.seed(1);
	std::uniform_real_distribution<float> ureal01(0.0f, 1.0f);

	std::vector<float> values01(dimension);
	std::vector<float> sampled(dimension);

	for(size_t i = 0; i != 1000; i++)
	{
		for(auto & j : values01)
			j = ureal01(generator);
		
		mveqfunc.transform(values01, sampled);
		
		for(const auto & j : sampled)
			std::cout << std::fixed << j << '\t';
		std::cout << std::endl;
	}
}
