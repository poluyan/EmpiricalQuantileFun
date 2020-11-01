#include <iostream>
#include <vector>
#include <random>
#include <mveqf/implicit.h>

int main()
{
	size_t dimension = 3;

	std::vector<size_t> grid = {5, 5, 5}; // grid

	// sample grid points
	std::vector<std::vector<int>> sample =
	{
		{1,0,0}, {2,0,0},	{4,0,0}, {0,2,0},	{4,4,0}, {4,3,0},	{3,3,0},
		{0,0,1}, {3,0,2}, {0,3,2}, {0,3,3}, {2,0,4}, {2,1,4}, {2,2,4}
	};

	std::vector<float> lb(dimension, -3.0f); // lower bound
	std::vector<float> ub(dimension, 3.0f);  // upper bound

	mveqf::ImplicitQuantile<int, float> mveqfunc(lb, ub, grid);
	mveqfunc.set_sample(sample);

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
