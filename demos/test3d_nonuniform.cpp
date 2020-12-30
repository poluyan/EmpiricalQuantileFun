#include <iostream>
#include <vector>
#include <random>
#include <mveqf/implicit_trie.h>

float threeExp(float x, float y, float z)
{
	std::vector<std::vector<float>> centers = { {3, 3, 3}, {-5, 0, 5}, {0, -5, -4} };
	float rez = 0;
	rez += std::exp(-(std::pow(x - centers[0][0], 2.0) + std::pow(y - centers[0][1], 2.0) + std::pow(z - centers[0][2], 2.0))*0.75);
	rez += std::exp(-(std::pow(x - centers[1][0], 2.0) + std::pow(y - centers[1][1], 2.0) + std::pow(z - centers[1][2], 2.0))*0.5)*0.75;
	rez += std::exp(-(std::pow(x - centers[2][0], 2.0) + std::pow(y - centers[2][1], 2.0) + std::pow(z - centers[2][2], 2.0))*0.25)*0.5;
	return rez;
}

int main()
{
	size_t dimension = 3;

	std::vector<size_t> grid = {50, 50, 50}; // grid

	// sample grid points
	auto sample = std::make_shared<mveqf::Trie<mveqf::NodeCount<int>,int>>();
	sample->set_dimension(dimension);

	std::vector<float> dx(dimension);
	std::vector<std::vector<float>> grids(dimension);
	for(size_t i = 0; i != grids.size(); i++)
	{
		std::vector<float> g(grid[i] + 1);
		float startp = -10.0;
		float endp = 10.0;
		float es = endp - startp;
		for(size_t j = 0; j != g.size(); j++)
		{
			g[j] = startp + j*es/static_cast<float>(grid[i]);
		}
		grids[i] = g;
		dx[i] = es/(static_cast<float>(grid[i])*2);
	}

	for(size_t i = 0; i != grid[0]; i++)
	{
		for(size_t j = 0; j != grid[1]; j++)
		{
			for(size_t k = 0; k != grid[2]; k++)
			{
				size_t value = static_cast<size_t>(1000.0f*threeExp(grids[0][i] + dx[0], grids[1][j] + dx[1], grids[2][k] + dx[2]));
				if(value)
				{
					std::vector<int> point = {int(i), int(j), int(k)};
					sample->insert(point, value);
				}
			}
		}
	}

	std::vector<float> lb(dimension, -10.0f); // lower bound
	std::vector<float> ub(dimension, 10.0f);  // upper bound

	mveqf::ImplicitTrieQuantile<int, float> mveqfunc(lb, ub, grid);
	mveqfunc.set_sample_shared(sample);

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
