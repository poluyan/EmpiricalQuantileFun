#include <iostream>
#include <vector>
#include <random>
#include <mveqf/implicit.h>

int main()
{
	using float_type = float; // float, double
	using grid_node_type = std::uint8_t; // char, unsigned char, int, std::uint32_t, ...

	size_t dimension = 1;

	std::vector<size_t> grid = {12}; // grid

	// sample grid points, non-uniform distribution
	typedef mveqf::TrieBased<mveqf::NodeCount<grid_node_type>,grid_node_type> trie_type;
	std::shared_ptr<trie_type> sample = std::make_shared<trie_type>();
	sample->set_dimension(dimension);
	sample->insert(std::vector<grid_node_type> {1});
	sample->insert(std::vector<grid_node_type> {2});
	sample->insert(std::vector<grid_node_type> {3});
	sample->insert(std::vector<grid_node_type> {7});
	sample->insert(std::vector<grid_node_type> {8});
	sample->insert(std::vector<grid_node_type> {9});
	sample->insert(std::vector<grid_node_type> {10});

	std::vector<float_type> lb(dimension, -2.0f); // lower bound
	std::vector<float_type> ub(dimension, 4.0f);  // upper bound

	mveqf::ImplicitQuantile<grid_node_type, float_type> mveqfunc(lb, ub, grid);
	mveqfunc.set_sample_shared_and_fill_count(sample);

	std::mt19937_64 generator;
	generator.seed(1);
	std::uniform_real_distribution<float_type> ureal01(0.0f, 1.0f);

	std::vector<float_type> values01(dimension);
	std::vector<float_type> sampled(dimension);

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
