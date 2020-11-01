#include <iostream>
#include <vector>
#include <random>
#include <mveqf/implicit.h>

int main()
{
	size_t dimension = 1;

	std::vector<size_t> grid = {12}; // grid

	// sample grid points, non-uniform distribution
	typedef mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<int>,int> trie_type;
	std::shared_ptr<trie_type> sample = std::make_shared<trie_type>();
	sample->set_dimension(dimension);
	sample->insert(std::vector<int> {1});
	sample->insert(std::vector<int> {2});
	sample->insert(std::vector<int> {3});
	sample->insert(std::vector<int> {7});
	sample->insert(std::vector<int> {8});
	sample->insert(std::vector<int> {9});
	sample->insert(std::vector<int> {10});

	std::vector<float> lb(dimension, -2.0f); // lower bound
	std::vector<float> ub(dimension, 4.0f);  // upper bound

	mveqf::ImplicitQuantile<int, float> mveqfunc(lb, ub, grid);
	mveqfunc.set_sample_shared_and_fill_count(sample);

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
