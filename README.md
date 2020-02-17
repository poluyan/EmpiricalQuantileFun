MVEQF
==========
Multivariate empirical continuous quantile function (grid-based). There are two approaches to quantile function evaluation depending on the type of sample storage. In the first case, the sample is presented in the explicit (real-valued) form and stored in the matrix. In the second case, the sample is presented in the implicit form and the trie-based structure. Here presented *header-based* library that allows you to perform quantile transforms based on given sample points. For more info and examples see: [poluyan.github.io/MVEQF](https://poluyan.github.io/MVEQF/)
### Installing
```sh
$ git clone https://github.com/poluyan/MVEQF
```
To perform quantile tranforms only header files from `src` are needed.
### Usage
Here is an example of two sample points and generated kernel density estimation with selected kernel type.
```cpp
#include <random>
#include <iostream>
#include <vector>

#include <mveqf.h> // the only header from MVEQF 

int main()
{
	size_t kernel_type = 0; // 0 - Gaussian, 1 - Epanechnikov, 2 - Uniform, 3 - Biweight,  4 - Triweight, 5 - Laplacian

	// 2 points example
	std::vector<std::vector<double>> sample = { {-1.13, -1.08}, {1.11, 1.27} };

	std::vector<double> lower_bound = {-3, -3};
	std::vector<double> upper_bound = {3, 3};

	mveqf::MVEQF<int, double> qf;
	qf.set_kernel_type(kernel_type); // 0 is default
	qf.set_sample_and_bounds(sample, lower_bound, upper_bound);

	// for generating points from [0,1]^n
	std::mt19937_64 generator;
	generator.seed(1);
	std::uniform_real_distribution<double> ureal01(0,1);

	for(size_t i = 0; i != 1000; i++) 
	{
		std::vector<double> t = {ureal01(generator), ureal01(generator)};
		auto sampled = qf.transform(t);
		std::cout << std::scientific << sampled.front() << '\t' << sampled.back() << std::endl;
	}
}
```
To build this example specify a directory to search for include files and thread support library, e.g. for `g++`
```
$ g++ -std=c++17 -IMVEQF/src main.cc -pthread
```
The figure shows the result of the example above for different number of sampled points. The initial sample points, initial density estimation with Gaussian kernel, and 10^3, 10^5, 10^7 sampled points with corresponding kernel density estimation are presented.
![Alt text](./maps/2dpdfexample1.png)
```cpp
// 5 points example
std::vector<std::vector<double>> sample = { {-1.13, -1.08}, {1.11, 1.27}, {0.86, 0.74}, {0.97, 1.31}, {0.79, -1.15} };
```
![Alt text](./maps/2dpdfexample2.png)
### License
The MVEQF is distributed under Apache License 2.0 and it is open-source software. Feel free to make a copy and modify the source code, but keep the copyright notice and license intact.
