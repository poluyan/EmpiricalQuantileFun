Multivariate empirical continuous quantile function (grid-based). There are two approaches to quantile function evaluation depending on the type of sample storage. In the first case, the sample is presented in the explicit (real-valued) form and stored in the matrix. In the second case, the sample is presented in the implicit form and the trie-based structure. Here presented header-only library that allows you to perform quantile transforms based on given sample points. For more info and examples see: [poluyan.github.io/mveqf](https://poluyan.github.io/mveqf/)

### Requirements

To compile from source, you need C++ 17 compiler and CMake for building examples.

### Installing

To use this library and perform quantile tranforms only header files from `mveqf` are needed. 

### Examples

Some examples of using `mveqf` to perform quantile transform presented in `demos` directory. Follow these steps to build and run the examples. After these steps all the binaries should be generated and presented in the `bin` directory.

##### Linux (gcc/clang)

```sh
$ git clone https://github.com/poluyan/mveqf
$ cd mveqf
$ cmake .
$ make
```

##### Windows (Visual Studio 2019+ with MSVC)

Clone the entire repository and build it locally. 

### Cite

S. V. Poluyan, N. M. Ershov, Quantile transform in structural bioinformatics problems // Computational nanotechnology, 2019, Vol. 6, no. 4, P. 29–43 DOI: [10.33693/2313-223X-2019-6-4-29-43](http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=cn&paperid=274&option_lang=eng)

```

```

### License

The `mveqf` library is distributed under Apache License 2.0 and it is open-source software. Feel free to make a copy and modify the source code, but keep the copyright notice and license intact.
