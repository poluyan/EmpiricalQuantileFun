mveqf
==========
Multivariate empirical continuous quantile function (grid-based). There are two approaches to quantile function evaluation depending on the type of sample storage. In the first case, the sample is presented in the explicit (real-valued) form and stored in the matrix. In the second case, the sample is presented in the implicit form and the trie-based structure. Here presented header-only library that allows you to perform quantile transforms based on given sample points. For more info and examples see: [poluyan.github.io/MVEQF](https://poluyan.github.io/MVEQF/)

### Installing
```sh
$ git clone https://github.com/poluyan/MVEQF
$ cmake .
$ make
```
To perform quantile tranforms only header files from `mveqf` are needed.

### License
The MVEQF is distributed under Apache License 2.0 and it is open-source software. Feel free to make a copy and modify the source code, but keep the copyright notice and license intact.
