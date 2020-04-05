/**************************************************************************

   Copyright © 2018 Sergey Poluyan <svpoluyan@gmail.com>

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
#include <test/test2d.h>
#include <utility/data_io.h>
#include <test/test.h>
#include <utility/timer.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>

void test_2d1()
{
  std::vector<size_t> grid_number = {9, 10};

  std::vector<std::vector<int>> sample_implicit =
  {
    {2,6},

    {3,2},
    {3,3},
    {3,5},
    {3,6},
    {3,7},

    {4,5},
    {4,6},
    {4,7},

    {5,3},
    {5,4},
    {5,5},
    {5,6},
    {5,7},

    {6,3},
    {6,4}
  };

/// multivariate quantile function [0,1]^n -> [-3,3]^n nrolls = 2e3
  explicit_quantile(-3, 3, grid_number, sample_implicit, 2e3);
  implicit_quantile_class(-3, 3, grid_number, sample_implicit, 2e3);
  implicit_quantile_class_sorted(-3, 3, grid_number, sample_implicit, 2e3);
}

void test_2d2()
{
  std::vector<size_t> grid_number;
  std::vector<std::vector<int>> sample_implicit;
  data_io::load_grid_and_sample("input/2d/grid.dat", "input/2d/points.dat", grid_number, sample_implicit);

  /// multivariate quantile function [0,1]^n -> [-3,3]^n
  timer::Timer time_cpp11;
  time_cpp11.reset();
  explicit_quantile(-3, 3, grid_number, sample_implicit, 2e+1);
  std::cout << "--------->   total time: " << time_cpp11.elapsed_seconds() << std::endl;
  time_cpp11.reset();
  implicit_quantile_class(-3, 3, grid_number, sample_implicit, 2e+1);
  std::cout << "--------->   total time: " << time_cpp11.elapsed_seconds() << std::endl;
  time_cpp11.reset();
  implicit_quantile_class_sorted(-3, 3, grid_number, sample_implicit, 2e+1);
  std::cout << "--------->   total time: " << time_cpp11.elapsed_seconds() << std::endl;
}




std::vector<double> test_2d_time(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nsamples, size_t Nrolls)
{
  std::mt19937_64 generator;
  generator.seed(1 + gridN.size() + gridN.front());
  std::uniform_real_distribution<float> ureal01(0.0,1.0);

  typedef mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<int>,int> sample_type;
  std::shared_ptr<sample_type> sample = std::make_shared<sample_type>();

  std::vector<std::vector<int>> sample_int;
  std::vector<int> temp(gridN.size());
  for(size_t i = 0; i != Nsamples;)
  {
    for(size_t j = 0; j != gridN.size(); j++)
    {
      temp[j] = static_cast<int>(std::round(ureal01(generator)*(gridN[j] - 1.0)));
    }
    if(!sample->search(temp))
    {
      sample->insert(temp);
      sample_int.push_back(temp);
      i++;
    }
  }

  std::vector<std::vector<float> > values01;

  std::vector<float> temp1(gridN.size());
  std::vector<float> temp2(temp1.size());

  for(size_t i = 0; i != Nrolls; ++i)
  {
    for(size_t j = 0; j != temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    values01.push_back(temp1);
  }

  for(size_t i = 0; i != temp1.size(); ++i)
  {
    for(size_t j = 0; j != temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    temp1[i] = 0.0;
    values01.push_back(temp1);
  }
  for(size_t i = 0; i != temp1.size(); ++i)
  {
    for(size_t j = 0; j != temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    temp1[i] = 1.0;
    values01.push_back(temp1);
  }

  for(size_t i = 0; i != temp1.size(); ++i)
  {
    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
    {
      temp1[j] = 0.0;
    }
    for(size_t j = i + 1; j < temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    values01.push_back(temp1);
  }
  for(size_t i = 0; i != temp1.size(); ++i)
  {
    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
    {
      temp1[j] = 1.0;
    }
    for(size_t j = i + 1; j < temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    values01.push_back(temp1);
  }

  for(size_t i = 0; i != temp1.size() - 1; ++i)
  {
    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    for(size_t j = i + 1; j < temp1.size(); j++)
    {
      temp1[j] = 0.0;
    }
    values01.push_back(temp1);
  }
  for(size_t i = 0; i != temp1.size() - 1; ++i)
  {
    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    for(size_t j = i + 1; j < temp1.size(); j++)
    {
      temp1[j] = 1.0;
    }
    values01.push_back(temp1);
  }

  std::vector<std::vector<float> > sampled(values01.size(), std::vector<float>(values01.front().size()));

  timer::Timer time_all_trans, full_time;
  std::vector<double> result;

//    empirical_quantile::ExplicitQuantile<int, float> quant_expl(lb, ub, gridN);
//    quant_expl.set_sample(sample_int);
//    full_time.reset();
//    time_all_trans.reset();
//    for(size_t i = 0; i != values01.size(); i++)
//        quant_expl.transform(values01[i], sampled[i]);
//    result.push_back(time_all_trans.elapsed_seconds());
//    result.push_back(result.back()/double(sampled.size()));
//    result.push_back(full_time.elapsed_seconds());
//    for(size_t i = 0; i != sampled.size(); i++)
//    {
//        for(size_t j = 0; j != sampled[i].size(); j++)
//        {
//            if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
//            {
//                std::cout << "beyond bounds" << std::endl;
//                std::cout << sampled[i][j] << std::endl;
//            }
//        }
//    }

  mveqf::ImplicitQuantile<int, float> quant_impl(lb, ub, gridN);
  full_time.reset();
  quant_impl.set_sample_shared_and_fill_count(sample);
  time_all_trans.reset();
  for(size_t i = 0; i != values01.size(); i++)
    quant_impl.transform(values01[i], sampled[i]);
  result.push_back(time_all_trans.elapsed_seconds());
  result.push_back(result.back()/double(sampled.size()));
  result.push_back(full_time.elapsed_seconds());
  for(size_t i = 0; i != sampled.size(); i++)
  {
    for(size_t j = 0; j != sampled[i].size(); j++)
    {
      if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
      {
        std::cout << "beyond bounds" << std::endl;
        std::cout << sampled[i][j] << std::endl;
      }
    }
  }

  mveqf::ImplicitQuantileSorted<int, float> quant_impls(lb, ub, gridN);
  full_time.reset();
  quant_impls.set_sample_shared_and_fill_count(sample);
  time_all_trans.reset();
  for(size_t i = 0; i != values01.size(); i++)
    quant_impls.transform(values01[i], sampled[i]);
  result.push_back(time_all_trans.elapsed_seconds());
  result.push_back(result.back()/double(sampled.size()));
  result.push_back(full_time.elapsed_seconds());
  for(size_t i = 0; i != sampled.size(); i++)
  {
    for(size_t j = 0; j != sampled[i].size(); j++)
    {
      if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
      {
        std::cout << "beyond bounds" << std::endl;
        std::cout << sampled[i][j] << std::endl;
      }
    }
  }

  mveqf::ImplicitQuantileSortedInterp<int, float> quant_implsi(lb, ub, gridN);
  full_time.reset();
  quant_implsi.set_sample_shared_and_fill_count(sample);
  time_all_trans.reset();
  for(size_t i = 0; i != values01.size(); i++)
    quant_implsi.transform(values01[i], sampled[i]);
  result.push_back(time_all_trans.elapsed_seconds());
  result.push_back(result.back()/double(sampled.size()));
  result.push_back(full_time.elapsed_seconds());
  for(size_t i = 0; i != sampled.size(); i++)
  {
    for(size_t j = 0; j != sampled[i].size(); j++)
    {
      if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
      {
        std::cout << "beyond bounds" << std::endl;
        std::cout << sampled[i][j] << std::endl;
      }
    }
  }
  return result;
}

void grid_test_2d()
{
  size_t N = 2;
  std::vector<size_t> g(N, 1440);
  std::vector<float> lb(N, -acos(-1.0));
  std::vector<float> ub(N, acos(-1.0));
  auto rez = test_2d_time(g, lb, ub, 1e6, 1e4);
  for(const auto &i : rez)
    std::cout << std::scientific << i << '\t';
  std::cout << std::endl;
}

double threeExp(double x, double y)
{
  std::vector<std::vector<double>> centers = { {3, 3}, {-5, 0}, {0, -5} };
  double rez = 0;
  rez += std::exp(-(std::pow(x - centers[0][0], 2.0) + std::pow(y - centers[0][1], 2.0))*0.75);
  rez += std::exp(-(std::pow(x - centers[1][0], 2.0) + std::pow(y - centers[1][1], 2.0))*0.5)*0.75;
  rez += std::exp(-(std::pow(x - centers[2][0], 2.0) + std::pow(y - centers[2][1], 2.0))*0.25)*0.5;
  return rez;
}

void test_2d_func()
{
  std::vector<size_t> gridn = {100, 100};
//    std::vector<size_t> gridn = {50, 50};

  auto sample = std::make_shared<mveqf::trie_based::Trie<mveqf::trie_based::NodeCount<int>,int>>();
//    auto sample = std::make_shared<trie_based::TrieBased<trie_based::NodeCount<int>,int>>();
  sample->set_dimension(gridn.size());
//    std::vector<std::vector<int>> sample;

  std::vector<double> dx(gridn.size());
  std::vector<std::vector<double>> grids(gridn.size());
  for(size_t i = 0; i != grids.size(); i++)
  {
    std::vector<double> grid(gridn[i] + 1);
    double startp = -10.0;
    double endp = 10.0;
    double es = endp - startp;
    for(size_t j = 0; j != grid.size(); j++)
    {
      grid[j] = startp + j*es/double(gridn[i]);
    }
    grids[i] = grid;
    dx[i] = es/(double(gridn[i])*2);
  }

  for(size_t i = 0; i != gridn[0]; i++)
  {
    for(size_t j = 0; j != gridn[1]; j++)
    {
      //std::cout << grids[0][i] + dx[0] << '\t' << grids[1][j] + dx[1] << std::endl;
      size_t v = size_t(1000.0*threeExp(grids[0][i] + dx[0], grids[1][j] + dx[1]));
      if(v)
      {
        std::vector<int> a = {int(i), int(j)};
        //std::cout << "add " << i << '\t' << j << std::endl;

//                sample->insert(a);
        sample->insert(a, v);


//                sample->insert(a, v);
//                for(size_t k = 0; k != v; k++)
//                    sample.push_back(a);
      }
    }
  }

  /*std::vector<int> a = {5,5};
  //sample->insert(a, 4);
  sample.push_back(a);
  sample.push_back(a);
  sample.push_back(a);
  sample.push_back(a);
  sample.push_back(a);

  a[0] = 4;
  a[1] = 5;
  //sample->insert(a, 2);
  sample.push_back(a);
  sample.push_back(a);
  sample.push_back(a);

  a[0] = 5;
  a[1] = 4;
  //sample->insert(a, 2);
  sample.push_back(a);*/


  size_t nrolls = 1e+5;

  std::mt19937_64 generator;
  generator.seed(1);
  std::uniform_real_distribution<float> ureal01(0.0,1.0);

  std::vector<std::vector<float> > sampled;
  std::vector<std::vector<float> > values01;

  float lb = -10.0, ub = 10.0;

//    empirical_quantile::ImplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
//    quant.set_sample_shared_and_fill_count(sample);

  mveqf::ImplicitTrieQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
  quant.set_sample_shared(sample);

//    empirical_quantile::ExplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
//    quant.set_sample(sample);

  timer::Timer time_cpp11;
  time_cpp11.reset();

  std::vector<float> temp1(gridn.size());
  std::vector<float> temp2(temp1.size());

  for(size_t i = 0; i != nrolls; ++i)
  {
    for(size_t j = 0; j != temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }

  for(size_t i = 0; i != temp1.size(); ++i)
  {
    for(size_t j = 0; j != temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    temp1[i] = 0.0;
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }
  for(size_t i = 0; i != temp1.size(); ++i)
  {
    for(size_t j = 0; j != temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    temp1[i] = 1.0;
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }

  for(size_t i = 0; i != temp1.size(); ++i)
  {
    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
    {
      temp1[j] = 0.0;
    }
    for(size_t j = i + 1; j < temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }
  for(size_t i = 0; i != temp1.size(); ++i)
  {
    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
    {
      temp1[j] = 1.0;
    }
    for(size_t j = i + 1; j < temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }

  for(size_t i = 0; i != temp1.size() - 1; ++i)
  {
    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    for(size_t j = i + 1; j < temp1.size(); j++)
    {
      temp1[j] = 0.0;
    }
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }
  for(size_t i = 0; i != temp1.size() - 1; ++i)
  {
    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    for(size_t j = i + 1; j < temp1.size(); j++)
    {
      temp1[j] = 1.0;
    }
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }
  data_io::write_default2d("maps/sampled_implicit.dat", sampled, 5);
}


void test_2d_uniform_vs_nonuniform()
{
  std::vector<size_t> gridn = {12, 12};

//    auto sample = std::make_shared<trie_based::TrieBased<trie_based::NodeCount<int>,int>>();
//    sample->set_dimension(gridn.size());
//    sample->insert(std::vector<int>{5,5});
//    sample->insert(std::vector<int>{4,5});
//    sample->insert(std::vector<int>{5,4});
//    sample->insert(std::vector<int>{1,3});



  auto sample = std::make_shared<mveqf::trie_based::Trie<mveqf::trie_based::NodeCount<int>,int>>();
  sample->set_dimension(gridn.size());
  sample->insert(std::vector<int> {5,5}, 1);
  sample->insert(std::vector<int> {4,5}, 1);
  sample->insert(std::vector<int> {5,4}, 1);
  sample->insert(std::vector<int> {1,3}, 1);


  size_t nrolls = 1e+4;

  std::mt19937_64 generator;
  generator.seed(1);
  std::uniform_real_distribution<float> ureal01(0.0,1.0);

  std::vector<std::vector<float> > sampled;
  std::vector<std::vector<float> > values01;

  float lb = -2.0, ub = 2.0;

//    empirical_quantile::ImplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
//    quant.set_sample_shared_and_fill_count(sample);

  mveqf::ImplicitTrieQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
  quant.set_sample_shared(sample);

//    empirical_quantile::ExplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
//    quant.set_sample(sample);

  timer::Timer time_cpp11;
  time_cpp11.reset();

  std::vector<float> temp1(gridn.size());
  std::vector<float> temp2(temp1.size());

  for(size_t i = 0; i != nrolls; ++i)
  {
    for(size_t j = 0; j != temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }

  for(size_t i = 0; i != temp1.size(); ++i)
  {
    for(size_t j = 0; j != temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    temp1[i] = 0.0;
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }
  for(size_t i = 0; i != temp1.size(); ++i)
  {
    for(size_t j = 0; j != temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    temp1[i] = 1.0;
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }

  for(size_t i = 0; i != temp1.size(); ++i)
  {
    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
    {
      temp1[j] = 0.0;
    }
    for(size_t j = i + 1; j < temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }
  for(size_t i = 0; i != temp1.size(); ++i)
  {
    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
    {
      temp1[j] = 1.0;
    }
    for(size_t j = i + 1; j < temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }

  for(size_t i = 0; i != temp1.size() - 1; ++i)
  {
    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    for(size_t j = i + 1; j < temp1.size(); j++)
    {
      temp1[j] = 0.0;
    }
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }
  for(size_t i = 0; i != temp1.size() - 1; ++i)
  {
    for(size_t j = 0; j != i + 1 && j < temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    for(size_t j = i + 1; j < temp1.size(); j++)
    {
      temp1[j] = 1.0;
    }
    quant.transform(temp1,temp2);
    values01.push_back(temp1);
    sampled.push_back(temp2);
  }
  data_io::write_default2d("maps/sampled_implicit.dat", sampled, 5);
}

int orientation(std::vector<double> p, std::vector<double> q, std::vector<double> r)
{
  double val = (q[1] - p[1]) * (r[0] - q[0]) -
               (q[0] - p[0]) * (r[1] - q[1]);
  if(val == 0) return 0;
  return (val > 0)? 1: 2;
}

void convexHull(const std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &hull)
{
  int n = points.size();
  if(n < 3) return;

  int l = 0;
  for(int i = 1; i < n; i++)
    if(points[i][0] < points[l][0])
      l = i;

  int p = l, q;
  do
  {
    std::vector<double> temp = {points[p][0], points[p][1]};
    hull.push_back(temp);

    q = (p+1)%n;
    for(int i = 0; i < n; i++)
    {
      if(orientation(points[p], points[i], points[q]) == 2)
        q = i;
    }
    p = q;
  }
  while(p != l);
  hull.push_back(hull.front());
}

void test_2d_uniform_vs_nonuniform_trie()
{
  std::vector<size_t> gridn = {9, 10};

//  auto sample = std::make_shared<mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<int>,int>>();
//  sample->set_dimension(gridn.size());
//  sample->insert(std::vector<int> {2,6});
//
//  sample->insert(std::vector<int> {3,2});
//  sample->insert(std::vector<int> {3,3});
//  sample->insert(std::vector<int> {3,5});
//  sample->insert(std::vector<int> {3,6});
//  sample->insert(std::vector<int> {3,7});
//
//  sample->insert(std::vector<int> {4,5});
//  sample->insert(std::vector<int> {4,6});
//  sample->insert(std::vector<int> {4,7});
//
//  sample->insert(std::vector<int> {5,3});
//  sample->insert(std::vector<int> {5,4});
//  sample->insert(std::vector<int> {5,5});
//  sample->insert(std::vector<int> {5,6});
//  sample->insert(std::vector<int> {5,7});
//
//  sample->insert(std::vector<int> {6,3});
//  sample->insert(std::vector<int> {6,4});



//  auto sample = std::make_shared<mveqf::trie_based::Trie<mveqf::trie_based::NodeCount<int>,int>>();
//  sample->set_dimension(gridn.size());
//  sample->insert(std::vector<int> {2,6}, 1);
//
//  sample->insert(std::vector<int> {3,2}, 1);
//  sample->insert(std::vector<int> {3,3}, 1);
//  sample->insert(std::vector<int> {3,5}, 1);
//  sample->insert(std::vector<int> {3,6}, 1);
//  sample->insert(std::vector<int> {3,7}, 1);
//
//  sample->insert(std::vector<int> {4,5}, 1);
//  sample->insert(std::vector<int> {4,6}, 9);
//  sample->insert(std::vector<int> {4,7}, 1);
//
//  sample->insert(std::vector<int> {5,3}, 1);
//  sample->insert(std::vector<int> {5,4}, 1);
//  sample->insert(std::vector<int> {5,5}, 1);
//  sample->insert(std::vector<int> {5,6}, 1);
//  sample->insert(std::vector<int> {5,7}, 1);
//
//  sample->insert(std::vector<int> {6,3}, 1);
//  sample->insert(std::vector<int> {6,4}, 1);

  auto sample = std::make_shared<mveqf::trie_based::Trie<mveqf::trie_based::NodeCount<int>,int>>();
  sample->set_dimension(gridn.size());
  sample->insert(std::vector<int> {2,6}, 342);

  sample->insert(std::vector<int> {3,2}, 242);
  sample->insert(std::vector<int> {3,3}, 776);
  sample->insert(std::vector<int> {3,5}, 342);
  sample->insert(std::vector<int> {3,6}, 293);
  sample->insert(std::vector<int> {3,7}, 650);

  sample->insert(std::vector<int> {4,5}, 195);
  sample->insert(std::vector<int> {4,6}, 327);
  sample->insert(std::vector<int> {4,7}, 391);

  sample->insert(std::vector<int> {5,3}, 194);
  sample->insert(std::vector<int> {5,4}, 399);
  sample->insert(std::vector<int> {5,5}, 256);
  sample->insert(std::vector<int> {5,6}, 195);
  sample->insert(std::vector<int> {5,7}, 216);

  sample->insert(std::vector<int> {6,3}, 306);
  sample->insert(std::vector<int> {6,4}, 269);

  std::vector<std::vector<double>> real_sample =
  {
    {-6.6666662693e-01,	-1.5000000000e+00},
    {-6.6666662693e-01,	-8.9999991655e-01},
    {8.9406967163e-08,	9.0000015497e-01},
    {8.9406967163e-08,	3.0000001192e-01},
    {-6.6666662693e-01,	3.0000001192e-01},
    {-6.6666662693e-01,	9.0000015497e-01},
    {-1.3333332539e+00,	9.0000015497e-01},
    {-6.6666662693e-01,	1.5000002384e+00},
    {8.9406967163e-08,	1.5000002384e+00},
    {6.6666686535e-01,	1.5000002384e+00},
    {6.6666686535e-01,	9.0000015497e-01},
    {6.6666686535e-01,	3.0000001192e-01},
    {6.6666686535e-01,	-2.9999989271e-01},
    {6.6666686535e-01,	-8.9999991655e-01},
    {1.3333333731e+00,	-8.9999991655e-01},
    {1.3333333731e+00,	-2.9999989271e-01}
  };

  std::vector<double> number(16,0);
  std::vector<std::vector<std::vector<double>>> values_n(16);
  std::vector<std::vector<std::vector<double>>> sampled_n(16);

  std::mt19937_64 generator;
  generator.seed(1);
  std::uniform_real_distribution<double> ureal01(0.0,1.0);

  std::vector<std::vector<double> > sampled;
  std::vector<std::vector<double> > values01;

  double lb = -3.0, ub = 3.0;

//  mveqf::ImplicitQuantile<int, double> quant(std::vector<double>(gridn.size(), lb), std::vector<double>(gridn.size(), ub), gridn);
//  quant.set_sample_shared_and_fill_count(sample);

  mveqf::ImplicitTrieQuantile<int, double> quant(std::vector<double>(gridn.size(), lb), std::vector<double>(gridn.size(), ub), gridn);
  quant.set_sample_shared(sample);

//    empirical_quantile::ExplicitQuantile<int, float> quant(std::vector<float>(gridn.size(), lb), std::vector<float>(gridn.size(), ub), gridn);
//    quant.set_sample(sample);

  timer::Timer time_cpp11;
  time_cpp11.reset();

  std::vector<double> temp1(gridn.size());
  std::vector<double> temp2(temp1.size());

//  for(size_t i = 0; i != nrolls; ++i)
//  {
//    for(size_t j = 0; j != temp1.size(); j++)
//    {
//      temp1[j] = ureal01(generator);
//    }
//    quant.transform(temp1,temp2);
//    values01.push_back(temp1);
//    sampled.push_back(temp2);
//  }


  size_t n = 2048;
  for(size_t x = 0; x != n; ++x)
  {
    for(size_t y = 0; y != n; ++y)
    {
//      temp1[0] = x / static_cast<double>(n) + 1.0/(2.0*n);
//      temp1[1] = y / static_cast<double>(n) + 1.0/(2.0*n);

      temp1[0] = x / static_cast<double>(n - 1);
      temp1[1] = y / static_cast<double>(n - 1);

//      std::cout << temp1[1] << std::endl;
//      std::cin.get();

      quant.transform(temp1,temp2);
//      values01.push_back(temp1);
//      sampled.push_back(temp2);

      size_t min_index = 0;
      double min_distance = std::sqrt(std::pow(temp2[0] - real_sample[0][0], 2.0) + std::pow(temp2[1] - real_sample[0][1], 2.0));
      for(size_t i = 1; i < real_sample.size(); i++)
      {
        double distance = std::sqrt(std::pow(temp2[0] - real_sample[i][0], 2.0) + std::pow(temp2[1] - real_sample[i][1], 2.0));
        if(distance < min_distance)
        {
          min_distance = distance;
          min_index = i;
        }
      }
      ++number[min_index];

      sampled_n[min_index].push_back(temp2);
      values_n[min_index].push_back(temp1);
    }
  }

  size_t nrolls = 1000;
  for(size_t k = 0; k != nrolls; ++k)
  {

    for(size_t j = 0; j != temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }

    quant.transform(temp1,temp2);
//    values01.push_back(temp1);
//    sampled.push_back(temp2);

    size_t min_index = 0;
    double min_distance = std::sqrt(std::pow(temp2[0] - real_sample[0][0], 2.0) + std::pow(temp2[1] - real_sample[0][1], 2.0));
    for(size_t i = 1; i < real_sample.size(); i++)
    {
      double distance = std::sqrt(std::pow(temp2[0] - real_sample[i][0], 2.0) + std::pow(temp2[1] - real_sample[i][1], 2.0));
      if(distance < min_distance)
      {
        min_distance = distance;
        min_index = i;
      }
    }
    ++number[min_index];

    sampled_n[min_index].push_back(temp2);
    values_n[min_index].push_back(temp1);
  }

  for(const auto & i : number)
  {
    std::cout << std::fixed << i << std::endl;
  }

//  for(size_t i = 0; i != 10000; ++i)
//  {
//    for(size_t j = 0; j != temp1.size(); j++)
//    {
//      temp1[j] = ureal01(generator);
//    }
//    quant.transform(temp1,temp2);
//    values01.push_back(temp1);
//    sampled.push_back(temp2);
//  }
  data_io::write_default2d("maps/sampled_implicit.dat", sampled, 10);

  for(size_t i = 0; i != sampled_n.size(); i++)
  {
    std::vector<std::vector<double>> hull;
    convexHull(sampled_n[i], hull);
    sampled_n[i].clear();
    sampled_n[i].shrink_to_fit();
    data_io::write_default2d("maps/sampled_lines" + std::to_string(i) + ".dat", hull, 10);

    //data_io::write_default2d("maps/sampled_lines" + std::to_string(i) + ".dat", sampled_n[i], 10);
  }
  for(size_t i = 0; i != values_n.size(); i++)
  {
    std::vector<std::vector<double>> hull;
    convexHull(values_n[i], hull);
    values_n[i].clear();
    values_n[i].shrink_to_fit();
    data_io::write_default2d("maps/value_lines" + std::to_string(i) + ".dat", hull, 10);

    //data_io::write_default2d("maps/sampled_lines" + std::to_string(i) + ".dat", sampled_n[i], 10);
  }
}

//


template <typename T>
bool increase(const std::vector<std::vector<T>> &v, std::vector<size_t> &it)
{
  for(size_t i = 0, size = it.size(); i != size; i++)
  {
    const size_t index = size - 1 - i;
    ++it[index];
    if(it[index] == v[index].size())
    {
      it[index] = 0;
    }
    else
    {
      return true;
    }
  }
  return false;
}


template <typename T>
std::vector<T> get_line(const std::vector<std::vector<T>> &v, std::vector<size_t> &it)
{
  std::vector<T> rez(v.size());
  for(size_t i = 0, size = v.size(); i != size; i++)
  {
    rez[i] = v[i][it[i]];
  }
  return rez;
}

template <typename T>
std::vector<std::vector<T>> iterate(const std::vector<std::vector<T>> &v)
{
  std::vector<size_t> it(v.size(), 0);
  std::vector<std::vector<T>> values;
  do
  {
    values.push_back(get_line(v, it));
  }
  while(increase(v, it));
  return values;
}


std::vector<double> worst_space(std::vector<size_t> gridN, std::vector<float> lb, std::vector<float> ub, size_t Nrolls, size_t seed_append)
{
  std::mt19937_64 generator;
  generator.seed(1 + gridN.size() + gridN.front() + seed_append);
  std::uniform_real_distribution<float> ureal01(0.0,1.0);

  typedef mveqf::trie_based::TrieBased<mveqf::trie_based::NodeCount<std::uint8_t>,std::uint8_t> sample_type;
  std::shared_ptr<sample_type> sample = std::make_shared<sample_type>();

  std::vector<std::vector<std::uint8_t>> variable_values(gridN.size());
  for(size_t i = 0; i != gridN.size(); i++)
  {
    variable_values[i].resize(gridN[i]);
    for(size_t j = 0; j != gridN[i]; j++)
    {
      variable_values[i][j] = j;
    }
  }
  std::vector<std::vector<std::uint8_t>> sample_int = iterate(variable_values);

  for(size_t i = 0; i != sample_int.size(); i++)
  {
    sample->insert(sample_int[i]);
  }

  std::vector<std::vector<float> > values01;

  std::vector<float> temp1(gridN.size());
  std::vector<float> temp2(temp1.size());

  for(size_t i = 0; i != Nrolls; ++i)
  {
    for(size_t j = 0; j != temp1.size(); j++)
    {
      temp1[j] = ureal01(generator);
    }
    values01.push_back(temp1);
  }

  std::vector<std::vector<float> > sampled(values01.size(), std::vector<float>(values01.front().size()));

  timer::Timer time_all_trans, full_time;
  std::vector<double> result;

  mveqf::ExplicitQuantile<std::uint8_t, float> quant_expl(lb, ub, gridN);
  quant_expl.set_sample(sample_int);
  full_time.reset();
  time_all_trans.reset();
  for(size_t i = 0; i != 10; i++)
    quant_expl.transform(values01[i], sampled[i]);
  result.push_back(time_all_trans.elapsed_seconds());
  result.push_back(result.back()/double(10));
  result.push_back(full_time.elapsed_seconds());
  for(size_t i = 0; i != sampled.size(); i++)
  {
    for(size_t j = 0; j != sampled[i].size(); j++)
    {
      if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
      {
        std::cout << "beyond bounds" << std::endl;
        std::cout << sampled[i][j] << std::endl;
      }
    }
  }

  mveqf::ImplicitQuantile<std::uint8_t, float> quant_impl(lb, ub, gridN);
  full_time.reset();
  quant_impl.set_sample_shared_and_fill_count(sample);
  time_all_trans.reset();
  for(size_t i = 0; i != values01.size(); i++)
    quant_impl.transform(values01[i], sampled[i]);
  result.push_back(time_all_trans.elapsed_seconds());
  result.push_back(result.back()/double(sampled.size()));
  result.push_back(full_time.elapsed_seconds());
  for(size_t i = 0; i != sampled.size(); i++)
  {
    for(size_t j = 0; j != sampled[i].size(); j++)
    {
      if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
      {
        std::cout << "beyond bounds" << std::endl;
        std::cout << sampled[i][j] << std::endl;
      }
    }
  }

  mveqf::ImplicitQuantileSorted<std::uint8_t, float> quant_impls(lb, ub, gridN);
  full_time.reset();
  quant_impls.set_sample_shared_and_fill_count(sample);
  time_all_trans.reset();
  for(size_t i = 0; i != values01.size(); i++)
    quant_impls.transform(values01[i], sampled[i]);
  result.push_back(time_all_trans.elapsed_seconds());
  result.push_back(result.back()/double(sampled.size()));
  result.push_back(full_time.elapsed_seconds());
  for(size_t i = 0; i != sampled.size(); i++)
  {
    for(size_t j = 0; j != sampled[i].size(); j++)
    {
      if(sampled[i][j] < (lb[j] - 0.001) || sampled[i][j] > (ub[j] + 0.001))
      {
        std::cout << "beyond bounds" << std::endl;
        std::cout << sampled[i][j] << std::endl;
      }
    }
  }
  return result;
}


void worst_space_test_grid_2d()
{
  size_t dim = 2, tries = 50, nrolls = 1e6;
  for(size_t grid_size = 1; grid_size != 100 + 1; grid_size++)
  {
    std::cout << grid_size << '\t';
    std::vector<double> times;
    for(size_t i = 0; i != tries; i++)
    {
      std::vector<size_t> g(dim, grid_size);
      std::vector<float> lb(dim, -1.0);
      std::vector<float> ub(dim, 1.0);
      auto rez = worst_space(g, lb, ub, nrolls, i);

      times.resize(rez.size());
      for(size_t j = 0; j != times.size(); j++)
      {
        times[j] += rez[j];
      }
    }
    for(size_t j = 0; j != times.size(); j++)
    {
      times[j] /= double(tries);
    }
    for(const auto &i : times)
      std::cout << std::scientific << i << '\t';
    std::vector<size_t> g(dim, grid_size);
    std::cout << std::accumulate(g.begin(), g.end(), 1, std::multiplies<size_t>());
    std::cout << '\t' << grid_size << '^' << dim << std::endl;
  }
}
