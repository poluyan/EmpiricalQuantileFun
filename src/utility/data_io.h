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

#ifndef DATA_IO_H
#define DATA_IO_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

namespace data_io
{

template<template<typename, typename...> class Container, class T, typename... Params>
void write_default1d(std::string fname, Container<T, Params...> const& u, size_t step = 1, size_t prec = 10)
{
    std::ofstream fOut;
    fOut.open(fname.c_str());
    if(!fOut.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }
    fOut.precision(prec);
    for(size_t i = 0; i != u.size(); i += step)
    {
        fOut << std::scientific << i << '\t' << u[i] << std::endl;
    }
    fOut.close();
    std::cout << fname << std::endl;
}

template<template<typename, typename...> class Container, class T, typename... Params>
void write_default2d(std::string fname, Container<T, Params...> const& u, size_t prec)
{
    std::ofstream fOut;
    fOut.open(fname.c_str());
    if(!fOut.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }
    fOut.precision(prec);
    for(const auto & i : u)
    {
        for(const auto & j : i)
        {
            fOut << std::scientific << j << '\t';
        }
        fOut << std::endl;
    }
    fOut.close();
    std::cout << fname << std::endl;
}

template<typename T, typename A>
void load_grid_and_sample(std::string fname_grid,
                          std::string fname_sample,
                          std::vector<size_t> &grid,
                          std::vector<std::vector<T, A>> &sample)
{
    std::ifstream gridIn;
    gridIn.open(fname_grid.c_str());
    if(!gridIn.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }

    while(!gridIn.eof())
    {
        size_t item = 0;
        gridIn >> item;
        grid.push_back(item);
    }
    grid.pop_back();
    gridIn.close();
    std::cout << fname_grid << '\t' << grid.size() << std::endl;

    std::ifstream sampleIn;
    sampleIn.open(fname_sample.c_str());
    if(!sampleIn.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }

    std::vector<T, A> push2data(grid.size());
    while(!sampleIn.eof())
    {
        size_t item = 0;
        for(size_t i = 0; i != push2data.size(); i++)
        {
            sampleIn >> item;
            push2data[i] = item;
        }
        sample.push_back(push2data);
    }
    sample.pop_back();

    sampleIn.close();
    std::cout << fname_sample << '\t' << sample.size() << std::endl;
}

}

#endif
