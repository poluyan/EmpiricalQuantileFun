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

#ifndef PRINT2FILE_H
#define PRINT2FILE_H

#include <string>
#include <vector>
#include <fstream>

template<typename T, typename A>
void print2file(std::string fname, std::vector<T, A> u, int step)
{
    std::ofstream fOut;
    fOut.open(fname.c_str());
    if(!fOut.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }
    fOut.precision(10);
    for(size_t i = 0; i < u.size(); i += step)
    {
        fOut << std::scientific << i << '\t' << u[i] << std::endl;
    }
    fOut.close();
    std::cout << fname << std::endl;
}

template<typename T, typename A>
void print2file2d(std::string fname, std::vector<std::vector<T, A> >& u)
{
    std::ofstream fOut;
    fOut.open(fname.c_str());
    if(!fOut.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }
    fOut.precision(10);
    for(size_t i = 0; i < u.size(); i++)
    {
        for(size_t j = 0; j < u[i].size(); j++)
        {
            fOut << std::scientific << u[i][j] << '\t';
        }
        fOut << std::endl;
    }
    fOut.close();
    std::cout << fname << std::endl;
}

#endif