/**************************************************************************

   Copyright © 2020 Sergey Poluyan <svpoluyan@gmail.com>

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
#ifndef SAMPLE_H
#define SAMPLE_H

namespace mveqf
{
	template <typename TIndex>
	struct Sample
	{
		virtual void set_dimension(size_t dim) = 0;
		virtual size_t get_dimension() const = 0;

		virtual void insert(const std::vector<TIndex> &key) = 0;
		virtual void insert(const std::vector<TIndex> &key, size_t number) = 0;
		virtual bool search(const std::vector<TIndex> &key) const = 0;

		virtual void fill_tree_count() = 0;

		virtual size_t get_link_count() const = 0;
		virtual size_t get_node_count() const = 0;

		virtual ~Sample() = default;
	};
}

#endif
