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
#ifndef CSTVECT_H
#define CSTVECT_H

#include <algorithm>
#include <iterator>

namespace mveqf
{
	namespace cst
	{
		template<class vect>
		class vector_iterator
		{
		public:
			using iterator_category = std::random_access_iterator_tag;
			using value_type = typename vect::value_type;
			using difference_type = std::ptrdiff_t;//typename vect::difference_type;
			using pointer = typename vect::pointer;
			using reference = typename vect::reference;
		public:
			vector_iterator()
			{
				ptr = nullptr;
			}
			vector_iterator(const vector_iterator &other)
			{
				ptr = other.ptr;
			}
			vector_iterator(pointer pt) : ptr(pt) {}
			reference operator*() const
			{
				return *ptr;
			}
			pointer operator->() const
			{
				return &(**this);
			}
			vector_iterator& operator++()
			{
				++ptr;
				return *this;
			}
			vector_iterator operator++(int)
			{
				vector_iterator t(*this);
				++(*this);
				return t;
			}
			vector_iterator& operator--()
			{
				--ptr;
				return *this;
			}
			vector_iterator operator--(int)
			{
				vector_iterator t(*this);
				--(*this);
				return t;
			}
			vector_iterator& operator+=(difference_type diff)
			{
				ptr += diff;
				return *this;
			}
			vector_iterator operator+(difference_type diff) const
			{
				vector_iterator tmp(*this);
				return tmp += diff;
			}
			vector_iterator& operator-=(difference_type diff)
			{
				ptr -= diff;
				return *this;
			}
			vector_iterator operator-(difference_type diff) const
			{
				vector_iterator tmp(*this);
				return tmp -= diff;
			}
			difference_type operator-(const vector_iterator& other) const
			{
				return ptr - other.ptr;
			}
			bool operator==(const vector_iterator& other) const
			{
				return ptr == other.ptr;
			}
			bool operator!=(const vector_iterator& other) const
			{
				return !(*this == other);
			}
			bool operator<(const vector_iterator& other) const
			{
				return ptr < other.ptr;
			}
			bool operator>(const vector_iterator& other) const
			{
				return other < *this;
			}
			bool operator<=(const vector_iterator& other) const
			{
				return !(other < *this);
			}
			bool operator>=(const vector_iterator& other) const
			{
				return !(*this < other);
			}
			vector_iterator& operator=(const vector_iterator& other)
			{
				ptr = other.ptr;
				return *this;
			}
		protected:
			pointer ptr;
		};

		template<class T>
		class vector
		{
		public:

			using difference_type = int;
			using value_type = T;
			using pointer = value_type*;
			using reference = value_type&;
			using const_reference = const value_type&;
			using size_type = std::size_t; // grid size!!
			using iter = vector_iterator<vector<value_type>>;

		public:
			vector();
			vector(const vector& cp);
			vector& operator=(const vector&) = delete;

			~vector();

			reference operator[](size_type idx);
			const_reference operator[](size_type idx) const;

			reference front();
			const_reference front() const;
			reference back();
			const_reference back() const;

			iter begin() const;
			iter end() const;

			void push_back(const value_type &val);
			void push_back(value_type &&val);

			template <class ... Args>
			void emplace_back(Args && ... args);

			void pop_back();

			inline void shrink_to_fit() const noexcept;

			size_type size() const noexcept;
			bool empty() const noexcept;
			void clear();

		protected:
			pointer	data;
			size_type	vec_sz;

			inline void reallocate_up();
		};

		template <typename T>
		inline void vector<T>::reallocate_up()
		{
			pointer tarr = new value_type [vec_sz];
			for(size_type i = 0; i < vec_sz - 1; ++i)
				::new(static_cast<void*>(&tarr[i])) value_type(std::move(data[i]));
			delete [] data;
			data = tarr;
		}

		template <typename T>
		vector<T>::vector(const vector& cp)
		{
			pointer tarr = new value_type [cp.vec_sz];
			for(size_type i = 0; i < cp.vec_sz; ++i)
				::new(static_cast<void*>(&tarr[i])) value_type(std::move(cp.data[i]));
			delete [] data;
			data = tarr;
		}

		template <typename T>
		inline void vector<T>::shrink_to_fit() const noexcept
		{
			// the main idea of this container is to not use capacity in any way
		}

		template <typename T>
		typename vector<T>::size_type vector<T>::size() const noexcept
		{
			return vec_sz;
		}

		template <typename T>
		bool vector<T>::empty() const noexcept
		{
			return vec_sz == 0;
		}

		template <typename T>
		void vector<T>::clear()
		{
			vec_sz = 0;
			delete [] data;
			data = nullptr;
		}

		template <typename T>
		typename vector<T>::reference vector<T>::operator[](vector<T>::size_type idx)
		{
			return *(this->data + idx);
		}

		template <typename T>
		typename vector<T>::const_reference vector<T>::operator[](vector<T>::size_type idx) const
		{
			return *(this->data + idx);
		}

		template <typename T>
		typename vector<T>::reference vector<T>::front()
		{
			return *begin();
		}

		template <typename T>
		typename vector<T>::const_reference vector<T>::front() const
		{
			return *begin();
		}

		template <typename T>
		typename vector<T>::reference vector<T>::back()
		{
			return *(end() - 1);
		}

		template <typename T>
		typename vector<T>::const_reference vector<T>::back() const
		{
			return *(end() - 1);
		}

		template <typename T>
		vector<T>::vector() : data(nullptr), vec_sz(0) {}

		template <typename T>
		void vector<T>::pop_back()
		{
			if(vec_sz > 1)
			{
				--vec_sz;
				pointer tarr = new value_type [vec_sz];
				for(size_type i = 0; i < vec_sz; ++i)
					::new(static_cast<void*>(&tarr[i])) value_type(std::move(data[i]));
				delete [] data;
				data = tarr;
			}
			else if(vec_sz == 1)
			{
				--vec_sz;
				delete [] data;
				data = nullptr;
			}
		}

		template <typename T>
		vector<T>::~vector()
		{
			delete [] data;
		}

		template <typename T>
		typename vector<T>::iter vector<T>::begin() const
		{
			return iter(data);
		}

		template <typename T>
		typename vector<T>::iter vector<T>::end() const
		{
			return iter(data + vec_sz);
		}

		template <typename T>
		void vector<T>::push_back(const vector<T>::value_type &val)
		{
			emplace_back(val);
		}

		template <typename T>
		void vector<T>::push_back(vector<T>::value_type &&val)
		{
			emplace_back(std::move(val));
		}

		template <typename T>
		template <class ... Args>
		void vector<T>::emplace_back(Args && ... args)
		{
			if(vec_sz < std::numeric_limits<size_type>::max())
			{
				++vec_sz;
				reallocate_up();
				data[vec_sz - 1] = std::move(value_type(std::forward<Args>(args) ...));
			}
		}
	}
}

#endif
