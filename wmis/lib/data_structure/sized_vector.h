/**
 * sized_vector.h
 * Purpose: Presized vector which allows fast clearing
 *
 *****************************************************************************/

#ifndef _SIZED_VECTOR_H
#define _SIZED_VECTOR_H

#include <vector>

template<class T>
class sized_vector {
private:
	std::vector<T> data;
	size_t counter = 0;

public:
	using iterator = std::vector<NodeID>::iterator;
	using const_iterator = std::vector<NodeID>::const_iterator;

	sized_vector(size_t min_capacity = 0) : data(min_capacity) { }
	sized_vector(size_t min_capacity, const T& value) : data(min_capacity, value) { }


	T& operator[] (size_t index) { return data[index]; }
	T& front() { return data[0]; }
	T& back() { return data[counter - 1]; }
	const T& operator[] (size_t index) const { return data[index]; }
	const T& front() const { return data[0]; }
	const T& back() const { return data[counter - 1]; }

	iterator begin() { return data.begin(); }
	iterator end() { return data.begin() + counter; }
	const_iterator begin() const { return data.begin(); }
	const_iterator end() const { return data.begin() + counter; }
	const_iterator cbegin() const { return data.cbegin(); }
	const_iterator cend() const { return data.cbegin() + counter; }

	size_t size() const noexcept { return counter; }
	size_t capacity() const { return data.capacity(); }
	bool empty() const noexcept { return counter == 0; }

	void clear() noexcept { counter = 0; }

	void resize(size_t capacity) {
		data.resize(capacity);

		if (counter >= capacity)
			counter = capacity - 1;
	}

	void resize(size_t capacity, const T& value) {
		data.resize(capacity, value);

		if (counter >= capacity)
			counter = capacity - 1;
	}

	void set_size(size_t size) {
		counter = size;
	}

	void push_back(const T& value) { data[counter++] = value; }
	void push_back(T&& value) { data[counter++] = std::move(value); }
	void pop_back() { --counter; }

	void remove(iterator iter) {
		std::swap(*iter, back());
		pop_back();
	}

	void remove(const_iterator iter) {
		std::swap(*iter, back());
		pop_back();
	}

	void remove(const size_t pos) {
		std::swap(data[pos], back());
		pop_back();
	}

	template <class... Args>
	void emplace_back(Args&&... args) {
		data[counter++] = T(std::forward<Args>(args)...);
	}

	void swap(sized_vector & other) {
		std::swap(data, other.data);
		std::swap(counter, other.counter);
	}
};

#endif
