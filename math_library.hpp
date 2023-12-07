#pragma once

#ifndef MATH_LIBRARY_HPP
#define MATH_LIBRARY_HPP

#include <iostream>
#include <regex>
#include <iomanip>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <map>
#include <tuple>
#include <numeric>
#include <conio.h>
#include <vector>
#include <algorithm>
#include <format>
#include <bitset>

template<long long p, long long n>
class GaloisField {
private:

	using field_element = std::int64_t;
	field_element _elem;

private:
	static field_element bin_pow(field_element element, long long pow) noexcept {
		field_element mult = 1;
		while (pow > 0) {
			if (pow & 1) {
				mult *= element;
			}
			element *= element;
			pow >>= 1;
		}
		return mult;
	}
public:
	GaloisField() : _elem(0)
	{}
	GaloisField(std::int64_t _field_elem) : _elem(_field_elem% p)
	{}
	~GaloisField() = default;
public:
	friend std::istream& operator >> (std::istream& in, GaloisField& rhs) noexcept {
		std::cout << "Input element of Galois Field\n";
		in >> rhs._elem;
		rhs._elem %= p;
		return in;
	}
	friend std::ostream& operator << (std::ostream& out, const GaloisField& rhs) noexcept {
		out << rhs._elem << "";
		return out;
	}
public:
	std::int64_t get_elem() const {
		return _elem;
	}
public:
	bool operator == (const GaloisField& rhs) {
		return _elem == rhs._elem;
	}
	bool operator != (const GaloisField& rhs) {
		return _elem != rhs._elem;
	}
public:
	constexpr decltype(auto) operator +(const GaloisField& rhs) const noexcept {
		auto copy = *this;
		copy._elem += rhs._elem;
		return copy._elem % p;
	}
	constexpr decltype(auto) operator -(const GaloisField& rhs) const noexcept {
		auto copy = *this;
		if (copy._elem > rhs._elem) {
			return GaloisField((copy._elem - rhs._elem) % p);
		}
		return GaloisField((copy._elem - rhs._elem + p) % p);
	}

	constexpr decltype(auto) operator *(const GaloisField& rhs) const noexcept {
		auto copy = *this;
		copy._elem *= rhs._elem;
		return copy._elem % p;
	}
	constexpr decltype(auto) operator /(const GaloisField& rhs) const noexcept {
		auto copy = *this;
		copy._elem /= rhs._elem;
		return copy;
	}
	constexpr decltype(auto) operator +=(const GaloisField& rhs) noexcept {
		_elem += rhs._elem;
		_elem %= p;
		return *this;
	}

	constexpr decltype(auto) operator -=(const GaloisField& rhs) noexcept {
		_elem -= rhs._elem;
		if (_elem < 0) {
			_elem += p;
		}
		return *this;
	}

	constexpr decltype(auto) operator *=(const GaloisField& rhs) noexcept {
		_elem *= rhs._elem;
		_elem %= p;
		return *this;
	}
	constexpr decltype(auto) operator /=(const GaloisField& rhs) noexcept {
		auto result_bin_pow = GaloisField::bin_pow(rhs._elem, p - 2);
		_elem *= result_bin_pow;
		_elem %= p;
		return *this;
	}
};

template <typename T>
struct is_galois_field : std::false_type {};

template <long long p, long long n>
struct is_galois_field<GaloisField<p, n>> : std::true_type {
	static constexpr long long _field_characteristic = p;
	static constexpr long long _vector_size = n;
};

template<class T>
class Polynom {
private:
	std::vector<T> _polynom;
public:
	Polynom() {}
	Polynom(std::initializer_list<T> _init_list) noexcept
		: _polynom(std::move_if_noexcept(_init_list))
	{}
	~Polynom() {}
private:
	T get_lead_degree() const {
		return _polynom.back();
	}
public:
	Polynom& operator = (const Polynom& rhs) {
		if (this != &rhs) {
			_polynom.assign(rhs._polynom.begin(), rhs._polynom.end());
		}
		return *this;
	}
public:
	decltype(auto) operator + (const Polynom& rhs) const {
		auto copy = *this;
		for (std::size_t i = 0; i < std::min(rhs._polynom.size(), copy._polynom.size()); ++i) {
			copy._polynom[i] += rhs._polynom[i];
		}
		return copy;
	}
	decltype(auto) operator +=(const Polynom& rhs) {
		for (std::size_t i = 0; i < std::min(rhs._polynom.size(), _polynom.size()); ++i) {
			_polynom[i] += rhs._polynom[i];
		}
		return *this;
	}
	decltype(auto) operator - (const Polynom& rhs) const {
		auto copy = *this;
		for (std::size_t i = 0; i < std::min(rhs._polynom.size(), copy._polynom.size()); ++i) {
			copy._polynom[i] -= rhs._polynom[i];
		}
		return copy;
	}
	decltype(auto) operator -=(const Polynom& rhs) {
		for (std::size_t i = 0; i < std::min(rhs._polynom.size(), _polynom.size()); ++i) {
			_polynom[i] -= rhs._polynom[i];
		}
		return *this;
	}
	decltype(auto) operator * (const Polynom& rhs) {
		auto copy = *this;
		std::vector<T> result(copy._polynom.size() + rhs._polynom.size() - 1, 0);
		for (std::size_t i = 0; i < copy._polynom.size(); ++i) {
			for (std::size_t j = 0; j < rhs._polynom.size(); ++j) {
				result[i + j] += copy._polynom[i] * rhs._polynom[j];
			}
		}
		copy._polynom = std::move(result);
		return copy;
	}
	decltype(auto) operator *= (const Polynom& rhs) {
		std::vector<T> result(_polynom.size() + rhs._polynom.size() - 1, 0);
		for (std::size_t i = 0; i < _polynom.size(); ++i) {
			for (std::size_t j = 0; j < rhs._polynom.size(); ++j) {
				result[i + j] += _polynom[i] * rhs._polynom[j];
			}
		}
		_polynom = std::move(result);
		return *this;
	}
	decltype(auto) operator / (const Polynom& rhs) {

		constexpr int smart_constant = 10;

		auto divident = *this;
		auto divide_into = rhs;

		auto remainder = Polynom();
		auto quitient = Polynom();

		quitient._polynom.resize(std::max(divident._polynom.size(), divide_into._polynom.size()) - 1);

		auto divident_sz = divident._polynom.size();
		auto divide_sz = divide_into._polynom.size();

		auto coefficient_of_correlation = 0;

		while (divident._polynom.size() >= divide_into._polynom.size()) {

			T first_polynom_divident_coefficient = divident._polynom.back();
			T second_polynom_divident_coefficient = divide_into._polynom.back();
			auto degree_of_quitient = divident._polynom.size() - divide_into._polynom.size();
			auto coefficient = first_polynom_divident_coefficient / second_polynom_divident_coefficient;
			quitient._polynom.insert(quitient._polynom.cbegin() + degree_of_quitient, coefficient);

			std::vector<T> _micro_polynom;
			_micro_polynom.resize(degree_of_quitient + 1);
			_micro_polynom.insert(_micro_polynom.begin() + degree_of_quitient, coefficient);
			auto new_polynom{ _micro_polynom };
			auto second_new_polynom = Polynom();
			second_new_polynom._polynom.resize(new_polynom.size() + smart_constant * 2);

			for (int i = divide_into._polynom.size() - 1; i >= 0; --i) {
				auto weight = 0;
				if (divide_into._polynom[i] != 0) {
					weight = i;
					second_new_polynom._polynom[degree_of_quitient + weight] = new_polynom[degree_of_quitient] * divide_into._polynom[i];
				}
			}
			divident -= second_new_polynom;
			while (divident._polynom.size() > 0 && divident._polynom.back() == 0) {
				divident._polynom.pop_back();
			}
			++coefficient_of_correlation;
		}

		--coefficient_of_correlation;
		auto _cpy_quitient = Polynom();
		_cpy_quitient._polynom.resize(quitient._polynom.size());

		for (int i = quitient._polynom.size() - 1; i >= 0; --i) {
			if (quitient._polynom[i] != 0) {
				_cpy_quitient._polynom[i - coefficient_of_correlation] = quitient._polynom[i];
				--coefficient_of_correlation;
			}
		}

		remainder = divident;

		return _cpy_quitient;
	}
	decltype(auto) operator % (const Polynom& rhs) {
		constexpr int smart_constant = 10;

		auto divident = *this;
		auto divide_into = rhs;

		auto remainder = Polynom();
		auto quitient = Polynom();

		quitient._polynom.resize(std::max(divident._polynom.size(), divide_into._polynom.size()) - 1);

		auto divident_sz = divident._polynom.size();
		auto divide_sz = divide_into._polynom.size();

		auto coefficient_of_correlation = 0;

		while (divident._polynom.size() >= divide_into._polynom.size()) {

			T first_polynom_divident_coefficient = divident._polynom.back();
			T second_polynom_divident_coefficient = divide_into._polynom.back();
			auto degree_of_quitient = divident._polynom.size() - divide_into._polynom.size();
			auto coefficient = first_polynom_divident_coefficient / second_polynom_divident_coefficient;
			quitient._polynom.insert(quitient._polynom.cbegin() + degree_of_quitient, coefficient);

			std::vector<T> _micro_polynom;
			_micro_polynom.resize(degree_of_quitient + 1);
			_micro_polynom.insert(_micro_polynom.begin() + degree_of_quitient, coefficient);
			auto new_polynom{ _micro_polynom };
			auto second_new_polynom = Polynom();
			second_new_polynom._polynom.resize(new_polynom.size() + smart_constant * 2);

			for (int i = divide_into._polynom.size() - 1; i >= 0; --i) {
				auto weight = 0;
				if (divide_into._polynom[i] != 0) {
					weight = i;
					second_new_polynom._polynom[degree_of_quitient + weight] = new_polynom[degree_of_quitient] * divide_into._polynom[i];
				}
			}
			divident -= second_new_polynom;
			while (divident._polynom.size() > 0 && divident._polynom.back() == 0) {
				divident._polynom.pop_back();
			}
			++coefficient_of_correlation;
		}

		--coefficient_of_correlation;
		auto _cpy_quitient = Polynom();
		_cpy_quitient._polynom.resize(quitient._polynom.size());

		for (int i = quitient._polynom.size() - 1; i >= 0; --i) {
			if (quitient._polynom[i] != 0) {
				_cpy_quitient._polynom[i - coefficient_of_correlation] = quitient._polynom[i];
				--coefficient_of_correlation;
			}
		}

		remainder = divident;

		return remainder;
	}
public:
	friend std::istream& operator>>(std::istream& in, Polynom<T>& rhs) {
		if constexpr (is_galois_field<T>::value) {
			rhs._polynom.resize(is_galois_field<T>::_vector_size);
			for (auto&& x : rhs._polynom) {
				std::cin >> x;
			}
		}
		else {
			try {
				std::int64_t sz;
				std::cin >> sz;
				if (sz < 0) {
					throw std::exception("Unnatural polynom size...\n");
				}
				rhs._polynom.resize(sz);
				for (auto&& x : rhs.polynom) {
					std::cin >> x;
				}
			}
			catch (const std::exception& ex) {
				std::cout << ex.what();
			}
		}
		return in;
	}
	friend std::ostream& operator<<(std::ostream& out, const Polynom& rhs) {
		if constexpr (is_galois_field<T>::value) {
			out << "Galois polynom: \n";
			for (std::size_t i = 0; i < rhs._polynom.size(); ++i) {
				if (i < rhs._polynom.size() - 1) {
					if (rhs._polynom[i].get_elem() != 0) {
						out << rhs._polynom[i] << "x^" << i << "+";
					}
				}
				else {
					if (rhs._polynom[i].get_elem() != 0) {
						out << rhs._polynom[i] << "x^" << i;
					}
				}
			}
		}
		else {
			std::cout << "Not Galois Field\n";
			for (std::size_t i = 0; i < rhs._polynom.size(); ++i) {
				if (i < rhs._polynom.size() - 1) {
					if (rhs._polynom[i]/*get_elem()*/ != 0) {
						out << rhs._polynom[i] << "x^" << i << "+";
					}
				}
				else {
					if (rhs._polynom[i] != 0) {
						out << rhs._polynom[i] << "x^" << i;
					}
				}
			}
		}
		return out;
	}
public:
	Polynom gcd(Polynom& a, Polynom& b) {
		while (b._polynom.size() > 1) {
			a = a % b;
			auto t = a;
			a = b;
			b = t;
		}
		return a;
	}
};

template<long long order_of_linear_recurrent_sequence>
class linear_recurrent_sequence_byte {
public:
private:
	std::map<int, int> mp;
	std::int64_t real_size{};
	std::int32_t _needs_block_size;
	std::vector<std::byte> _stock;
	std::vector<std::byte> _initial_padding;
	std::vector<std::int64_t> _masks;
public:
	linear_recurrent_sequence_byte() = default;
	linear_recurrent_sequence_byte(std::initializer_list<std::byte> _initializer_list) noexcept :
		_stock(std::move_if_noexcept(_initializer_list))
	{
		_needs_block_size = linear_recurrent_sequence_byte::calculate_byte(0);
	};
	~linear_recurrent_sequence_byte() = default;
private:
	static long long calculate_byte(long long x) {
		while (x * 8 < order_of_linear_recurrent_sequence) {
			++x;
		}
		return x * 8;
	}
	static long long my_abs(long long x) noexcept {
		if (x < 0) {
			return x * -1;
		}
		return x;
	}
	static std::byte result_of_xor(std::byte x) {
		x ^= (x >> 4);
		x ^= (x >> 2);
		x ^= (x >> 1);
		return (x & static_cast<std::byte>(1));
	}
public:
	void asd() {
		for (std::size_t i = 0; i < _masks.size(); ++i) {
			_masks[i] -= 500;
		}
	}
	void _set_value(std::byte& _byte) {
		_initial_padding.clear();
		_initial_padding.push_back(_byte);
		_needs_block_size = linear_recurrent_sequence_byte::calculate_byte(0);
	}
	friend std::istream& operator >> (std::istream& in, linear_recurrent_sequence_byte& rhs) {

		auto _check_block = 0;

		std::cout << "Input initiall padding of linear_recurrent_sequence\n";
		std::string sequence{};
		try {

			in >> sequence;

			if (sequence.size() != order_of_linear_recurrent_sequence) {
				throw std::runtime_error("size of sequence is greater than order of lfsr");
			}

			for (std::size_t i = order_of_linear_recurrent_sequence; i < _needs_block_size; ++i) {
				sequence[i] += '0';
			}

			std::byte _value_type_byte{};
			std::size_t _checker = 7;

			for (std::size_t i = 0; i < sequence.size(); ++i) {
				if (sequence[i] == '1') {
					auto _one_bit = sequence[i] - '0';
					_one_bit <<= _checker;
					_value_type_byte |= static_cast<std::byte>(_one_bit);
				}
				if ((i + 1) % 8 == 0) {
					_checker = 7;
					rhs._stock.push_back(static_cast<std::byte>(_value_type_byte));
					_value_type_byte ^= _value_type_byte;
					++_check_block;
				}
				else {
					--_checker;
				}
			}
			if (static_cast<std::int32_t>(_value_type_byte) != 0x0) {
				rhs._stock.push_back(static_cast<std::byte>(_value_type_byte));
			}
		}
		catch (std::runtime_error& e_r) {
			std::cout << e_r.what();
		}
		return in;
	}

	friend std::ostream& operator << (std::ostream& out, const linear_recurrent_sequence_byte& rhs) {
		int i = 1;
		for (auto&& x : rhs._stock) {
			std::bitset<8> _bits(static_cast<unsigned int>(x));
			out << std::format("Block number {}, check bits {}\n", i, _bits.to_string());
			++i;
		}
		return out;
	}


	void points_of_removal() {
		std::cout << "Input amount of removal poits\n";
		std::int32_t _amount{};
		std::cin >> _amount;
		if (_amount < 0) {
			throw std::runtime_error("amount of removal points can`t be negative");
		}
		_masks.reserve(_amount);
		short removal_point{};
		for (std::size_t i = 0; i < _amount; ++i) {
			std::cin >> removal_point;
			_masks.push_back(removal_point);
		}
	}

	template<long long N>
	void generate_sequence_with_help_of_points() {

		decltype(auto) how_many_iterations = N;
		decltype(auto) how_many_bits_not_in_the_block = _needs_block_size - order_of_linear_recurrent_sequence;
		decltype(auto) how_many_blocks = static_cast<int>(std::ceil(_needs_block_size / 8));

		real_size = (8 - ((how_many_iterations - how_many_bits_not_in_the_block) % 8));

		_initial_padding.assign(_stock.begin(), _stock.end());
		if (how_many_bits_not_in_the_block != 0) {
			--how_many_bits_not_in_the_block;
		}


		std::byte xor_sum = static_cast<std::byte>(0b0);

		while (how_many_iterations-- > 0) {
			long long state_block = how_many_blocks;
			if (state_block != 0) {
				--state_block;
			}
			for (std::size_t i = 0; i < _masks.size(); ++i) {
				auto _bit =
					(
						_masks[i] < 8 ? 7 - _masks[i] : 7 - (_masks[i] % 8)
					);
				auto _detect_block = std::ceil(_masks[i] / 8);
				auto _get_removal = static_cast<std::byte>(1 << _bit);
				auto _current_byte = _stock[_detect_block];
				xor_sum ^= static_cast<std::byte>(_current_byte & _get_removal);
			}

			xor_sum = (result_of_xor(xor_sum) & static_cast<std::byte>(1));

			_stock[state_block] |= static_cast<std::byte>(xor_sum << how_many_bits_not_in_the_block);

			xor_sum ^= xor_sum;

			std::transform(_masks.begin(), _masks.end(), _masks.begin(), [](auto& x) {return ++x; });

			if (how_many_bits_not_in_the_block == 0) {
				++how_many_blocks;
				how_many_bits_not_in_the_block = 7;
				_stock.resize(how_many_blocks);
			}
			else {
				--how_many_bits_not_in_the_block;
			}

			if (how_many_iterations == 0 && how_many_bits_not_in_the_block != 0) {
				_stock[state_block] |= static_cast<std::byte>(xor_sum << how_many_bits_not_in_the_block);
			}
		}
	}

	int find_pattern(const std::string& text, const std::string& pattern) {
		int textLength = text.length();
		int patternLength = pattern.length();
		//std::cout << "\n\n\n\n" << text;
		for (int i = 0; i < textLength - patternLength; ++i) {
			bool match = true;
			for (int j = 0; j < patternLength; ++j) {
				if (text[i + j] != pattern[j]) {
					match = false;
					break;
				}
			}
			if (match) {
				return i;
			}
		}

		return -1;
	}
	void print_map() {
		for (const auto& x : mp) {
			std::cout << x.first << " " << x.second << "\n";
		}
	}
	void find_period_of_linear_recurrent_sequence() {

		std::string text{};
		std::string pattern{};

		decltype(auto) get_p = (static_cast<int>(std::ceil(_needs_block_size / 8)) * 8) - order_of_linear_recurrent_sequence;
		/*pattern = std::bitset<8>(i).to_string();*/
		for (std::size_t i = 0; i < _initial_padding.size(); ++i) {
			pattern += std::bitset<8>(static_cast<int>(_initial_padding[i])).to_string();
		}

		pattern.erase(pattern.begin() + pattern.size() - get_p, pattern.end());

		for (std::size_t i = 0; i < 500; ++i) { //magic constant 500 as in asd method
			text += std::bitset<8>(static_cast<int>(_stock[0])).to_string();
		}
		text.erase(text.begin() + text.size() - real_size, text.end());
		text.erase(text.begin(), text.begin() + pattern.size());

		std::cout << "Period: " << find_pattern(text, pattern) + pattern.size() << "\t" << pattern << "\n";

		mp[find_pattern(text, pattern) + pattern.size()]++;
		//std::cout << "\n";
		//print_map();
	}
};

class BooleanFunction {
private:

	long long count_degree_of_two;
	std::vector<int> _boolean;

public:
	BooleanFunction() = default;
	BooleanFunction(std::initializer_list<int> _init) noexcept :
		_boolean(std::move_if_noexcept(_init))
	{
		count_degree_of_two = 0;
		while (_boolean.size() > (1 << count_degree_of_two)) {
			++count_degree_of_two;
		}
	}
private:
	std::string binary_view(unsigned int n, const int& degree)
	{
		std::string buffer{};
		buffer.reserve(std::numeric_limits<unsigned int>::digits);

		do
		{
			buffer += static_cast<char>('0' + n % 2);
			n = n / 2;
		} while (n > 0);

		while (buffer.size() < degree) {
			buffer += '0';
		}
		return std::string(buffer.crbegin(), buffer.crend());
	}
	int degree_of_polynom(const std::string& _polynom) {

		std::regex regular("([x0-9]+)");
		std::sregex_iterator begin(_polynom.begin(), _polynom.end(), regular);
		std::sregex_iterator end;
		std::unordered_map<std::string, int> unordered{};


		for (; begin != end; ++begin) {
			auto s = begin->str(1);
			unordered[s]++;
		}

		std::string new_string{};

		for (auto& x : unordered) {
			if (x.second == 1) {
				new_string += x.first;
				new_string += "+";
			}
		}
	
		std::sregex_iterator begin_(new_string.begin(), new_string.end(), regular);
		std::sregex_iterator end_;

		std::ptrdiff_t degree = 0;

		for (; begin_ != end_; ++begin_) {
			auto s = begin_->str(1);
			auto _x = count(s.begin(), s.end(), 'x');
			degree = std::max(degree, _x);
		}


		return degree;
	}
	decltype(auto) polynom_zhegalkin_pascal_(std::vector<int>& like_boolean) {
		std::vector<int> _polynom;
		while (like_boolean.size() != 0) {
			_polynom.push_back(like_boolean[0]);
			for (int i = 0; i < like_boolean.size() - 1; ++i) {
				like_boolean[i] ^= like_boolean[i + 1];
			}
			like_boolean.pop_back();
		}
		return _polynom;
	}

public:
	int weight_of_boolean_function() {
		int weight{};
		weight = std::count(_boolean.begin(), _boolean.end(), true);
		return weight;
	}

	void full_weight_structure() {

		std::map<std::string, int> calculating_bits{};

		std::cout << "a1a2a3a4   f\n";

		for (int i = 0; i < (1 << count_degree_of_two); ++i) {
			calculating_bits.insert(std::make_pair(binary_view(i, count_degree_of_two), _boolean[i]));
			std::cout << binary_view(i, count_degree_of_two) << std::setw(8) << _boolean[i];
			std::cout << "\n";
		}

		std::cout << "\nHow many query?\n";
		int t;
		std::cin >> t;
		while (t-- > 0) {

			std::cout << "\nInput how many fixed points result do you want to see, remember iteration start with zero\n";
			long long how_many_fixed_points{};
			std::cin >> how_many_fixed_points;
			std::vector<long long> fixed_points(how_many_fixed_points);
			std::cout << "Input fixed points\n";
			for (auto&& x : fixed_points) {
				std::cin >> x;
			}
			std::cout << "\n";
			auto degree = 0;
			for (const auto& [fixed_string, value] : calculating_bits) {
				std::string new_string = "";
				for (std::size_t i = 0; i < fixed_points.size(); ++i) {
					new_string += fixed_string[fixed_points[i]];
				}
				for (const auto& [secondary_string, boolean_value] : calculating_bits) {
					std::string your_string = "";
					for (std::size_t i = 0; i < fixed_points.size(); ++i) {
						your_string += secondary_string[fixed_points[i]];
					}
					if (new_string == your_string) {
						degree += boolean_value;
					}
				}
				std::cout << "weight of function on: " << fixed_string << " and " << new_string << " is " << degree << "\n";
				degree = 0;
			}
		}
	}

	void full_analytic_structure() {
		std::map<std::string, int> calculating_bits{};

		std::cout << "a1a2a3a4   f\n";

		for (int i = 0; i < (1 << count_degree_of_two); ++i) {
			calculating_bits.insert(std::make_pair(binary_view(i, count_degree_of_two), _boolean[i]));
			std::cout << binary_view(i, count_degree_of_two) << std::setw(8) << _boolean[i];
			std::cout << "\n";
		}
		std::vector<int>  zhegalkin;
		std::cout << "\nHow many query?\n";
		int t;
		std::cin >> t;
		while (t-- > 0) {

			std::cout << "\nInput how many fixed points result do you want to see, remember iteration start with zero\n";
			long long how_many_fixed_points{};
			std::cin >> how_many_fixed_points;
			std::vector<long long> fixed_points(how_many_fixed_points);
			std::cout << "Input fixed points\n";
			for (auto&& x : fixed_points) {
				std::cin >> x;
			}
			std::cout << "\n";
			auto degree = 0;
			for (const auto& [fixed_string, value] : calculating_bits) {
				std::string new_string = "";
				for (std::size_t i = 0; i < fixed_points.size(); ++i) {
					new_string += fixed_string[fixed_points[i]];
				}

				for (const auto& [secondary_string, boolean_value] : calculating_bits) {
					std::string your_string = "";
					for (std::size_t i = 0; i < fixed_points.size(); ++i) {
						your_string += secondary_string[fixed_points[i]];
					}
					if (new_string == your_string) {
						zhegalkin.push_back(boolean_value);
					}
				}

				auto _polynom = polynom_zhegalkin_pascal_(zhegalkin);
				long long k = 1;

				while (_polynom.size() > (1 << k)) {
					++k;
				}

				std::vector<std::string> binary_iteration{};
				for (int i = 0; i < _polynom.size(); ++i) {
					binary_iteration.push_back(binary_view(i, k));
					//std::cout << binary_view(i, k) << std::setw(5) << std::setw(5) << "->" << std::setw(5) << _polynom[i] << "\n";
				}

				std::string polynom_zhegalkin{};

				for (std::size_t i = 0; i < _polynom.size(); ++i) {
					if (_polynom[i] == 1) {
						for (std::size_t j = 0; j < k; ++j) {
							if (binary_iteration[i][j] == '1') {
								polynom_zhegalkin += "x" + std::to_string(j + 1 + fixed_points.size());
							}
						}
						if (i < _polynom.size() - 1) {
							polynom_zhegalkin += "+";
						}
						if (polynom_zhegalkin.size() == 1) {
							polynom_zhegalkin[0] = '1';
							polynom_zhegalkin += "+";
						}
					}
				}

				std::cout << "\nPolynom zhegalkin:\n" << polynom_zhegalkin << " degree: " << degree_of_polynom(polynom_zhegalkin) << "\n";

			}
		}
	}
private:
	decltype(auto) return_sub_arrays(std::vector<int>& a, int mid) {
		std::pair<std::vector<int>, std::vector<int>> pr;
		for (std::size_t i = 0; i < mid; ++i) {
			pr.first.push_back(a[i]);
		}
		for (std::size_t i = mid; i < a.size(); ++i) {
			pr.second.push_back(a[i]);
		}
		return pr;
	}
	void help_fft_for_fourier(std::vector<int>& a, int left, int right, std::vector<int>& result) {
		if (a.size() <= 1) {
			return;
		}

		auto mid = (right - left) / 2;

		auto [first_half, second_half] = return_sub_arrays(a, mid);

		std::vector<int> first_result(first_half.size()), second_result(second_half);

		for (std::size_t i = 0; i < second_half.size(); ++i) {
			first_result[i] = first_half[i] + second_half[i];
		}

		for (std::size_t i = 0; i < second_half.size(); ++i) {
			second_result[i] = first_half[i] - second_half[i];
		}

		help_fft_for_fourier(first_result, 0, mid, result);
		help_fft_for_fourier(second_result, mid, a.size(), result);
		if (first_half.size() == 1) {
			result.push_back(first_result[0]);
			result.push_back(second_result[0]);
		}
	}
	void help_fft_for_zhegalkin(std::vector<int>& a, int left, int right, std::vector<int>& result) {
		if (a.size() <= 1) {
			return;
		}

		auto mid = (right - left) / 2;

		auto [first_half, second_half] = return_sub_arrays(a, mid);

		for (std::size_t i = 0; i < second_half.size(); ++i) {
			second_half[i] ^= first_half[i];
		}

		help_fft_for_zhegalkin(first_half, 0, mid, result);
		help_fft_for_zhegalkin(second_half, mid, a.size(), result);
		if (first_half.size() == 1) {
			result.push_back(first_half[0]);
			result.push_back(second_half[0]);
		}

	}
public:
	void fft_for_fourier_algorithm() {
		std::vector<int> a, b, result;
		for (std::size_t i = _boolean.size() / 2, j = 0; i < _boolean.size(); ++i, ++j) {
			a.push_back(_boolean[i] + _boolean[j]);
		}
		for (std::size_t i = _boolean.size() / 2, j = 0; i < _boolean.size(); ++i, ++j) {
			b.push_back(_boolean[j] - _boolean[i]);
		}
		help_fft_for_fourier(a, 0, a.size(), result);
		help_fft_for_fourier(b, 0, b.size(), result);
		for (const auto& x : result) {
			std::cout << x << ' ';
		}
	}
	void fft() {
		std::vector<int> a(_boolean.size() / 2), b(_boolean.size() / 2), result;
		for (std::size_t i = 0; i < _boolean.size() / 2; ++i) {
			a[i] = _boolean[i];
		}
		for (std::size_t i = _boolean.size() / 2, j = 0; i < _boolean.size(); ++i, ++j) {
			b[j] = _boolean[i] ^ _boolean[j];
		}
		help_fft_for_zhegalkin(a, 0, a.size(), result);
		help_fft_for_zhegalkin(b, 0, b.size(), result);
		for (const auto& x : result) {
			std::cout << x << ' ';
		}
	}
	void walsh_hadamard() {
		std::vector<int> a, b, result, cpy{_boolean.begin(),_boolean.end()};
		std::transform(cpy.begin(), cpy.end(), cpy.begin(), [](auto&& x)
		{
				if (x == 0) {
					return 1;
				}
				else {
					return -1;
				}
		});
		for (std::size_t i = cpy.size() / 2, j = 0; i < cpy.size(); ++i, ++j) {
			a.push_back(cpy[i] + cpy[j]);
		}
		for (std::size_t i = cpy.size() / 2, j = 0; i < cpy.size(); ++i, ++j) {
			b.push_back(cpy[j] - cpy[i]);
		}
		help_fft_for_fourier(a, 0, a.size(), result);
		help_fft_for_fourier(b, 0, b.size(), result);
		for (const auto& x : result) {
			std::cout << x << ' ';
		}

	}

	decltype(auto) polynom_zhegalkin_pascal() noexcept {

		std::vector<int> _polynom;
		std::vector<int> like_boolean{ _boolean.begin(),_boolean.end() };

		while (like_boolean.size() != 0) {
			_polynom.push_back(like_boolean[0]);
			for (int i = 0; i < like_boolean.size() - 1; ++i) {
				like_boolean[i] ^= like_boolean[i + 1];
			}
			like_boolean.pop_back();
		}

		std::cout << "Polynom of zhegalkin and lexgraph boolean function\n\n";

		for (int i = 1; i <= count_degree_of_two; ++i) {
			std::cout << "x" << i;
		}
		std::cout << std::setw(6) << "f";
		std::cout << "\n";

		std::vector<std::string> binary_iteration{};
		for (int i = 0; i < _polynom.size(); ++i) {
			binary_iteration.push_back(binary_view(i, count_degree_of_two));
			std::cout << binary_view(i, count_degree_of_two) << std::setw(5) << std::setw(5) << "->" << std::setw(5) << _polynom[i] << "\n";
		}

		std::string polynom_zhegalkin{};

		for (std::size_t i = 0; i < _polynom.size(); ++i) {
			if (_polynom[i] == 1) {
				for (std::size_t j = 0; j < count_degree_of_two; ++j) {
					if (binary_iteration[i][j] == '1') {
						polynom_zhegalkin += "x" + std::to_string(j + 1);
					}
				}
				if (i < _polynom.size() - 1) {
					polynom_zhegalkin += "+";
				}
				if (polynom_zhegalkin.size() == 1) {
					polynom_zhegalkin[0] = '1';
					polynom_zhegalkin += "+";
				}
			}
		}
		std::cout << "\nPolynom zhegalkin:\n" << polynom_zhegalkin << " degree: " << degree_of_polynom(polynom_zhegalkin) << "\n";
		return polynom_zhegalkin;
	}
	friend std::istream& operator >> (std::istream& in, BooleanFunction& rhs) {
		std::cout << "Input size of boolean function:\n";
		long long n;
		try {
			std::cin >> n;
			if (n < 0) {
				throw std::exception("Unnatural sz of boolean function");
			}
			std::cout << "Input boolean function\n";
			rhs._boolean.resize(n);
			std::string s;
			std::cin >> s;
			if (s.size() > n) {
				throw std::exception("bad length of boolean function");
			}
			for (std::size_t i = 0; i < s.size(); ++i) {
				if (s[i] == '0') {
					rhs._boolean[i] = 0;
				}
				else {
					rhs._boolean[i] = 1;
				}
			}
			while (rhs._boolean.size() > (1 << rhs.count_degree_of_two)) {
				++rhs.count_degree_of_two;
			}
		}
		catch (const std::exception& ex) {
			std::cout << ex.what();
		}
		return in;
	}
	friend std::ostream& operator << (std::ostream& out, const BooleanFunction& rhs) {
		std::cout << "Boolean function is: \n";
		for (const auto& x : rhs._boolean) {
			std::cout << x;
		}
		return out;
	}
	~BooleanFunction() = default;
};

template<typename T>
class LinearReccurentSequence {
private:
	std::vector<T> _sequence;
	std::vector<std::int64_t> _feed_back_polynom;
public:

	LinearReccurentSequence() = default;

	LinearReccurentSequence(std::initializer_list<T> __sequence, std::initializer_list<std::int64_t> __feed_back_polynom) noexcept :
		_sequence(std::move_if_noexcept(__sequence)), _feed_back_polynom(std::move_if_noexcept(__feed_back_polynom))
	{}

	LinearReccurentSequence(std::initializer_list<T> __sequence) noexcept :
		_sequence(std::move_if_noexcept(__sequence))
	{}
	
	auto get_feed_back_polynom() noexcept {
		return _feed_back_polynom;
	}

	auto get_sequence() noexcept {
		return _sequence;
	}

	template<int N>
	void generate_sequence() noexcept {
		auto value = N;
		while (value-- > 0) {
			T sum{};
			for (std::size_t i = 0; i < _feed_back_polynom.size(); ++i) {
				sum += _sequence[_feed_back_polynom[i]];
			}
			std::transform(_feed_back_polynom.begin(), _feed_back_polynom.end(), _feed_back_polynom.begin(), [](auto&& x) {return ++x; });
			_sequence.push_back(sum);
		}
		for (auto&& x : _feed_back_polynom) {
			x %= N;
		}
	}

	void set_feed_back_polynom(std::vector<std::int64_t>& new_polynom) noexcept {
		_feed_back_polynom.assign(new_polynom.begin(), new_polynom.end());
	}

	auto get_size_of_current_lfsr() noexcept {
		return _sequence.size();
	}

	LinearReccurentSequence& operator = (const LinearReccurentSequence& rhs) {
		if (this != &rhs) {
			_sequence.assign(rhs._sequence.begin(), rhs._sequence.end());
		}
		return *this;
	}

	friend std::ostream& operator << (std::ostream& out, const LinearReccurentSequence<T>& rhs) noexcept {
		std::cout << "Linier Reccurent Sequence is:\n";
		for (const auto& x : rhs._sequence) {
			out << x;
		}
		return out;
	}
	friend std::istream& operator >> (std::istream& in, LinearReccurentSequence& rhs) noexcept  {
		std::cout << "register length\n";
		std::size_t n;
		std::cin >> n;
		rhs._sequence.resize(n);
		std::cout << "Input initiall padding of lfsr\n";
		for (auto&& x : rhs._sequence) {
			in >> x;
		}
		std::cout << "Input size of feed back polynom\n";
		std::size_t sz;
		rhs._feed_back_polynom.resize(sz);
		std::cout << "Input removal points\n";
		for (auto&& x : rhs._feed_back_polynom) {
			std::cin >> x;
		}
		return in;
	}

	~LinearReccurentSequence() = default;
};

template<typename T>
class Matrix {
private:
	std::vector<std::vector<T>> _mtrx;
public:
	friend std::istream& operator << (std::istream& in, Matrix& mtrx) {

	}
};

template<typename T>
class VectorSpace {
private:
	std::vector<T> _scalars;
public:
	friend std::ostream& operator << (std::ostream& out, const VectorSpace& rhs) {
		
		return out;
	}
	friend std::istream& operator >> (std::istream& in, VectorSpace& rhs) {

		return in;
	}
};

class Permutation {
private:
	std::vector<std::int64_t> _permutation;
public:
	Permutation() = default;

	Permutation(std::size_t order_of_permutation) noexcept {
		_permutation.resize(order_of_permutation);
		std::iota(_permutation.begin(), _permutation.end(), 1);
	}

	template<int N>
	void set_permutation_order() {
		_permutation.resize(N);
		std::iota(_permutation.begin(), _permutation.end(), 1);
	}

	auto get_next_lexgraph_permutation() {
		return std::next_permutation(_permutation.begin(), _permutation.end());
	}

	auto get_permutation_without_last_element() {
		auto cpy = _permutation;
		cpy.pop_back();
		return cpy;
	}

	auto get_permutation() {
		return _permutation;
	}

	friend std::ostream& operator << (std::ostream& out, const Permutation& rhs) {
		for (const auto& x : rhs._permutation) {
			std::cout << x;
		}
		std::cout << "\n";
		return out;
	}

	~Permutation() = default;
};

#endif


