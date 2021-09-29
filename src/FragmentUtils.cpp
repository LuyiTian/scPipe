#include <string>
#include "FragmentThread.hpp"
#include "FragmentUtils.hpp"


std::string
FragToString(FragmentStruct frag, bool chrom=true, bool start=true, bool end=true, bool bc=true) {
	std::stringstream ss;
	ss  << (chrom ? frag.chromosome : "")
		<< ((start && chrom) ? "|" : "")
		<< (start ? std::to_string(frag.start) : "")
		<< ((end && (start || chrom)) ? "|" : "")
		<< (end ? std::to_string(frag.end) : "")
		<< ((bc && (end || start || chrom)) ? "|" : "")
		<< (bc ? frag.cell_barcode : "");

	return ss.str();
}


// be very careful using this on strings that don't contain the full fragment
FragmentStruct
StringToFrag(std::string s) {
	int l_pos = s.find('|');
	std::string chrom = s.substr(0, l_pos);
	s = s.substr(l_pos + 1);

	l_pos = s.find('|');
	int start = stoi(s.substr(0, l_pos));
	s = s.substr(l_pos + 1);

	l_pos = s.find('|');
	int end = stoi(s.substr(0, l_pos));
	s = s.substr(l_pos + 1);

	return {
		chrom,
		start,
		end,
		s,
		true
	};
}

////////////////////////////////
// Generic typedefs
////////////////////////////////
// create a dictionary from an iterable
template <typename T, typename U>
std::map<U, int> *id_lookup(std::map<std::string, T> &list, std::function<U(T)> f) {
	std::map<U, int> *lookup = new std::map<U, int>;
	int x = 1;
	for (auto el : list) {
		U cur = f(el.second);

		if ((*lookup)[cur] == 0) {
			(*lookup)[cur] = x++;
		}
	}

	return lookup;
}

template <typename T, typename U>
std::map<U, T> invertMap(std::map<T, U> *map) {
	std::map<U, T> i_map;

	for (auto item : map) {
		i_map[item->second] = item->first;
	}

	return i_map;
}
