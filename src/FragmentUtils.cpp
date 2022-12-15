#include "Fragments.h"

#include <string>
#include <map>
#include <functional>
#include <sstream>

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
		true,
		0
	};
}

// create a dictionary from an iterable
// where each unique item is key, value is the order it was inserted into the map
std::map<std::string, int> *id_lookup(FragmentMap &list, std::function<std::string(FragmentStruct)> f) {
	std::map<std::string, int> *lookup = new std::map<std::string, int>;
	int x = 1;
	for (FragmentMap::iterator el = list.begin(); el != list.end(); el++) {
		std::string cur = f(el->second);

		if ((*lookup)[cur] == 0) {
			(*lookup)[cur] = x++;
		}
	}

	for (auto it = lookup->begin(); it != lookup->end(); it++) {
		it->second--;
	}

	return lookup;
}

std::map<int, std::string> invertMap(std::map<std::string, int> *map) {
	std::map<int, std::string> i_map;

	for (std::map<std::string, int>::iterator item = map->begin(); item != map->end(); item++) {
		i_map[item->second] = item->first;
	}

	return i_map;
}

std::string FragMapToString(FragmentMap &fm) {
	std::stringstream ss;

	for (auto it : fm) {
		ss << it.first << ": " << FragToString(it.second, true, true, true, true) << "\n";
	}

	return ss.str();
}