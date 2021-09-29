#include <string>
#include <vector>
#include "FragmentThread.hpp"
#include "FragmentUtils.hpp"

////////////
// This contains all static utility functions for the FragmentThread clas
//////////////


// Done Not Tested
/// create a map where key is subsetted fragment coords (no start or no stop)
/// value is the list of full fragments that share the coordinates in the key
/// only entries where >1 full fragments share the same coordinate are retained.
/// @param fragments is the FragmentMap containing all fragments
/// @param start indicates if start should be retained. If False, end is retained
std::map<std::string, std::vector<FragmentStruct>> 
FragmentThread::createPositionLookup(FragmentMap &fragments, bool start){
	// basic function to convert a fragment to a string including either
	// start or end position
	std::function<std::string(FragmentStruct &, bool, bool, bool, bool)> convertToString = 
		[&start](FragmentStruct &f, bool _1, bool _2, bool _3, bool _4)->std::string {
			return FragToString(f, true, start, !start, true);
		};

	std::map<std::string, int> counts = 
		FragmentThread::CounterMapFragment(fragments,convertToString);

	std::map<std::string, std::vector<FragmentStruct>> starts;

	// iterate over all the fragments
	// for each string version which has >1 counts
	// add it to the forward list at that location.
	for (auto frag : fragments) {
		std::string key = convertToString(frag.second, true, true, true, true);
		if (counts[key] > 1) {
			starts[key].push_back(frag.second);
		}
	}

	return starts;
}

// Done Not Tested
/// General purpose function for generating a Counter map such that
/// each value is the number of times the given key appeared in the given map
/// Keys are generated using the provided function f.
std::map<std::string, int>
FragmentThread::CounterMapFragment(FragmentMap &fragments, std::function<std::string(FragmentStruct&, bool, bool, bool, bool)> f) {
	std::map<std::string, int> counts;
	for (auto frag : fragments) {
		counts[f(frag.second, true, true, true, true)]++;
	}
	return counts;
}