#ifndef FRAGMENT_UTILS
#define FRAGMENT_UTILS

#include <string>
#include <map>
#include <functional>
#include "Fragments.h"

std::string
FragToString(FragmentStruct, bool, bool, bool, bool);

FragmentStruct
StringToFrag(std::string);

std::map<std::string, int> *id_lookup(std::map<std::string, FragmentStruct> &, std::function<std::string(FragmentStruct)> );

std::map<int, std::string> invertMap(std::map<std::string, int> *);

std::string FragMapToString (FragmentMap &);
#endif
