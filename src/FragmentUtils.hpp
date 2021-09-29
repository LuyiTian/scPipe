#include <functional>
#include <map>

#include "FragmentThread.hpp"


#ifndef FRAGMENT_UTILS
#define FRAGMENT_UTILS

std::string FragToString(FragmentStruct, bool, bool, bool, bool);
FragmentStruct StringToFrag(std::string);

template <typename T, typename U>
std::map<U, int> *id_lookup(std::map<std::string, T> &, std::function<U(T)>);

template <typename T, typename U>
std::map<U, T> invertMap(std::map<T, U> *);
#endif
