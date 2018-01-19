#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include "Interval.h"

// for position, we dont store chromosome info because we will put genes from the same chromosome
// together in a hashmap
class Gene: public Interval
{
public:
    std::string gene_id;
    std::vector<Interval> exon_vec;

    Gene(std::string id, int st, int en, int snd);
    Gene(std::string id, int snd);
    Gene();

    void set_ID(std::string id);

    void add_exon(Interval it);

    bool in_exon(Interval it);
    bool in_exon(Interval it, bool check_strand);

    int distance_to_end(Interval it);
    void sort_exon();

    friend std::ostream& operator<< (std::ostream& out, const Gene& obj);
};