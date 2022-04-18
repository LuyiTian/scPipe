#include <iostream>
#include <vector>
#include <algorithm>
#include "Interval.h"

#ifndef GENE_H
#define GENE_H

// for position, we dont store chromosome info because we will put genes from the same chromosome
// together in a hashmap
class Gene: public Interval
{
public:
    std::string gene_id;
    std::vector<Interval> exon_vec;

    Gene(std::string id, int st, int en, int snd);
    Gene(std::string id, int snd);
    Gene(std::string id);
    Gene();

    void set_ID(std::string id);
    
    int distance_to_end(Interval it);

    void add_exon(Interval it);

    bool in_exon(const Interval &it);
    bool in_exon(const Interval &it, const bool check_strand);

    // sort exons by starting position
    void sort_exon();
    // flattens exons so that overlapping exons are merged
    void flatten_exon();

    friend std::ostream& operator<< (std::ostream& out, const Gene& obj);
};

#endif
