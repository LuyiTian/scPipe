#include "Gene.h"

using std::string;

Gene::Gene(string id, int st, int en, int snd): Interval(st, en, snd), gene_id(id) {}
Gene::Gene(string id, int snd): Interval(-1, -1, snd), gene_id(id) {}
Gene::Gene(): Interval(-1, -1, 0), gene_id("") {}

void Gene::set_ID(string id)
{
    gene_id = id;
}

void Gene::add_exon(Interval it)
{
        exon_vec.push_back(it);
        // expand the gene interval to include the new exon
        if (st>it.st || st<0)
        {
            st = it.st;
        }
        if (en < it.en || en<0)
        {
            en = it.en;
        }
        if (snd == 0)
        {
            snd = it.snd;
        }
}

int Gene::distance_to_end(Interval it)
{
    int distance = 0;
    int tmp_en = 0;
    auto iter = std::lower_bound(exon_vec.begin(), exon_vec.end(), it);
    if (snd == 1)
    {
        distance += iter->en - ((iter->st)>it.st?(iter->st):it.st);
        tmp_en = iter->en;
        for (auto i = iter+1; i != exon_vec.end(); ++i)
        {
            if (tmp_en < i->st)
            {
                distance += i->en - i->st;
                tmp_en = i->en;
            }

        }

    }
    else if (snd == -1)
    {
        for (auto i = exon_vec.begin(); i != iter; ++i)
        {
            if (tmp_en < i->st)
            {
                distance += i->en - i->st;
                tmp_en = i->en;
            }
        }
        if (tmp_en < iter->st)
        {
            distance += ((iter->en)<it.en?(iter->en):it.en) - iter->st;
        }
    }

    return distance;
}

bool Gene::in_exon(const Interval &it)
{
    auto search_result = std::find(exon_vec.begin(), exon_vec.end(), it);
    return search_result != exon_vec.end();
}

bool Gene::in_exon(const Interval &it, const bool check_strand)
{
    if (check_strand && (it.snd*snd == -1))
    {
        return false;
    }
    else
    {
        auto search_result = std::find(exon_vec.begin(), exon_vec.end(), it);
        return search_result != exon_vec.end();
    }
}

void Gene::sort_exon()
{
    std::sort(exon_vec.begin(), exon_vec.end());
}

std::ostream& operator<< (std::ostream& out, const Gene& obj)
{
    out << "Gene ID:   " << obj.gene_id  << "\n";
    out << "\t" << "start/end:   " << obj.st  << "/" << obj.en << "\n";
    out << "\t" << "strand:   " << obj.snd  << "\n";
    out << "\t" << "number of exons:   " << obj.exon_vec.size()  << "\n";
    for (int i = 0; i < obj.exon_vec.size(); ++i)
    {
        out << "\t" << "exon[" << i+1 << "]: (" << obj.exon_vec[i].st << ", " << obj.exon_vec[i].en << ")" << "\n";
    }
    return out;
}
