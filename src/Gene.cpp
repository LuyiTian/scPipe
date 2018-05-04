#include "Gene.h"

using std::string;

Gene::Gene(string id, int st, int en, int snd) : Interval(st, en, snd), gene_id(id) {}
Gene::Gene(string id, int snd)                 : Interval(-1, -1, snd), gene_id(id) {}
Gene::Gene(string id)                          : Interval(-1, -1, 0), gene_id(id) {}
Gene::Gene()                                   : Interval(-1, -1, 0), gene_id("") {}

void Gene::set_ID(string id)
{
    gene_id = id;
}

void Gene::add_exon(Interval it)
{
        exon_vec.push_back(it);
        // expand the gene interval to include the new exon
        if (st > it.st || st < 0)
        {
            st = it.st;
        }
        if (en < it.en || en < 0)
        {
            en = it.en;
        }
        if (snd == 0)
        {
            snd = it.snd;
        }
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
    std::sort(exon_vec.begin(), exon_vec.end(),
        [] (const Interval &a, const Interval &b) { return a.st < b.st; }
    );
}

void Gene::flatten_exon() {
    std::vector<Interval> merged_exons;
    merged_exons.reserve(exon_vec.size());

    merged_exons.push_back(exon_vec[0]);

    for (auto i = 1; i < exon_vec.size(); i++)
    {
        const auto exon = exon_vec[i];
        const auto &last_merged_exon = merged_exons.back();
        if (exon.st > last_merged_exon.en) {
            // if new exon does not overlap last merged exon
            merged_exons.push_back(exon);
        }
        else if (exon.en > last_merged_exon.en) {
            // if new exon does overlap last merged exon and ends later
            auto temp_exon = exon;
            temp_exon.st = last_merged_exon.st;
            merged_exons.back() = temp_exon;
        }
    }

    exon_vec = merged_exons;
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
