#ifndef INTERVAL_H
#define INTERVAL_H

// the general interval class
class Interval
{
public:
    int st;
    int en;
    int snd; // strand, -1 for negative, 1 for positive, 0 for unknown
    
    Interval(int start, int end);
    Interval(int start, int end, int strand);

    // check if this interval overlaps with (start, end)
    // return -1 if (start, end) on the left side of this interval
    // return 1 if (start, end) on the right side of this interval
    // return 0 if overlap
    int overlap(int start, int end);

    friend inline bool operator < (const Interval &L, const Interval &R) {
        return L.en < R.st;
    };
    friend inline bool operator > (const Interval &L, const Interval &R) {
        return L.st > R.en;
    };
    friend inline bool operator == (const Interval &L, const Interval &R) {
        return !(L<R) && !(L>R);
    };
};

#endif
