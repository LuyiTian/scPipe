#pragma once
// the general interval class
class Interval
{
public:
    int st;
    int en;
    int snd; // strand, -1 for negative, 1 for positive, 0 for unknown
    Interval(int s, int e);
    Interval(int s, int e, int sd);

    // check if interval (st1, en1) overlaps with (st2, en2)
    // return -1 if (st2, en2) on the left side of (st1, en1)
    // return 1 if (st2, en2) on the right side of (st1, en1)
    // return 0 if overlap
    int overlap(int st2, int en2);

    friend inline bool operator == (const Interval &L, const Interval &R);
    friend inline bool operator != (const Interval &L, const Interval &R);
    friend inline bool operator <  (const Interval &L, const Interval &R);
    friend inline bool operator >  (const Interval &L, const Interval &R);
    friend inline bool operator <= (const Interval &L, const Interval &R);
    friend inline bool operator >= (const Interval &L, const Interval &R);
};