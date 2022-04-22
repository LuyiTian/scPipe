#ifndef RESIZEARRAY_H
#define RESIZEARRAY_H

class ResizeArray {
        void Expand(void);
        long *arr;
        int size;
    public:
        ResizeArray (int);
        int Increment(int);
        int length(void);
        long operator [] (int);
        int Max(long *);
        void Delete(void);
        
        static void Print(ResizeArray *);
};

#endif