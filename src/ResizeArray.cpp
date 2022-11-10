#include <iostream>
#include <Rcpp.h>
#include "ResizeArray.h"

ResizeArray::ResizeArray(int len) {
    size = len;
    arr = new long[size];
    for (int i = 0; i < size; i++) {
        arr[i] = 0;
    }
}
void ResizeArray::Expand(void) {
    long *new_array = new long[size * 2];
    for (int i = 0; i < (size * 2); i++) {
        if (i < size) {
            new_array[i] = arr[i];
        } else {
            new_array[i] = 0;
        }
    }
    size = size * 2;
    delete [] arr;
    arr = new_array;
}
int ResizeArray::Increment(int position) {
    while (position >= size) {
        Expand();
    }
    arr[position]++;
    //Rprintf("\tincrement %d, now %ld\n", position, arr[position]);
    return size;
}
int ResizeArray::length(void) {    return size;
}
long ResizeArray::operator[] (int index) {
    if (index >= 0 && index < size) {
        return arr[index];
    } else {
        throw std::out_of_range("Index must be in range of resizable array size\n");
    }
}
int ResizeArray::Max (long *max_value) {
    long curr_max = 0;
    int position = -1;
    for (int i = 0; i < size; i++) {
        if (arr[i] > curr_max) {
            curr_max = arr[i];
            position = i;
        }
    }
    *max_value = curr_max;
    return position;
}
void ResizeArray::Delete (void) {
    delete [] arr;
}

void ResizeArray::Print(ResizeArray *positions) {
    for (int i = 0; i < positions->length(); i++) {
        Rprintf("Pos %d, val: %ld\t", i, (*positions)[i]);
        if (i % 5 == 4) {
            Rprintf("\n");
        }
    }
}
