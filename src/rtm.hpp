#ifndef _RTM_HPP_
#define _RTM_HPP_
#include "table.h"
#include "util.hpp"
typedef uint8_t Bitmap;
typedef __m256i AlignType;
/**
 * @param deg number of neighbors
 * @return size of the bitmap vector (bytes)
 * */
template <typename T>
int get_vector_size(int deg)
{
   const int bitwidth = sizeof(T) * 8;
   return (deg % bitwidth) ? (deg / bitwidth + 1) * sizeof(T): (deg / bitwidth) * sizeof(T);
};
/**
 * @param bitmap triangle intersection vector to be expanded
 * @param out output index place
 * @param vector_size of length of the bitmap
 * */
int expandToIndex(Bitmap *bitmap, int *out, int vector_size);

/**
 * @param bitmap triangle intersection vector to be expanded
 * @param out output vertex id place
 * @param vector_size of length of the bitmap
 * @param id_list start of the adjacancy list
 * */
int expandToID(Bitmap *bitmap, int *out, int vector_size, int *id_list);

/**
 * 
 * @note set intersection of two list, vec_a and vec_b
 * @return number common vertex
 * @return bitvec the position of common vertex in bitmap format relative to vec_a
 * @param vec_a first list
 * @param size_a size of first list
 * @param vec_b second lsit
 * @param size_b size of the second list
 * */
int mark_intersect(int *set_a, int size_a, int *set_b, int size_b, Bitmap *bitvec);
int mark_intersect_simd8x(int *set_a, int size_a, int *set_b, int size_b, Bitmap *bitvec);

/**
 * 
 * @note bitwise and operation of two bitstream
 * @return out the result of two intersection
 * */
void intersect(Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *out, int vector_size);
bool intersect_allzero(Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *out, int vector_size);
void bitwise_not(Bitmap *bitmap_a, Bitmap *out, int vector_size);
void bitwise_and(Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *out, int vector_size);
void bitwise_nand(Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *out, int vector_size);
void bitwise_andn(Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *out, int vector_size);

/**
 * @note out = bitmap_a \ bitmap_b
 * @return number of zeros in out
 * */
int bitwise_andn_count(Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *out, int vector_size);
void fill_with_one(Bitmap * bitmap, int num_one);
void mark_as_one(Bitmap * bitmap, int index);
void mark_as_zero(Bitmap * bitmap, int index);
bool all_zero(Bitmap *bitmap, int vector_size);
/**
 * @param bitmap the coming bitstream
 * @param out place that holds the indices in the bitmap
 * @param start starting point of the bitmap (bits)
 * @param end end place of the bitmap (bits)
 * @return number of indices
 * */
int expandToIndex(Bitmap *bitmap, int *out, int start, int end);
/**
 * @param bitmap the coming bitstream
 * @param out place that holds the vertex ids in the bitmap
 * @param start starting point of the bitmap (bits)
 * @param end end place of the bitmap (bits)
 * @param id_list the id_list corresponds to the bitmap
 * @return number of indices
 * */
int expandToID(Bitmap *bitmap, int *out, int *id_list, int start, int end);

/** 
 * @param bitmap coming bitstream
 * @param vector_size length of the bitstream
 * @return number of 1s in the bitstream
 * */
int count_bitmap(Bitmap *bitmap, int vector_size);
/** 
 * @param bitmap coming bitstream
 * @param vector_size length of the bitstream
 * @param start starting point of the bitmap (bits)
 * @param end end place of the bitmap (bits)
 * @return number of 1s in the bitstream
 * */
int count_bitmap(Bitmap *bitmap, int start, int end);

/**
 * @param bitmap coming bitstream
 * @param vector_size number of bitmap to exam
 * @return true if all the bits are zero, false otherwise
 * */
int find_first_index(Bitmap * bitmap, int vector_size);
int find_last_index(Bitmap * bitmap, int vector_size);

#endif