#ifndef _RTM_H_
#define _RTM_H_
#include "util.hpp"
#include "table.h"

namespace RTM_AVX2
{

    typedef uint8_t Bitmap;
    const Bitmap BITMASK = 0xff;
    /**
 * @param deg number of neighbors
 * @return size of the bitmap vector (bytes)
 * */
    int get_vector_size(int deg);

    /**
 * @param bitmap triangle intersection vector to be expanded
 * @param out output index place
 * @param vector_size of length of the bitmap
 * */
    int expandToIndex(uint8_t *bitmap, int *out, int vector_size);

    /**
 * @param bitmap triangle intersection vector to be expanded
 * @param out output vertex id place
 * @param vector_size of length of the bitmap
 * @param id_list start of the adjacancy list
 * */
    int expandToID(uint8_t *bitmap, int *out, int vector_size, int *id_list);

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
    int mark_intersect(int *vec_a, int size_a, int *vec_b, int size_b, uint8_t *bitvec);

    /**
 * 
 * @note bitwise and operation of two bitstream
 * @return out the result of two intersection
 * */
    void intersect(uint8_t *bitmap_a, uint8_t *bitmap_b, uint8_t *out, int vector_size);            // bitwise and operation
    void mask_intersect(uint8_t *bitmap_a, uint8_t *bitmap_b, uint8_t mask, int pos, uint8_t *out); // bitwise and operation with mask applied

    /**
 * @param bitmap the coming bitstream
 * @param out place that holds the indices in the bitmap
 * @param start starting point of the bitmap (bits)
 * @param end end place of the bitmap (bits)
 * @return number of indices
 * */
    int expandToIndex(uint8_t *bitmap, int *out, int start, int end);
    /**
 * @param bitmap the coming bitstream
 * @param out place that holds the vertex ids in the bitmap
 * @param start starting point of the bitmap (bits)
 * @param end end place of the bitmap (bits)
 * @param id_list the id_list corresponds to the bitmap
 * @return number of indices
 * */
    int expandToID(uint8_t *bitmap, int *out, int *id_list, int start, int end);

    /** 
 * @param bitmap coming bitstream
 * @param vector_size length of the bitstream
 * @return number of 1s in the bitstream
 * */
    int count_bitmap(uint8_t *bitmap, int vector_size);

    /** 
 * @param bitmap coming bitstream
 * @param vector_size length of the bitstream
 * @param start starting point of the bitmap (bits)
 * @param end end place of the bitmap (bits)
 * @return number of 1s in the bitstream
 * */
    int count_bitmap(uint8_t *bitmap, int start, int end);

    /**
 * @param bitmap coming bitstream
 * @param vector_size number of bitmap to exam
 * @return true if all the bits are zero, false otherwise
 * */
    bool all_zero(uint8_t *bitmap, int vector_size);
}

namespace RTM_64
{

    typedef uint64_t Bitmap;
    const Bitmap BITMASK = 0xffffffffffffffff;
    /**
 * @param deg number of neighbors
 * @return size of the bitmap vector (bytes)
 * */
    int get_vector_size(int deg);

    /**
 * @param bitmap triangle intersection vector to be expanded
 * @param out output index place
 * @param vector_size of length of the bitmap
 * */
    int expandToIndex(uint64_t *bitmap, int *out, int vector_size);

    /**
 * @param bitmap triangle intersection vector to be expanded
 * @param out output vertex id place
 * @param vector_size of length of the bitmap
 * @param id_list start of the adjacancy list
 * */
    int expandToID(uint64_t *bitmap, int *out, int vector_size, int *id_list);

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
    int mark_intersect(int *vec_a, int size_a, int *vec_b, int size_b, uint64_t *bitvec);
    int mark_intersect_simd8x(int *vec_a, int size_a, int *vec_b, int size_b, uint64_t *bitvec);

    /**
 * 
 * @note bitwise and operation of two bitstream
 * @return out the result of two intersection
 * */
    void intersect(uint64_t *bitmap_a, uint64_t *bitmap_b, uint64_t *out, int vector_size);            // bitwise and operation
    bool intersect_allzero(uint64_t *bitmap_a, uint64_t *bitmap_b, uint64_t *out, int vector_size);            // bitwise and operation
    void mask_intersect(uint64_t *bitmap_a, uint64_t *bitmap_b, uint64_t mask, int pos, uint64_t *out); // bitwise and operation with mask applied
    void mask_intersect_simd4x(uint64_t *bitmap_a, uint64_t *bitmap_b, uint64_t mask, int pos, uint64_t *out); // bitwise and operation with mask applied

    /**
 * @param bitmap the coming bitstream
 * @param out place that holds the indices in the bitmap
 * @param start starting point of the bitmap (bits)
 * @param end end place of the bitmap (bits)
 * @return number of indices
 * */
    int expandToIndex(uint64_t *bitmap, int *out, int start, int end);
    /**
 * @param bitmap the coming bitstream
 * @param out place that holds the vertex ids in the bitmap
 * @param start starting point of the bitmap (bits)
 * @param end end place of the bitmap (bits)
 * @param id_list the id_list corresponds to the bitmap
 * @return number of indices
 * */
    int expandToID(uint64_t *bitmap, int *out, int *id_list, int start, int end);

    /** 
 * @param bitmap coming bitstream
 * @param vector_size length of the bitstream
 * @return number of 1s in the bitstream
 * */
    int count_bitmap(uint64_t *bitmap, int vector_size);

    /** 
 * @param bitmap coming bitstream
 * @param vector_size length of the bitstream
 * @param start starting point of the bitmap (bits)
 * @param end end place of the bitmap (bits)
 * @return number of 1s in the bitstream
 * */
    int count_bitmap(uint64_t *bitmap, int start, int end);

    /**
 * @param bitmap coming bitstream
 * @param vector_size number of bitmap to exam
 * @return true if all the bits are zero, false otherwise
 * */
    bool all_zero(uint64_t *bitmap, int vector_size);
}
#endif