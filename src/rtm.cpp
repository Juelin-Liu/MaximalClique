#include "rtm.hpp"
#include "rtm_table.h"
// #define __USING_AVX2__
#define _CYCLIC_SHIFT1_ 57  //_MM_SHUFFLE(0,3,2,1); //rotating right
#define _CYCLIC_SHIFT2_ 147 //_MM_SHUFFLE(2,1,0,3); //rotating left
#define _CYCLIC_SHIFT3_ 78  //_MM_SHUFFLE(1,0,3,2);
#define _MASK_SHIFT_ 27     //_MM_SHUFFLE(0,1,2,3);
const Bitmap BITMASK = 0xff;
const __m256i add8 = _mm256_set1_epi32(8);
const __m256i add64 = _mm256_set1_epi32(64);
const int base[8] = {0, 1, 2, 3, 4, 5, 6, 7};

int expand_avx2(Bitmap *bitmap, int *out, int vector_size)
{
   int *initout = out;
   __m256i baseVec = _mm256_set1_epi32(-1);
   uint8_t mask;
   int i = 0;
   while (i < vector_size)
   {
      if (*(uint64_t *)(bitmap + i) == 0)
      {
         baseVec = _mm256_add_epi32(baseVec, add64);
         i += 8;
      }
      else
      {
         mask = bitmap[i];
         __m256i vecA = _mm256_load_si256((const __m256i *)vecDecodeTable[mask]);
         vecA = _mm256_add_epi32(baseVec, vecA);
         _mm256_storeu_si256((__m256i *)out, vecA);
         out += lengthTable[mask];
         baseVec = _mm256_add_epi32(baseVec, add8);
         i++;
      }
   }
   return out - initout;
};

int expand_avx2_compress(Bitmap *bitmap, int *out, int vector_size)
{
   int *init_out = out;
   const __m256i zero_vec = _mm256_set1_epi8(char(0));
   for (int i = 0; i < vector_size; i += 32)
   {
      __m256i cmp_mask = _mm256_cmpeq_epi8(*(__m256i *)&bitmap[i], zero_vec);
      uint32_t bitset = ~_mm256_movemask_epi8(cmp_mask);
      while (bitset != 0)
      {
         uint32_t t = bitset & -bitset;
         int r = __builtin_ctz(bitset);
         __m256i baseVec = _mm256_set1_epi32((i + r) * 8 - 1);
         uint8_t mask = bitmap[r + i];
         __m256i vecA = _mm256_load_si256((const __m256i *)vecDecodeTable[mask]);
         vecA = _mm256_add_epi32(baseVec, vecA);
         _mm256_storeu_si256((__m256i *)out, vecA);
         out += lengthTable[mask];
         bitset ^= t;
      };
   }
   return out - init_out;
};

int maskz_expand_avx2_compress(Bitmap *bitmask, Bitmap *bitmap, Bitmap *buffer, int *out, int vector_size)
{
   int *init_out = out;
   const __m256i zero_vec = _mm256_set1_epi8(char(0));
   for (int i = 0; i < vector_size; i += 32)
   {
      *(__m256i *)buffer = _mm256_andnot_si256(*(__m256i *)&bitmask[i], *(__m256i *)&bitmap[i]);
      __m256i cmp_mask = _mm256_cmpeq_epi8(*(__m256i *)buffer, zero_vec);
      uint32_t bitset = ~_mm256_movemask_epi8(cmp_mask);
      while (bitset != 0)
      {
         uint32_t t = bitset & -bitset;
         int r = __builtin_ctz(bitset);
         __m256i baseVec = _mm256_set1_epi32((i + r) * 8 - 1); // large overhead
         uint8_t mask = buffer[r];
         __m256i vecA = _mm256_load_si256((const __m256i *)vecDecodeTable[mask]);
         vecA = _mm256_add_epi32(baseVec, vecA);
         _mm256_storeu_si256((__m256i *)out, vecA); // large overhead
         out += lengthTable[mask];
         bitset ^= t;
      };
   }
   return out - init_out;
};

int maskzor_expand_avx2_compress(Bitmap *bitmask, Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *buffer, int *out, int vector_size)
{
   int *init_out = out;
   const __m256i zero_vec = _mm256_set1_epi8(char(0));
   for (int i = 0; i < vector_size; i += 32)
   {
      *(__m256i *)buffer = _mm256_or_si256(*(__m256i *)&bitmap_a[i], *(__m256i *)&bitmap_b[i]);
      *(__m256i *)buffer = _mm256_andnot_si256(*(__m256i *)&bitmask[i], *(__m256i *)buffer);
      __m256i cmp_mask = _mm256_cmpeq_epi8(*(__m256i *)buffer, zero_vec);
      uint32_t bitset = ~_mm256_movemask_epi8(cmp_mask);
      while (bitset != 0)
      {
         uint32_t t = bitset & -bitset;
         int r = __builtin_ctz(bitset);
         __m256i baseVec = _mm256_set1_epi32((i + r) * 8 - 1); // large overhead
         uint8_t mask = buffer[r];
         __m256i vecA = _mm256_load_si256((const __m256i *)vecDecodeTable[mask]);
         vecA = _mm256_add_epi32(baseVec, vecA);
         _mm256_storeu_si256((__m256i *)out, vecA); // large overhead
         out += lengthTable[mask];
         bitset ^= t;
      };
   }
   return out - init_out;
};

int expand_ctz(Bitmap *bitmap, int *out, int vector_size)
{
   int pos = 0;
   for (int i = 0; i < vector_size; i += 8)
   {
      uint64_t bitset = *(uint64_t *)(bitmap + i);
      while (bitset != 0)
      {
         uint64_t t = bitset & -bitset;
         int r = __builtin_ctzll(bitset);
         out[pos++] = i * 8 + r;
         bitset ^= t;
      }
   }
   return pos;
};

int maskz_expand_ctz(Bitmap *bitmask, Bitmap *bitmap_a, int *out, int vector_size)
{
   int pos = 0;
   for (int i = 0; i < vector_size; i += 8)
   {
      uint64_t bitset = *(uint64_t *)(bitmap_a + i) & ~*(uint64_t *)(bitmask + i);
      while (bitset != 0)
      {
         uint64_t t = bitset & -bitset;
         int r = __builtin_ctzll(bitset);
         out[pos++] = i * 8 + r;
         bitset ^= t;
      }
   }
   return pos;
};

int maskzor_expand_ctz(Bitmap *bitmask, Bitmap *bitmap_a, Bitmap *bitmap_b, int *out, int vector_size)
{
   int pos = 0;
   for (int i = 0; i < vector_size; i += 8)
   {
      uint64_t bitset = (*(uint64_t *)(bitmap_a + i) | *(uint64_t *)(bitmap_b + i)) & (~(*(uint64_t *)(bitmask + i)));
      while (bitset != 0)
      {
         uint64_t t = bitset & -bitset;
         int r = __builtin_ctzll(bitset);
         out[pos++] = i * 8 + r;
         bitset ^= t;
      }
   }
   return pos;
};

int expandToID(Bitmap *bitmap, int *out, int vector_size, int *id_list)
{
   int p = expand_avx2(bitmap, out, vector_size);
   for (int i = 0; i < p; i++)
   {
      out[i] = id_list[out[i]];
   }
   return p;
};

int mark_intersect(int *set_a, int size_a, int *set_b, int size_b, Bitmap *bitvec)
{
   int i = 0, j = 0, cnt = 0;

   while (i < size_a && j < size_b)
   {
      if (set_a[i] == set_b[j])
      {
         bitvec[i / 8] |= 1 << (i % 8);
         i++;
         j++;
         cnt++;
      }
      else if (set_a[i] < set_b[j])
      {
         i++;
      }
      else
      {
         j++;
      }
   }
   return cnt;
};

int mark_intersect_simd8x(int *set_a, int size_a, int *set_b, int size_b, Bitmap *bitvec)
{
   int i = 0, j = 0, k = 0, cnt = 0;
   int qs_a = size_a - (size_a & 7);
   int qs_b = size_b - (size_b & 7);
   //   const __m256i reverse_mask = _mm256_loadu_si256((__m256i *)reverese);

   while (i < qs_a && j < qs_b)
   {
      __m256i v_a = _mm256_lddqu_si256((__m256i *)(set_a + i));
      __m256i v_b = _mm256_lddqu_si256((__m256i *)(set_b + j));
      k = i / 8;
      int a_max = set_a[i + 7];
      int b_max = set_b[j + 7];
      if (a_max == b_max)
      {
         i += 8;
         j += 8;
         _mm_prefetch((const char *)set_a + i, _MM_HINT_NTA);
         _mm_prefetch((const char *)set_b + j, _MM_HINT_NTA);
      }
      else if (a_max > b_max)
      {
         j += 8;
         _mm_prefetch((const char *)set_b + j, _MM_HINT_NTA);
      }
      else
      {
         i += 8;
         _mm_prefetch((const char *)set_a + i, _MM_HINT_NTA);
      }

      __m256i cmp_mask0 = _mm256_cmpeq_epi32(v_a, v_b);
      __m256 rot1 = _mm256_permute_ps((__m256)v_b, _CYCLIC_SHIFT1_);
      __m256i cmp_mask1 = _mm256_cmpeq_epi32(v_a, (__m256i)rot1);
      __m256 rot2 = _mm256_permute_ps((__m256)v_b, _CYCLIC_SHIFT2_);
      __m256i cmp_mask2 = _mm256_cmpeq_epi32(v_a, (__m256i)rot2);
      __m256 rot3 = _mm256_permute_ps((__m256)v_b, _CYCLIC_SHIFT3_);
      __m256i cmp_mask3 = _mm256_cmpeq_epi32(v_a, (__m256i)rot3);

      __m256 rot4 = _mm256_permute2f128_ps((__m256)v_b, (__m256)v_b, 1);
      __m256i cmp_mask4 = _mm256_cmpeq_epi32(v_a, (__m256i)rot4);
      __m256 rot5 = _mm256_permute_ps(rot4, _CYCLIC_SHIFT1_);
      __m256i cmp_mask5 = _mm256_cmpeq_epi32(v_a, (__m256i)rot5);
      __m256 rot6 = _mm256_permute_ps(rot4, _CYCLIC_SHIFT2_);
      __m256i cmp_mask6 = _mm256_cmpeq_epi32(v_a, (__m256i)rot6);
      __m256 rot7 = _mm256_permute_ps(rot4, _CYCLIC_SHIFT3_);
      __m256i cmp_mask7 = _mm256_cmpeq_epi32(v_a, (__m256i)rot7);

      __m256i cmp_mask = _mm256_or_si256(
          _mm256_or_si256(
              _mm256_or_si256(cmp_mask0, cmp_mask1),
              _mm256_or_si256(cmp_mask2, cmp_mask3)),
          _mm256_or_si256(
              _mm256_or_si256(cmp_mask4, cmp_mask5),
              _mm256_or_si256(cmp_mask6, cmp_mask7)));
      // __m256 shifted_mask = _mm256_permutevar8x32_ps((__m256)cmp_mask, reverse_mask);
      uint8_t mask = _mm256_movemask_ps((__m256)cmp_mask);
      bitvec[k] |= mask;
      cnt += lengthTable[mask];
   }

   while (i < size_a && j < size_b)
   {
      if (set_a[i] == set_b[j])
      {
         bitvec[i / 8] |= 1 << (i % 8);
         i++;
         j++;
         cnt++;
      }
      else if (set_a[i] < set_b[j])
      {
         i++;
      }
      else
      {
         j++;
      }
   }
   return cnt;
};

int expand_avx2(Bitmap *bitmap, int *out, int start, int end)
{
   if (start == end)
      return 0;
   int *initout = out;
   uint8_t mask;
   int start_vec = start / 8;
   int end_vec = get_vector_size<Bitmap>(end);
   __m256i baseVec = _mm256_set1_epi32(start_vec * 8 - 1);
   mask = bitmap[start_vec];
   mask = mask & (BITMASK << (start % 8));
   __m256i vecA = _mm256_load_si256((const __m256i *)vecDecodeTable[mask]);
   vecA = _mm256_add_epi32(baseVec, vecA);
   _mm256_storeu_si256((__m256i *)out, vecA);
   out += lengthTable[mask];
   baseVec = _mm256_add_epi32(baseVec, add8);
   int i = start_vec + 1;
   while (i < end_vec)
   {
      if (*(uint64_t *)(bitmap + i) == 0)
      {
         baseVec = _mm256_add_epi32(baseVec, add64);
         i += 8;
      }
      else
      {
         mask = bitmap[i];
         __m256i vecA = _mm256_load_si256((const __m256i *)vecDecodeTable[mask]);
         vecA = _mm256_add_epi32(baseVec, vecA);
         _mm256_storeu_si256((__m256i *)out, vecA);
         out += lengthTable[mask];
         baseVec = _mm256_add_epi32(baseVec, add8);
         i++;
      }
   }
   return out - initout;
};

int expandToID(Bitmap *bitmap, int *out, int *id_list, int start, int end)
{
   int p = expand_avx2(bitmap, out, start, end);
   for (int i = 0; i < p; i++)
   {
      out[i] = id_list[out[i]];
   }
   return p;
};

int count_bitmap(Bitmap *bitmap, int vector_size)
{
   int p = 0;
   for (int i = 0; i < vector_size - 8; i += 8)
   {
      p += __builtin_popcountll(*(uint64_t *)&bitmap[i]);
   }
   for (int i = vector_size - vector_size % 8; i < vector_size; i++)
   {
      p += lengthTable[bitmap[i]];
   }
   return p;
};

int count_bitmap(Bitmap *bitmap, int start, int end)
{
   if (start == end)
      return 0;
   int vector_end = get_vector_size<Bitmap>(end);
   int vector_start = start / 8;
   auto bitset = bitmap[vector_start] & (0xff << (start % 8));
   int vector_size = vector_end - ++vector_start;
   if (vector_size > 0)
   {
      return lengthTable[bitset] + count_bitmap(bitmap + vector_start, vector_size);
   }
   else
   {
      return lengthTable[bitset];
   }
};

bool all_zero(Bitmap *bitmap, int vector_size)
{
   #ifdef __USING_AVX2__
      const __m256i zero_vec = _mm256_set1_epi8(char(0));
      for (int i = 0; i < vector_size; i += 32)
      {
         if (~_mm256_movemask_epi8(_mm256_cmpeq_epi8(*(__m256i *)&bitmap[i], zero_vec)))
         {
            return false;
         }
      }
      return true;
   #else
   for (int i = 0; i < vector_size; i += 8)
   {
      if (*(uint64_t *)(bitmap + i))
      {
         return false;
      }
   }
   return true;
   #endif
};

void bitwise_nand(Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *out, int vector_size)
{
   for (int i = 0; i < vector_size; i++)
   {
      out[i] = ~(bitmap_a[i] & bitmap_b[i]);
   }
};

void bitwise_andn(Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *out, int vector_size)
{
#ifdef __USING_AVX2__
   for (int i = 0; i < vector_size; i += 32)
   {
      // order of bitmap_a and bitmap_b matters!
      _mm256_store_si256((__m256i *)&out[i], _mm256_andnot_si256(*(__m256i *)&bitmap_b[i], *(__m256i *)&bitmap_a[i]));
   }
#else
   for (int i = 0; i < vector_size; i++)
   {
      out[i] = bitmap_a[i] & ~bitmap_b[i];
   }
#endif
};

// large overhead
int bitwise_andn_count(Bitmap *bitmap_a, Bitmap *bitmap_b, int vector_size)
{
#ifdef __USING_AVX2__
   int cnt = 0;
   for (int i = 0; i < vector_size; i += 32)
   {
      // order of bitmap_a and bitmap_b matters!
      __m256i and_not_mask = _mm256_andnot_si256(*(__m256i *)&bitmap_b[i], *(__m256i *)&bitmap_a[i]);
      cnt += __builtin_popcountll(_mm256_extract_epi64(and_not_mask, 0)) + __builtin_popcountll(_mm256_extract_epi64(and_not_mask, 1)) + __builtin_popcountll(_mm256_extract_epi64(and_not_mask, 2)) + __builtin_popcountll(_mm256_extract_epi64(and_not_mask, 3));
   }
   return cnt;
#else
   int cnt = 0;
   for (int i = 0; i < vector_size; i += 8)
   {
      uint64_t bitset = *(uint64_t *)(bitmap_a + i) & ~ *(uint64_t *)(bitmap_b + i);
      cnt += __builtin_popcountll(bitset);
   }
#endif
};
int bitwise_andn_count(Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *out, int vector_size)
{
#ifdef __USING_AVX2__
   int cnt = 0;
   for (int i = 0; i < vector_size; i += 32)
   {
      // order of bitmap_a and bitmap_b matters!
      __m256i and_not_mask = _mm256_andnot_si256(*(__m256i *)&bitmap_b[i], *(__m256i *)&bitmap_a[i]);
      _mm256_store_si256((__m256i *)(out + i), and_not_mask);
      cnt += __builtin_popcountll(_mm256_extract_epi64(and_not_mask, 0)) + __builtin_popcountll(_mm256_extract_epi64(and_not_mask, 1)) + __builtin_popcountll(_mm256_extract_epi64(and_not_mask, 2)) + __builtin_popcountll(_mm256_extract_epi64(and_not_mask, 3));
   }
   return cnt;
#else
   int cnt = 0;
   for (int i = 0; i < vector_size; i += 8)
   {
      uint64_t bitset = *(uint64_t *)(bitmap_a + i) & ~*(uint64_t *)(bitmap_b + i);
      *(uint64_t *)(out + i) = bitset;
      cnt += __builtin_popcountll(*(uint64_t *)bitset);
   }
   return cnt;
#endif
};

int bitwise_and_count(Bitmap *bitmap_a, Bitmap *bitmap_b, int vector_size)
{
#ifdef __USING_AVX2__
   int cnt = 0;
   for (int i = 0; i < vector_size; i += 32)
   {
      __m256i andmask = _mm256_and_si256(*(__m256i *)&bitmap_b[i], *(__m256i *)&bitmap_a[i]);
      cnt += __builtin_popcountll(_mm256_extract_epi64(andmask, 0)) + __builtin_popcountll(_mm256_extract_epi64(andmask, 1)) + __builtin_popcountll(_mm256_extract_epi64(andmask, 2)) + __builtin_popcountll(_mm256_extract_epi64(andmask, 3));
   }
   return cnt;
#else
   int cnt = 0;
   for (int i = 0; i < vector_size; i += 8)
   {
      cnt += __builtin_popcountll(*(uint64_t *)(bitmap_a + i) & *(uint64_t *)(bitmap_b + i));
   }
   return cnt;
#endif
};

int bitwise_and_count(Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *out, int vector_size)
{
#ifdef __USING_AVX2__
   int cnt = 0;
   for (int i = 0; i < vector_size; i += 32)
   {
      __m256i andmask = _mm256_and_si256(*(__m256i *)&bitmap_b[i], *(__m256i *)&bitmap_a[i]);
      _mm256_store_si256((__m256i *)(out+i), andmask);
      cnt += __builtin_popcountll(_mm256_extract_epi64(andmask, 0)) + __builtin_popcountll(_mm256_extract_epi64(andmask, 1)) + __builtin_popcountll(_mm256_extract_epi64(andmask, 2)) + __builtin_popcountll(_mm256_extract_epi64(andmask, 3));
   }
   return cnt;
#else
   int cnt = 0;
   for (int i = 0; i < vector_size; i += 8)
   {
      uint64_t bitset = *(uint64_t *)(bitmap_a + i) & *(uint64_t *)(bitmap_b + i);
      *(uint64_t *)(out+i) = bitset;
      cnt += __builtin_popcountll(bitset);
   }
   return cnt;
#endif
};

void bitwise_and(Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *out, int vector_size)
{
#ifdef __USING_AVX2__
   for (int i = 0; i < vector_size; i += 32)
   {
      _mm256_store_si256((__m256i *)&out[i], _mm256_and_si256(*(__m256i *)&bitmap_a[i], *(__m256i *)&bitmap_b[i]));
   }
#else
   for (int i = 0; i < vector_size; i++)
   {
      *(uint64_t *)&out[i] = *(uint64_t *)&bitmap_a[i] & *(uint64_t *)&bitmap_b[i];
   }
#endif
};
void double_bitwise_and(Bitmap *bitmask, Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *a_out, Bitmap *b_out, int vector_size)
{
#ifdef __USING_AVX2__
   for (int i = 0; i < vector_size; i += 32)
   {
      _mm256_store_si256((__m256i *)&a_out[i], _mm256_and_si256(*(__m256i *)&bitmask[i], *(__m256i *)&bitmap_a[i]));
      ;
      _mm256_store_si256((__m256i *)&b_out[i], _mm256_and_si256(*(__m256i *)&bitmask[i], *(__m256i *)&bitmap_b[i]));
      ;
   }
#else
   for (int i = 0; i < vector_size; i++)
   {
      *(uint64_t *)&a_out[i] = *(uint64_t *)&bitmap_a[i] & *(uint64_t *)&bitmask[i];
      *(uint64_t *)&b_out[i] = *(uint64_t *)&bitmap_b[i] & *(uint64_t *)&bitmask[i];
   }
#endif
};

void bitwise_or(Bitmap *bitmap_a, Bitmap *bitmap_b, Bitmap *out, int vector_size)
{
#ifdef __USING_AVX2__
   for (int i = 0; i < vector_size; i += 32)
   {
      _mm256_store_si256((__m256i *)&out[i], _mm256_or_si256(*(__m256i *)&bitmap_a[i], *(__m256i *)&bitmap_b[i]));
   }
#else
   for (int i = 0; i < vector_size; i += 8)
   {
      *(uint64_t *)&out[i] = *(uint64_t *)&bitmap_a[i] & *(uint64_t *)&bitmap_b[i];
   }
#endif
};

void bitwise_not(Bitmap *bitmap_a, Bitmap *out, int vector_size)
{
   for (int i = 0; i < vector_size; i += 8)
   {
      *(uint64_t *)&out[i] = ~*(uint64_t *)&bitmap_a[i];
   }
};

int find_first_index(Bitmap *bitmap, int vector_size)
{
#ifdef __USING_AVX2__
   const __m256i zero_vec = _mm256_set1_epi8(char(0));
   for (int i = 0; i < vector_size; i += 32)
   {
      __m256i cmp_mask = _mm256_cmpeq_epi8(*(__m256i *)&bitmap[i], zero_vec);
      uint32_t bitset = ~_mm256_movemask_epi8(cmp_mask);
      if (bitset)
      {
         const int r = __builtin_ctz(bitset);
         uint8_t mask = bitmap[r + i];
         int offset = vecDecodeTable[mask][0] - 1;
         return i * 8 + r * 8 + offset;
      }
   }
   return -1;
#else
   for (int i = 0; i < vector_size; i += 8)
   {
      uint64_t bitset = *(uint64_t *)(bitmap + i);
      if (bitset)
      {
         return i * 8 + __builtin_ctzll(bitset);
      }
   }
   return -1;
#endif
};


int maskzor_find_first_index(Bitmap * bitmask, Bitmap *bitmap_a, Bitmap * bitmap_b, Bitmap *simd_buffer, int vector_size){
   const __m256i zero_vec = _mm256_set1_epi8(char(0));
   for (int i = 0; i < vector_size; i += 32)
   {
      *(__m256i *)simd_buffer = _mm256_or_si256(*(__m256i *)(bitmap_a + i), *(__m256i *)(bitmap_b + i));
      *(__m256i *)simd_buffer = _mm256_andnot_si256(*(__m256i *)simd_buffer, zero_vec);
      uint32_t bitset = ~_mm256_movemask_epi8(*(__m256i *)simd_buffer);
      if (bitset)
      {
         const int r = __builtin_ctz(bitset);
         uint8_t mask = simd_buffer[r];
         int offset = vecDecodeTable[mask][0] - 1;
         return i * 8 + r * 8 + offset;
      }
   }
   return -1;
};

int maskzor_find_first_index(Bitmap * bitmask, Bitmap *bitmap_a, Bitmap * bitmap_b, int vector_size){

   for (int i = 0; i < vector_size; i += 8)
   {
      uint64_t bitset = *(uint64_t *)(bitmap_a + i) | *(uint64_t *)(bitmap_b + i) ;
      bitset &= ~*(uint64_t *)(bitmask + i);
      if (bitset)
      {
         return i * 8 + __builtin_ctzll(bitset);
      }
   }
   return -1;
}
int find_last_index(Bitmap *bitmap, int vector_size)
{
   for (int i = vector_size - 8; i >= 0; i-=8)
   {
      uint64_t bitset = *(uint64_t *)(bitmap + i);
      if (bitset)
      {
         return i * 8 + 63 - __builtin_clzll(bitset);
      }
   }
   return -1;
};

void fill_with_one(Bitmap *bitmap, int num_one)
{
   if (num_one < 8)
   {
      *bitmap = BITMASK >> (8 - num_one);
      return;
   }
   int num_bytes = num_one / 8, offset = num_one % 8;
   memset(bitmap, 0xff, num_bytes);
   if (offset)
   {
      bitmap[num_bytes] = BITMASK >> (8 - offset);
   }
};

void mark_as_one(Bitmap *bitmap, int index)
{
   bitmap[index / 8] |= 1 << (index & 7);
};

void mark_as_zero(Bitmap *bitmap, int index)
{
   bitmap[index / 8] &= ~(1 << (index & 7));
};