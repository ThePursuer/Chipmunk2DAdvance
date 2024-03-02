#ifndef FIXED_PRECISION14_H
#define FIXED_PRECISION14_H

#include <stdint.h>
#include <sys/cdefs.h>

// #include <gba_base.h>
// #include <gba_systemcalls.h>

#include "division_table_14_32.h"

typedef int32_t fix14_t;  // signed 18.14 fixed point
extern const uint32_t gSinCosTable[4096];
__always_inline fix14_t int_to_fix14(int a){
    return a << 14;
};
inline int fix14_to_int(fix14_t a){
    return a >> 14;
};
inline int fix14_to_int_rounded(fix14_t a){
    return (a + 8192) >> 14;
};
inline int32_t fix14_mul(int32_t a, int32_t b){
    return ((a >> 7) * (b >> 7));
};
inline fix14_t fix14_div(fix14_t a, fix14_t b){
    if(b < 0)
        return -(((int64_t)(a >> 7) * division_table_14_32[DIV_TABLE_14_32_INDEX(b)] >> 4) >> (3 + DIVISION_TABLE_14_32_PRECISION));
    return ((int64_t)(a >> 7) * division_table_14_32[DIV_TABLE_14_32_INDEX(b)] >> 4) >> (3 + DIVISION_TABLE_14_32_PRECISION);
};
inline fix14_t fix14_inverse(fix14_t a){
    division_table_14_32[DIV_TABLE_14_32_INDEX(a)] >> (DIVISION_TABLE_14_32_PRECISION);
}
inline fix14_t fix14_add(fix14_t a, fix14_t b){
    return a + b;
};
inline fix14_t fix14_sub(fix14_t a, fix14_t b){
    return a - b;
};
inline fix14_t fix14_div_div(fix14_t a, fix14_t b){
    return (int64_t)(a << 7) / (b >> 7);
}
inline fix14_t fix14_sqrt(fix14_t a){
    return (int16_t)(Sqrt((uint32_t)(a))) << 7;
};
inline fix14_t fix14_arctan(fix14_t a){
    a = (fix14_t)(ArcTan((int16_t)(fix14_to_int_rounded(a))));
    return int_to_fix14(a);
};
#define sincos(x,s,c) {\
    uint32 sc = gSinCosTable[uint32(x << 16) >> 20];\
    s = int32(sc) >> 16;\
    c = int32(sc) << 16 >> 16;\
}

#endif // FIXED_PRECISION14_H