/* Copyright (c) 2013 Scott Lembcke and Howling Moon Software
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef CHIPMUNK_TRANSFORM_H
#define CHIPMUNK_TRANSFORM_H

#include "chipmunk_types.h"
#include "cpVect.h"
#include "cpBB.h"

#include "fix14.h"

/// Identity transform matrix.
static const cpTransform cpTransformIdentity = {(1 << 14), 0, 0, (1 << 14), 0, 0};

/// Construct a new transform matrix.
/// (a, b) is the x basis vector.
/// (c, d) is the y basis vector.
/// (tx, ty) is the translation.
static inline cpTransform
cpTransformNew(cpFloat a, cpFloat b, cpFloat c, cpFloat d, cpFloat tx, cpFloat ty)
{
	cpTransform t = {int_to_fix14(a), int_to_fix14(b), int_to_fix14(c), int_to_fix14(d), int_to_fix14(tx), int_to_fix14(ty)};
	return t;
}

/// Construct a new transform matrix in transposed order.
static inline cpTransform
cpTransformNewTranspose(cpFloat a, cpFloat c, cpFloat tx, cpFloat b, cpFloat d, cpFloat ty)
{
	cpTransform t = {int_to_fix14(a), int_to_fix14(b), int_to_fix14(c), int_to_fix14(d), int_to_fix14(tx), int_to_fix14(ty)};
	return t;
}

/// Get the inverse of a transform matrix.
static inline cpTransform
cpTransformInverse(cpTransform t)
{
    fix14_t inv_det = fix14_inverse(fix14_mul(t.a, t.d) - fix14_mul(t.c, t.b));
    return cpTransformNewTranspose(
        fix14_mul(t.d, inv_det), fix14_mul(-t.c, inv_det), fix14_mul((fix14_mul(t.c, t.ty) - fix14_mul(t.tx, t.d)), inv_det),
       fix14_mul(-t.b, inv_det),  fix14_mul(t.a, inv_det), fix14_mul((fix14_mul(t.tx, t.b) - fix14_mul(t.a, t.ty)), inv_det)
    );
}

/// Multiply two transformation matrices.
static inline cpTransform
cpTransformMult(cpTransform t1, cpTransform t2)
{
    return cpTransformNewTranspose(
        fix14_mul(t1.a, t2.a) + fix14_mul(t1.c, t2.b), fix14_mul(t1.a, t2.c) + fix14_mul(t1.c, t2.d), fix14_mul(t1.a, t2.tx) + fix14_mul(t1.c, t2.ty) + t1.tx,
        fix14_mul(t1.b, t2.a) + fix14_mul(t1.d, t2.b), fix14_mul(t1.b, t2.c) + fix14_mul(t1.d, t2.d), fix14_mul(t1.b, t2.tx) + fix14_mul(t1.d, t2.ty) + t1.ty
    );
}

/// Transform an absolute point. (i.e. a vertex)
static inline cpVect
cpTransformPoint(cpTransform t, cpVect p)
{
    return cpv(fix14_mul(t.a, p.x) + fix14_mul(t.c, p.y) + t.tx, fix14_mul(t.b, p.x) + fix14_mul(t.d, p.y) + t.ty);
}

/// Transform a vector (i.e. a normal)
static inline cpVect
cpTransformVect(cpTransform t, cpVect v)
{
    return cpv(fix14_mul(t.a, v.x) + fix14_mul(t.c, v.y), fix14_mul(t.b, v.x) + fix14_mul(t.d, v.y));
}

/// Transform a cpBB.
static inline cpBB
cpTransformbBB(cpTransform t, cpBB bb)
{
	cpVect center = cpBBCenter(bb);
	fix14_t hw = (bb.r - bb.l) >> 1;
	fix14_t hh = (bb.t - bb.b) >> 1;
	
	fix14_t a = fix14_mul(t.a, hw), b = fix14_mul(t.c, hh), d = fix14_mul(t.b, hw), e = fix14_mul(t.d, hh);
	fix14_t hw_max = cpfmax(cpfabs(a + b), cpfabs(a - b));
	fix14_t hh_max = cpfmax(cpfabs(d + e), cpfabs(d - e));
	return cpBBNewForExtents(cpTransformPoint(t, center), hw_max, hh_max);
}

/// Create a transation matrix.
static inline cpTransform
cpTransformTranslate(cpVect translate)
{
    return cpTransformNewTranspose(
        int_to_fix14(1), 0, translate.x,
        0, int_to_fix14(1), translate.y
    );
}

/// Create a scale matrix.
static inline cpTransform
cpTransformScale(cpFloat scaleX, cpFloat scaleY)
{
	return cpTransformNewTranspose(
		int_to_fix14(scaleX), 0, 0,
		   0, int_to_fix14(scaleY), 0
	);
}

/// Create a rotation matrix.
static inline cpTransform
cpTransformRotate(cpFloat radians)
{
	cpVect rot = cpvforangle(radians);
	return cpTransformNewTranspose(
		rot.x, fix14_mul(int_to_fix14(-1), rot.y), 0,
		rot.y,  rot.x, 0
	);
}

/// Create a rigid transformation matrix. (transation + rotation)
static inline cpTransform
cpTransformRigid(cpVect translate, cpFloat radians)
{
	cpVect rot = cpvforangle(radians);
	return cpTransformNewTranspose(
		rot.x, -rot.y, translate.x,
		rot.y,  rot.x, translate.y
	);
}

/// Fast inverse of a rigid transformation matrix.
static inline cpTransform
cpTransformRigidInverse(cpTransform t)
{
  return cpTransformNewTranspose(
     t.d, -t.c, (fix14_mul(t.c, t.ty) - fix14_mul(t.tx, t.d)),
    -t.b,  t.a, (fix14_mul(t.tx, t.b) - fix14_mul(t.a, t.ty))
  );
}

//MARK: Miscellaneous (but useful) transformation matrices.
// See source for documentation...

static inline cpTransform
cpTransformWrap(cpTransform outer, cpTransform inner)
{
    return cpTransformMult(cpTransformInverse(outer), cpTransformMult(inner, outer));
}

static inline cpTransform
cpTransformWrapInverse(cpTransform outer, cpTransform inner)
{
    return cpTransformMult(outer, cpTransformMult(inner, cpTransformInverse(outer)));
}

static inline cpTransform
cpTransformOrtho(cpBB bb)
{
    return cpTransformNewTranspose(
        fix14_div(int_to_fix14(2), (bb.r - bb.l)), 0, fix14_div(-(bb.r + bb.l), (bb.r - bb.l)),
        0, fix14_div(int_to_fix14(2), (bb.t - bb.b)), fix14_div(-(bb.t + bb.b), (bb.t - bb.b))
    );
}

static inline cpTransform
cpTransformBoneScale(cpVect v0, cpVect v1)
{
    cpVect d = cpvsub(v1, v0); 
    return cpTransformNewTranspose(
        d.x, fix14_mul(int_to_fix14(-1), d.y), v0.x,
        d.y,  d.x, v0.y
    );
}

static inline cpTransform
cpTransformAxialScale(cpVect axis, cpVect pivot, cpFloat scale)
{
    fix14_t A = fix14_mul(fix14_mul(axis.x, axis.y), scale - int_to_fix14(1));
    fix14_t B = fix14_mul(cpvdot(axis, pivot), int_to_fix14(1) - scale);
    
    return cpTransformNewTranspose(
        fix14_mul(scale, fix14_mul(axis.x, axis.x)) + fix14_mul(axis.y, axis.y), A, fix14_mul(axis.x, B),
        A, fix14_mul(axis.x, axis.x) + fix14_mul(scale, fix14_mul(axis.y, axis.y)), fix14_mul(axis.y, B)
    );
}

#endif
