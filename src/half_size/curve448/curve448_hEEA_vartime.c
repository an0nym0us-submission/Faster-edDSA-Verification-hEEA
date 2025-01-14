/*
 * Copyright (c) 2024 xxxxxxxxxx (xxxxx@xxxxx). All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following disclaimer
 *    in the documentation and/or other materials provided with the
 *    distribution.
 *  * Neither the name of Google Inc. nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Inspired by https://github.com/pornin/curve9767/blob/master/src/scalar_amd64.c
*/

#include "curve448_hEEA_vartime.h"

/*
 * We use the intrinsic functions:
 *   _lzcnt_u32(), _lzcnt_u64(), _addcarry_u64(), _subborrow_u64().
 */
#include <immintrin.h>
// #include <stdio.h>



/*
 * This macro defines a function with prototype:
 *  static inline void name(uint64_t *a, const uint64_t *b, unsigned s)
 * It adds lshift(b,s) to a (if intrinsic_op is _addcarry_u64), or
 * subtracts lshift(b,s) from a (if intrinsic_op is _subborrow_u64).
 * Both values consist of 'size' limbs. Truncation happens if the result
 * does not fit in the output.
 */
#define DEF_OP_LSHIFT(name, intrinsic_op, size) \
static void \
name(uint64_t *a, const uint64_t *b, unsigned s) \
{ \
	uint64_t b2[(size)]; \
	unsigned char cc; \
	unsigned long long w, w2, e; \
	size_t i; \
 \
 	if (s >= 64) { \
		unsigned k; \
 \
		k = s >> 6; \
		s &= 63; \
		if (k >= (size)) { \
			return; \
		} \
		memset(b2, 0, k * sizeof b2[0]); \
		memcpy(b2 + k, b, ((size) - k) * sizeof b2[0]); \
		b = b2; \
	} \
	if (s == 0) { \
		cc = 0; \
		for (i = 0; i < (size); i ++) { \
			cc = intrinsic_op(cc, a[i], b[i], &w); \
			a[i] = w; \
		} \
	} else { \
		cc = 0; \
		e = 0; \
		for (i = 0; i < (size); i ++) { \
			w = b[i]; \
			cc = intrinsic_op(cc, a[i], (w << s) | e, &w2); \
			e = w >> (64 - s); \
			a[i] = w2; \
		} \
	} \
}

DEF_OP_LSHIFT(add_lshift_7, _addcarry_u64, 7)
DEF_OP_LSHIFT(sub_lshift_7, _subborrow_u64, 7)
DEF_OP_LSHIFT(add_lshift_6, _addcarry_u64, 6)
DEF_OP_LSHIFT(sub_lshift_6, _subborrow_u64, 6)
DEF_OP_LSHIFT(add_lshift_5, _addcarry_u64, 5)
DEF_OP_LSHIFT(sub_lshift_5, _subborrow_u64, 5)
DEF_OP_LSHIFT(add_lshift_4, _addcarry_u64, 4)
DEF_OP_LSHIFT(sub_lshift_4, _subborrow_u64, 4)



/* return c0 and c1 s.t. c1 * v = c0 mod el */
void
curve448_hEEA_vartime(
	uint64_t *c0, uint64_t *c1, const uint64_t *v)
{

	unsigned long long r2_0, r2_1, r2_2, r2_3, r2_4, r2_5, r2_6;
	unsigned long long r1_0, r1_1, r1_2, r1_3, r1_4, r1_5, r1_6;
	unsigned long long t2_0, t2_1, t2_2, t2_3;
	unsigned long long t1_0, t1_1, t1_2, t1_3;

	unsigned char cc;


	/*
	el = 2**446 - 13818066809895115352007386748515426880336692474882178609894547503885
	   = 181709681073901722637330951972001133588410340171829515070372549795146003961539585716195755291692375963310293709091662304773755859649779
	   = 0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb49aed63690216cc2728dc58f552378c292ab5844f3
	L = 1 * el
	L = 0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb49aed63690216cc2728dc58f552378c292ab5844f3
	446 = L.bit_lenght() 
	*/
	/* L = 1* el */
	static const uint64_t L[] = {
		0x2378c292ab5844f3,
		0x216cc2728dc58f55,
		0xc44edb49aed63690,
		0xffffffff7cca23e9,
		0xffffffffffffffff,
		0xffffffffffffffff,
		0x3fffffffffffffff
	};


	/*
	 * Algorithm:
	 * r2 = L, r1 = v, t2 = 0, t1 = 1
	 * lenr1 = len(r2)
	 * while True
	 * 	lenr2 = lenr1
	 *	lenr1 = len(r1)
	 * 	if lenr1 <= target
	 * 		return c0 = r1, c1 = t1
	 * 	if (sign(r2) = sign(r1))
	 *		r = r2 - (r1 << s)
	 *		t = t2 - (t1 << s)
	 * 	else
	 * 		r = r2 + (r1 << s)
	 *		t = t2 + (t1 << s)
	 * 	r2 = r1
	 * 	r1 = r
	 * 	t2 = t1
	 * 	t1 = t
	 */
	r2_0 = L[0];
	r2_1 = L[1];
	r2_2 = L[2];
	r2_3 = L[3];
	r2_4 = L[4];
	r2_5 = L[5];
	r2_6 = L[6];

	r1_0 = v[0];
	r1_1 = v[1];
	r1_2 = v[2];
	r1_3 = v[3];
	r1_4 = v[4];
	r1_5 = v[5];
	r1_6 = v[6];

	t2_0 = 0;
	t2_1 = 0;
	t2_2 = 0;
	t2_3 = 0;

	t1_0 = 1;
	t1_1 = 0;
	t1_2 = 0;
	t1_3 = 0;

	
	unsigned bl_r2, bl_r1, bl_r;
	bl_r2 = 446;

#define BITLENGTH(size, bb)   do { \
		unsigned bitlength_acc = 7; \
		int bitlength_flag = 1; \
		unsigned long long bitlength_mask = -(bb ## 6 >> 63); \
		unsigned long long bitlength_top = 0; \
		unsigned long long bitlength_word; \
		bitlength_word = bb ## 6 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 5 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 4 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 3 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 2 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 1 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 0 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		(size) = (bitlength_acc << 6) - _lzcnt_u64(bitlength_top); \
	} while (0)

	BITLENGTH(bl_r1, r1_);

	for (;;) {
		unsigned s;		
		/*
		 * If r2 & r1 are small enough to fit into 6 limbs, jump to the shrink implementation.
		 */
		if ((bl_r2 < 384))
			break;


		/*
		 * Compute shift amount.
		 */
		s = bl_r2 - bl_r1;

		/*
		 * It is very rare that s > 31. We handle it with some
		 * generic code; branch prediction will soon learn that
		 * this path is normally not taken.
		 */
		if (s > 31) {
			uint64_t tr1[7], tr[7];
			uint64_t tt1[4], tt[4];
			
			tr[0] = r2_0;
			tr[1] = r2_1;
			tr[2] = r2_2;
			tr[3] = r2_3;
			tr[4] = r2_4;
			tr[5] = r2_5;
			tr[6] = r2_6;

			tr1[0] = r1_0;
			tr1[1] = r1_1;
			tr1[2] = r1_2;
			tr1[3] = r1_3;
			tr1[4] = r1_4;
			tr1[5] = r1_5;
			tr1[6] = r1_6;

			tt[0] = t2_0;
			tt[1] = t2_1;
			tt[2] = t2_2;
			tt[3] = t2_3;

			tt1[0] = t1_0;
			tt1[1] = t1_1;
			tt1[2] = t1_2;
			tt1[3] = t1_3;

			if ((r2_6 >> 63) == (r1_6 >> 63)){
				sub_lshift_7(tr, tr1, s);
				sub_lshift_4(tt, tt1, s);
			} else {
				add_lshift_7(tr, tr1, s);
				add_lshift_4(tt, tt1, s);
			}
			unsigned long long r_0, r_1, r_2, r_3, r_4, r_5, r_6;
			r_0 = tr[0];
			r_1 = tr[1];
			r_2 = tr[2];
			r_3 = tr[3];
			r_4 = tr[4];
			r_5 = tr[5];
			r_6 = tr[6];
			BITLENGTH(bl_r, r_);
			if (bl_r > bl_r1){
				r2_0 = r_0;
				r2_1 = r_1;
				r2_2 = r_2;
				r2_3 = r_3;
				r2_4 = r_4;
				r2_5 = r_5;
				r2_6 = r_6;
				t2_0 = tt[0];
				t2_1 = tt[1];
				t2_2 = tt[2];
				t2_3 = tt[3];
				bl_r2 = bl_r;

			}else{
				r2_0 = r1_0;
				r2_1 = r1_1;
				r2_2 = r1_2;
				r2_3 = r1_3;
				r2_4 = r1_4;
				r2_5 = r1_5;
				r2_6 = r1_6;
				r1_0 = tr[0];
				r1_1 = tr[1];
				r1_2 = tr[2];
				r1_3 = tr[3];
				r1_4 = tr[4];
				r1_5 = tr[5];
				r1_6 = tr[6];
				t2_0 = t1_0;
				t2_1 = t1_1;
				t2_2 = t1_2;
				t2_3 = t1_3;
				t1_0 = tt[0];
				t1_1 = tt[1];
				t1_2 = tt[2];
				t1_3 = tt[3];
				bl_r2 = bl_r1;
				bl_r1 = bl_r;
			}
			
			continue;
		}

	unsigned long long r_0, r_1, r_2, r_3, r_4, r_5, r_6, t_0, t_1, t_2, t_3;
		if ((r2_6 >> 63) == (r1_6 >> 63)){
			if (s == 0) {
				cc = _subborrow_u64(0, r2_0, r1_0, &r_0);
				cc = _subborrow_u64(cc, r2_1, r1_1, &r_1);
				cc = _subborrow_u64(cc, r2_2, r1_2, &r_2);
				cc = _subborrow_u64(cc, r2_3, r1_3, &r_3);
				cc = _subborrow_u64(cc, r2_4, r1_4, &r_4);
				cc = _subborrow_u64(cc, r2_5, r1_5, &r_5);
				(void)_subborrow_u64(cc, r2_6, r1_6, &r_6);

				cc = _subborrow_u64(0, t2_0, t1_0, &t_0);
				cc = _subborrow_u64(cc, t2_1, t1_1, &t_1);
				cc = _subborrow_u64(cc, t2_2, t1_2, &t_2);
				(void)_subborrow_u64(cc, t2_3, t1_3, &t_3);

			} else {
				unsigned rs;
				rs = 64 - s;

				cc = _subborrow_u64(0, r2_0,
					(r1_0 << s), &r_0);
				cc = _subborrow_u64(cc, r2_1,
					(r1_1 << s) | (r1_0 >> rs), &r_1);
				cc = _subborrow_u64(cc, r2_2,
					(r1_2 << s) | (r1_1 >> rs), &r_2);
				cc = _subborrow_u64(cc, r2_3,
					(r1_3 << s) | (r1_2 >> rs), &r_3);
				cc = _subborrow_u64(cc, r2_4,
					(r1_4 << s) | (r1_3 >> rs), &r_4);
				cc = _subborrow_u64(cc, r2_5,
					(r1_5 << s) | (r1_4 >> rs), &r_5);
				(void)_subborrow_u64(cc, r2_6,
					(r1_6 << s) | (r1_5 >> rs), &r_6);

				cc = _subborrow_u64(0, t2_0,
					(t1_0 << s), &t_0);
				cc = _subborrow_u64(cc, t2_1,
					(t1_1 << s) | (t1_0 >> rs), &t_1);
				cc = _subborrow_u64(cc, t2_2,
					(t1_2 << s) | (t1_1 >> rs), &t_2);
				(void)_subborrow_u64(cc, t2_3,
					(t1_3 << s) | (t1_2 >> rs), &t_3);
			}
		} else {
			if (s == 0) {
				cc = _addcarry_u64(0, r2_0, r1_0, &r_0);
				cc = _addcarry_u64(cc, r2_1, r1_1, &r_1);
				cc = _addcarry_u64(cc, r2_2, r1_2, &r_2);
				cc = _addcarry_u64(cc, r2_3, r1_3, &r_3);
				cc = _addcarry_u64(cc, r2_4, r1_4, &r_4);
				cc = _addcarry_u64(cc, r2_5, r1_5, &r_5);
				(void)_addcarry_u64(cc, r2_6, r1_6, &r_6);

				cc = _addcarry_u64(0, t2_0, t1_0, &t_0);
				cc = _addcarry_u64(cc, t2_1, t1_1, &t_1);
				cc = _addcarry_u64(cc, t2_2, t1_2, &t_2);
				(void)_addcarry_u64(cc, t2_3, t1_3, &t_3);

			} else {
				unsigned rs;
				rs = 64 - s;

				cc = _addcarry_u64(0, r2_0,
					(r1_0 << s), &r_0);
				cc = _addcarry_u64(cc, r2_1,
					(r1_1 << s) | (r1_0 >> rs), &r_1);
				cc = _addcarry_u64(cc, r2_2,
					(r1_2 << s) | (r1_1 >> rs), &r_2);
				cc = _addcarry_u64(cc, r2_3,
					(r1_3 << s) | (r1_2 >> rs), &r_3);
				cc = _addcarry_u64(cc, r2_4,
					(r1_4 << s) | (r1_3 >> rs), &r_4);
				cc = _addcarry_u64(cc, r2_5,
					(r1_5 << s) | (r1_4 >> rs), &r_5);
				(void)_addcarry_u64(cc, r2_6,
					(r1_6 << s) | (r1_5 >> rs), &r_6);

				cc = _addcarry_u64(0, t2_0,
					(t1_0 << s), &t_0);
				cc = _addcarry_u64(cc, t2_1,
					(t1_1 << s) | (t1_0 >> rs), &t_1);
				cc = _addcarry_u64(cc, t2_2,
					(t1_2 << s) | (t1_1 >> rs), &t_2);
				(void)_addcarry_u64(cc, t2_3,
					(t1_3 << s) | (t1_2 >> rs), &t_3);
			}
		}
		BITLENGTH(bl_r, r_);
		if(bl_r > bl_r1){
			r2_0 = r_0;
			r2_1 = r_1;
			r2_2 = r_2;
			r2_3 = r_3;
			r2_4 = r_4;
			r2_5 = r_5;
			r2_6 = r_6;
			t2_0 = t_0;
			t2_1 = t_1;
			t2_2 = t_2;
			t2_3 = t_3;
			bl_r2 = bl_r;
		}else{
			r2_0 = r1_0;
			r2_1 = r1_1;
			r2_2 = r1_2;
			r2_3 = r1_3;
			r2_4 = r1_4;
			r2_5 = r1_5;
			r2_6 = r1_6;
			r1_0 = r_0;
			r1_1 = r_1;
			r1_2 = r_2;
			r1_3 = r_3;
			r1_4 = r_4;
			r1_5 = r_5;
			r1_6 = r_6;
			t2_0 = t1_0;
			t2_1 = t1_1;
			t2_2 = t1_2;
			t2_3 = t1_3;
			t1_0 = t_0;
			t1_1 = t_1;
			t1_2 = t_2;
			t1_3 = t_3;
			bl_r2 = bl_r1;
			bl_r1 = bl_r;
		}
		
	}

	/*
	 * The part below is reached when r1 and r2 are small engouh to fit in 6 limbs.
	 */

	for (;;) {
		unsigned s;

#define BITLENGTH_SHRUNK6(size, bb)   do { \
		unsigned bitlength_acc = 6; \
		int bitlength_flag = 1; \
		unsigned long long bitlength_mask = -(bb ## 5 >> 63); \
		unsigned long long bitlength_top = 0; \
		unsigned long long bitlength_word; \
		bitlength_word = bb ## 5 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 4 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 3 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 2 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 1 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 0 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		(size) = (bitlength_acc << 6) - _lzcnt_u64(bitlength_top); \
	} while (0)

		/*
		 * If r2 & r1 are small enough to fit into 5 limbs, jump to the shrink implementation.
		 */
		if ((bl_r2 < 320))
			break;


		/*
		 * Compute shift amount.
		 */
		s = bl_r2 - bl_r1;

		/*
		 * It is very rare that s > 31. We handle it with some
		 * generic code; branch prediction will soon learn that
		 * this path is normally not taken.
		 */
		if (s > 31) {
			uint64_t tr1[6], tr[6];
			uint64_t tt1[4], tt[4];
			
			tr[0] = r2_0;
			tr[1] = r2_1;
			tr[2] = r2_2;
			tr[3] = r2_3;
			tr[4] = r2_4;
			tr[5] = r2_5;
			// tr[6] = r2_6;

			tr1[0] = r1_0;
			tr1[1] = r1_1;
			tr1[2] = r1_2;
			tr1[3] = r1_3;
			tr1[4] = r1_4;
			tr1[5] = r1_5;
			// tr1[6] = r1_6;

			tt[0] = t2_0;
			tt[1] = t2_1;
			tt[2] = t2_2;
			tt[3] = t2_3;

			tt1[0] = t1_0;
			tt1[1] = t1_1;
			tt1[2] = t1_2;
			tt1[3] = t1_3;

			if ((r2_5 >> 63) == (r1_5 >> 63)){
				sub_lshift_6(tr, tr1, s);
				sub_lshift_4(tt, tt1, s);
			} else {
				add_lshift_6(tr, tr1, s);
				add_lshift_4(tt, tt1, s);
			}

			unsigned long long r_0, r_1, r_2, r_3, r_4, r_5;//, r_6;
			r_0 = tr[0];
			r_1 = tr[1];
			r_2 = tr[2];
			r_3 = tr[3];
			r_4 = tr[4];
			r_5 = tr[5];
			// r_6 = tr[6];
			BITLENGTH_SHRUNK6(bl_r, r_);
			if (bl_r > bl_r1){
				r2_0 = r_0;
				r2_1 = r_1;
				r2_2 = r_2;
				r2_3 = r_3;
				r2_4 = r_4;
				r2_5 = r_5;
				// r2_6 = r_6;
				t2_0 = tt[0];
				t2_1 = tt[1];
				t2_2 = tt[2];
				t2_3 = tt[3];
				bl_r2 = bl_r;

			}else{
				r2_0 = r1_0;
				r2_1 = r1_1;
				r2_2 = r1_2;
				r2_3 = r1_3;
				r2_4 = r1_4;
				r2_5 = r1_5;
				// r2_6 = r1_6;
				r1_0 = tr[0];
				r1_1 = tr[1];
				r1_2 = tr[2];
				r1_3 = tr[3];
				r1_4 = tr[4];
				r1_5 = tr[5];
				// r1_6 = tr[6];
				t2_0 = t1_0;
				t2_1 = t1_1;
				t2_2 = t1_2;
				t2_3 = t1_3;
				t1_0 = tt[0];
				t1_1 = tt[1];
				t1_2 = tt[2];
				t1_3 = tt[3];
				bl_r2 = bl_r1;
				bl_r1 = bl_r;
			}
			continue;
		}

	unsigned long long r_0, r_1, r_2, r_3, r_4, r_5, t_0, t_1, t_2, t_3;
		if ((r2_5 >> 63) == (r1_5 >> 63)){
			if (s == 0) {
				cc = _subborrow_u64(0, r2_0, r1_0, &r_0);
				cc = _subborrow_u64(cc, r2_1, r1_1, &r_1);
				cc = _subborrow_u64(cc, r2_2, r1_2, &r_2);
				cc = _subborrow_u64(cc, r2_3, r1_3, &r_3);
				cc = _subborrow_u64(cc, r2_4, r1_4, &r_4);
				(void)_subborrow_u64(cc, r2_5, r1_5, &r_5);

				cc = _subborrow_u64(0, t2_0, t1_0, &t_0);
				cc = _subborrow_u64(cc, t2_1, t1_1, &t_1);
				cc = _subborrow_u64(cc, t2_2, t1_2, &t_2);
				(void)_subborrow_u64(cc, t2_3, t1_3, &t_3);

			} else {
				unsigned rs;
				rs = 64 - s;

				cc = _subborrow_u64(0, r2_0,
					(r1_0 << s), &r_0);
				cc = _subborrow_u64(cc, r2_1,
					(r1_1 << s) | (r1_0 >> rs), &r_1);
				cc = _subborrow_u64(cc, r2_2,
					(r1_2 << s) | (r1_1 >> rs), &r_2);
				cc = _subborrow_u64(cc, r2_3,
					(r1_3 << s) | (r1_2 >> rs), &r_3);
				cc = _subborrow_u64(cc, r2_4,
					(r1_4 << s) | (r1_3 >> rs), &r_4);
				(void)_subborrow_u64(cc, r2_5,
					(r1_5 << s) | (r1_4 >> rs), &r_5);

				cc = _subborrow_u64(0, t2_0,
					(t1_0 << s), &t_0);
				cc = _subborrow_u64(cc, t2_1,
					(t1_1 << s) | (t1_0 >> rs), &t_1);
				cc = _subborrow_u64(cc, t2_2,
					(t1_2 << s) | (t1_1 >> rs), &t_2);
				(void)_subborrow_u64(cc, t2_3,
					(t1_3 << s) | (t1_2 >> rs), &t_3);
			}
		} else {
			if (s == 0) {
				cc = _addcarry_u64(0, r2_0, r1_0, &r_0);
				cc = _addcarry_u64(cc, r2_1, r1_1, &r_1);
				cc = _addcarry_u64(cc, r2_2, r1_2, &r_2);
				cc = _addcarry_u64(cc, r2_3, r1_3, &r_3);
				cc = _addcarry_u64(cc, r2_4, r1_4, &r_4);
				(void)_addcarry_u64(cc, r2_5, r1_5, &r_5);

				cc = _addcarry_u64(0, t2_0, t1_0, &t_0);
				cc = _addcarry_u64(cc, t2_1, t1_1, &t_1);
				cc = _addcarry_u64(cc, t2_2, t1_2, &t_2);
				(void)_addcarry_u64(cc, t2_3, t1_3, &t_3);

			} else {
				unsigned rs;
				rs = 64 - s;

				cc = _addcarry_u64(0, r2_0,
					(r1_0 << s), &r_0);
				cc = _addcarry_u64(cc, r2_1,
					(r1_1 << s) | (r1_0 >> rs), &r_1);
				cc = _addcarry_u64(cc, r2_2,
					(r1_2 << s) | (r1_1 >> rs), &r_2);
				cc = _addcarry_u64(cc, r2_3,
					(r1_3 << s) | (r1_2 >> rs), &r_3);
				cc = _addcarry_u64(cc, r2_4,
					(r1_4 << s) | (r1_3 >> rs), &r_4);
				(void)_addcarry_u64(cc, r2_5,
					(r1_5 << s) | (r1_4 >> rs), &r_5);

				cc = _addcarry_u64(0, t2_0,
					(t1_0 << s), &t_0);
				cc = _addcarry_u64(cc, t2_1,
					(t1_1 << s) | (t1_0 >> rs), &t_1);
				cc = _addcarry_u64(cc, t2_2,
					(t1_2 << s) | (t1_1 >> rs), &t_2);
				(void)_addcarry_u64(cc, t2_3,
					(t1_3 << s) | (t1_2 >> rs), &t_3);
			}
		}

		BITLENGTH_SHRUNK6(bl_r, r_);
		if(bl_r > bl_r1){
			r2_0 = r_0;
			r2_1 = r_1;
			r2_2 = r_2;
			r2_3 = r_3;
			r2_4 = r_4;
			r2_5 = r_5;
			// r2_6 = r_6;
			t2_0 = t_0;
			t2_1 = t_1;
			t2_2 = t_2;
			t2_3 = t_3;
			bl_r2 = bl_r;
		}else{
			r2_0 = r1_0;
			r2_1 = r1_1;
			r2_2 = r1_2;
			r2_3 = r1_3;
			r2_4 = r1_4;
			r2_5 = r1_5;
			// r2_6 = r1_6;
			r1_0 = r_0;
			r1_1 = r_1;
			r1_2 = r_2;
			r1_3 = r_3;
			r1_4 = r_4;
			r1_5 = r_5;
			// r1_6 = r_6;
			t2_0 = t1_0;
			t2_1 = t1_1;
			t2_2 = t1_2;
			t2_3 = t1_3;
			t1_0 = t_0;
			t1_1 = t_1;
			t1_2 = t_2;
			t1_3 = t_3;
			bl_r2 = bl_r1;
			bl_r1 = bl_r;
		}
	}

	/*
	 * The part below is reached when r1 and r2 are small engouh to fit in 5 limbs.
	 */

	for (;;) {
		unsigned s;

#define BITLENGTH_SHRUNK5(size, bb)   do { \
		unsigned bitlength_acc = 5; \
		int bitlength_flag = 1; \
		unsigned long long bitlength_mask = -(bb ## 4 >> 63); \
		unsigned long long bitlength_top = 0; \
		unsigned long long bitlength_word; \
		bitlength_word = bb ## 4 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 3 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 2 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 1 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 0 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		(size) = (bitlength_acc << 6) - _lzcnt_u64(bitlength_top); \
	} while (0)
		
		/*
		 * If r2 & r1 are small enough to fit into 4 limbs, jump to the shrink implementation.
		 */

		if ((bl_r2 < 256))
			break;


		/*
		 * Compute shift amount.
		 */
		s = bl_r2 - bl_r1;

		/*
		 * It is very rare that s > 31. We handle it with some
		 * generic code; branch prediction will soon learn that
		 * this path is normally not taken.
		 */
		if (s > 31) {
			uint64_t tr1[5], tr[5];
			uint64_t tt1[4], tt[4];
			
			tr[0] = r2_0;
			tr[1] = r2_1;
			tr[2] = r2_2;
			tr[3] = r2_3;
			tr[4] = r2_4;
			// tr[5] = r2_5;
			// tr[6] = r2_6;

			tr1[0] = r1_0;
			tr1[1] = r1_1;
			tr1[2] = r1_2;
			tr1[3] = r1_3;
			tr1[4] = r1_4;
			// tr1[5] = r1_5;
			// tr1[6] = r1_6;

			tt[0] = t2_0;
			tt[1] = t2_1;
			tt[2] = t2_2;
			tt[3] = t2_3;

			tt1[0] = t1_0;
			tt1[1] = t1_1;
			tt1[2] = t1_2;
			tt1[3] = t1_3;

			if ((r2_4 >> 63) == (r1_4 >> 63)){
				sub_lshift_5(tr, tr1, s);
				sub_lshift_4(tt, tt1, s);
			} else {
				add_lshift_5(tr, tr1, s);
				add_lshift_4(tt, tt1, s);
			}

			unsigned long long r_0, r_1, r_2, r_3, r_4;//, r_5;//, r_6;
			r_0 = tr[0];
			r_1 = tr[1];
			r_2 = tr[2];
			r_3 = tr[3];
			r_4 = tr[4];
			// r_5 = tr[5];
			// r_6 = tr[6];
			BITLENGTH_SHRUNK5(bl_r, r_);
			if (bl_r > bl_r1){
				r2_0 = r_0;
				r2_1 = r_1;
				r2_2 = r_2;
				r2_3 = r_3;
				r2_4 = r_4;
				// r2_5 = r_5;
				// r2_6 = r_6;
				t2_0 = tt[0];
				t2_1 = tt[1];
				t2_2 = tt[2];
				t2_3 = tt[3];
				bl_r2 = bl_r;

			}else{
				r2_0 = r1_0;
				r2_1 = r1_1;
				r2_2 = r1_2;
				r2_3 = r1_3;
				r2_4 = r1_4;
				// r2_5 = r1_5;
				// r2_6 = r1_6;
				r1_0 = tr[0];
				r1_1 = tr[1];
				r1_2 = tr[2];
				r1_3 = tr[3];
				r1_4 = tr[4];
				// r1_5 = tr[5];
				// r1_6 = tr[6];
				t2_0 = t1_0;
				t2_1 = t1_1;
				t2_2 = t1_2;
				t2_3 = t1_3;
				t1_0 = tt[0];
				t1_1 = tt[1];
				t1_2 = tt[2];
				t1_3 = tt[3];
				bl_r2 = bl_r1;
				bl_r1 = bl_r;
			}
			continue;
		}

	unsigned long long r_0, r_1, r_2, r_3, r_4, t_0, t_1, t_2, t_3;
		if ((r2_4 >> 63) == (r1_4 >> 63)){
			if (s == 0) {
				cc = _subborrow_u64(0, r2_0, r1_0, &r_0);
				cc = _subborrow_u64(cc, r2_1, r1_1, &r_1);
				cc = _subborrow_u64(cc, r2_2, r1_2, &r_2);
				cc = _subborrow_u64(cc, r2_3, r1_3, &r_3);
				(void)_subborrow_u64(cc, r2_4, r1_4, &r_4);

				cc = _subborrow_u64(0, t2_0, t1_0, &t_0);
				cc = _subborrow_u64(cc, t2_1, t1_1, &t_1);
				cc = _subborrow_u64(cc, t2_2, t1_2, &t_2);
				(void)_subborrow_u64(cc, t2_3, t1_3, &t_3);

			} else {
				unsigned rs;
				rs = 64 - s;

				cc = _subborrow_u64(0, r2_0,
					(r1_0 << s), &r_0);
				cc = _subborrow_u64(cc, r2_1,
					(r1_1 << s) | (r1_0 >> rs), &r_1);
				cc = _subborrow_u64(cc, r2_2,
					(r1_2 << s) | (r1_1 >> rs), &r_2);
				cc = _subborrow_u64(cc, r2_3,
					(r1_3 << s) | (r1_2 >> rs), &r_3);
				(void)_subborrow_u64(cc, r2_4,
					(r1_4 << s) | (r1_3 >> rs), &r_4);

				cc = _subborrow_u64(0, t2_0,
					(t1_0 << s), &t_0);
				cc = _subborrow_u64(cc, t2_1,
					(t1_1 << s) | (t1_0 >> rs), &t_1);
				cc = _subborrow_u64(cc, t2_2,
					(t1_2 << s) | (t1_1 >> rs), &t_2);
				(void)_subborrow_u64(cc, t2_3,
					(t1_3 << s) | (t1_2 >> rs), &t_3);
			}
		} else {
			if (s == 0) {
				cc = _addcarry_u64(0, r2_0, r1_0, &r_0);
				cc = _addcarry_u64(cc, r2_1, r1_1, &r_1);
				cc = _addcarry_u64(cc, r2_2, r1_2, &r_2);
				cc = _addcarry_u64(cc, r2_3, r1_3, &r_3);
				(void)_addcarry_u64(cc, r2_4, r1_4, &r_4);

				cc = _addcarry_u64(0, t2_0, t1_0, &t_0);
				cc = _addcarry_u64(cc, t2_1, t1_1, &t_1);
				cc = _addcarry_u64(cc, t2_2, t1_2, &t_2);
				(void)_addcarry_u64(cc, t2_3, t1_3, &t_3);

			} else {
				unsigned rs;
				rs = 64 - s;

				cc = _addcarry_u64(0, r2_0,
					(r1_0 << s), &r_0);
				cc = _addcarry_u64(cc, r2_1,
					(r1_1 << s) | (r1_0 >> rs), &r_1);
				cc = _addcarry_u64(cc, r2_2,
					(r1_2 << s) | (r1_1 >> rs), &r_2);
				cc = _addcarry_u64(cc, r2_3,
					(r1_3 << s) | (r1_2 >> rs), &r_3);
				(void)_addcarry_u64(cc, r2_4,
					(r1_4 << s) | (r1_3 >> rs), &r_4);

				cc = _addcarry_u64(0, t2_0,
					(t1_0 << s), &t_0);
				cc = _addcarry_u64(cc, t2_1,
					(t1_1 << s) | (t1_0 >> rs), &t_1);
				cc = _addcarry_u64(cc, t2_2,
					(t1_2 << s) | (t1_1 >> rs), &t_2);
				(void)_addcarry_u64(cc, t2_3,
					(t1_3 << s) | (t1_2 >> rs), &t_3);
			}
		}

		BITLENGTH_SHRUNK5(bl_r, r_);
		if(bl_r > bl_r1){
			r2_0 = r_0;
			r2_1 = r_1;
			r2_2 = r_2;
			r2_3 = r_3;
			r2_4 = r_4;
			// r2_5 = r_5;
			// r2_6 = r_6;
			t2_0 = t_0;
			t2_1 = t_1;
			t2_2 = t_2;
			t2_3 = t_3;
			bl_r2 = bl_r;
		}else{
			r2_0 = r1_0;
			r2_1 = r1_1;
			r2_2 = r1_2;
			r2_3 = r1_3;
			r2_4 = r1_4;
			// r2_5 = r1_5;
			// r2_6 = r1_6;
			r1_0 = r_0;
			r1_1 = r_1;
			r1_2 = r_2;
			r1_3 = r_3;
			r1_4 = r_4;
			// r1_5 = r_5;
			// r1_6 = r_6;
			t2_0 = t1_0;
			t2_1 = t1_1;
			t2_2 = t1_2;
			t2_3 = t1_3;
			t1_0 = t_0;
			t1_1 = t_1;
			t1_2 = t_2;
			t1_3 = t_3;
			bl_r2 = bl_r1;
			bl_r1 = bl_r;
		}
	}


	/*
	 * The part below is reached when r1 and r2 are small engouh to fit in 4 limbs.
	 */

	for (;;) {
		unsigned s;

#define BITLENGTH_SHRUNK4(size, bb)   do { \
		unsigned bitlength_acc = 4; \
		int bitlength_flag = 1; \
		unsigned long long bitlength_mask = -(bb ## 3 >> 63); \
		unsigned long long bitlength_top = 0; \
		unsigned long long bitlength_word; \
		bitlength_word = bb ## 3 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 2 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 1 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 0 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		(size) = (bitlength_acc << 6) - _lzcnt_u64(bitlength_top); \
	} while (0)

		/*
		 * If r1 is small enough, return.
		 */

		if (bl_r1 <= 223) {
			c0[0] = r1_0;
			c0[1] = r1_1;
			c0[2] = r1_2;
			c0[3] = r1_3;
			c1[0] = t1_0;
			c1[1] = t1_1;
			c1[2] = t1_2;
			c1[3] = t1_3;
			return;
		} 
	

		/*
		 * Compute shift amount.
		 */
		s = bl_r2 - bl_r1;

		/*
		 * It is very rare that s > 31. We handle it with some
		 * generic code; branch prediction will soon learn that
		 * this path is normally not taken.
		 */
		if (s > 31) {
			uint64_t tr1[4], tr[4];
			uint64_t tt1[4], tt[4];
			
			tr[0] = r2_0;
			tr[1] = r2_1;
			tr[2] = r2_2;
			tr[3] = r2_3;
			// tr[4] = r2_4;
			// tr[5] = r2_5;
			// tr[6] = r2_6;

			tr1[0] = r1_0;
			tr1[1] = r1_1;
			tr1[2] = r1_2;
			tr1[3] = r1_3;
			// tr1[4] = r1_4;
			// tr1[5] = r1_5;
			// tr1[6] = r1_6;

			tt[0] = t2_0;
			tt[1] = t2_1;
			tt[2] = t2_2;
			tt[3] = t2_3;

			tt1[0] = t1_0;
			tt1[1] = t1_1;
			tt1[2] = t1_2;
			tt1[3] = t1_3;

			if ((r2_3 >> 63) == (r1_3 >> 63)){
				sub_lshift_4(tr, tr1, s);
				sub_lshift_4(tt, tt1, s);
			} else {
				add_lshift_4(tr, tr1, s);
				add_lshift_4(tt, tt1, s);
			}

			unsigned long long r_0, r_1, r_2, r_3;//, r_4;//, r_5;//, r_6;
			r_0 = tr[0];
			r_1 = tr[1];
			r_2 = tr[2];
			r_3 = tr[3];
			// r_4 = tr[4];
			// r_5 = tr[5];
			// r_6 = tr[6];
			BITLENGTH_SHRUNK4(bl_r, r_);
			if (bl_r > bl_r1){
				r2_0 = r_0;
				r2_1 = r_1;
				r2_2 = r_2;
				r2_3 = r_3;
				// r2_4 = r_4;
				// r2_5 = r_5;
				// r2_6 = r_6;
				t2_0 = tt[0];
				t2_1 = tt[1];
				t2_2 = tt[2];
				t2_3 = tt[3];
				bl_r2 = bl_r;

			}else{
				r2_0 = r1_0;
				r2_1 = r1_1;
				r2_2 = r1_2;
				r2_3 = r1_3;
				// r2_4 = r1_4;
				// r2_5 = r1_5;
				// r2_6 = r1_6;
				r1_0 = tr[0];
				r1_1 = tr[1];
				r1_2 = tr[2];
				r1_3 = tr[3];
				// r1_4 = tr[4];
				// r1_5 = tr[5];
				// r1_6 = tr[6];
				t2_0 = t1_0;
				t2_1 = t1_1;
				t2_2 = t1_2;
				t2_3 = t1_3;
				t1_0 = tt[0];
				t1_1 = tt[1];
				t1_2 = tt[2];
				t1_3 = tt[3];
				bl_r2 = bl_r1;
				bl_r1 = bl_r;
			}
			continue;
		}

	unsigned long long r_0, r_1, r_2, r_3, t_0, t_1, t_2, t_3;
		if ((r2_3 >> 63) == (r1_3 >> 63)){
			if (s == 0) {
				cc = _subborrow_u64(0, r2_0, r1_0, &r_0);
				cc = _subborrow_u64(cc, r2_1, r1_1, &r_1);
				cc = _subborrow_u64(cc, r2_2, r1_2, &r_2);
				(void)_subborrow_u64(cc, r2_3, r1_3, &r_3);

				cc = _subborrow_u64(0, t2_0, t1_0, &t_0);
				cc = _subborrow_u64(cc, t2_1, t1_1, &t_1);
				cc = _subborrow_u64(cc, t2_2, t1_2, &t_2);
				(void)_subborrow_u64(cc, t2_3, t1_3, &t_3);

			} else {
				unsigned rs;
				rs = 64 - s;

				cc = _subborrow_u64(0, r2_0,
					(r1_0 << s), &r_0);
				cc = _subborrow_u64(cc, r2_1,
					(r1_1 << s) | (r1_0 >> rs), &r_1);
				cc = _subborrow_u64(cc, r2_2,
					(r1_2 << s) | (r1_1 >> rs), &r_2);
				(void)_subborrow_u64(cc, r2_3,
					(r1_3 << s) | (r1_2 >> rs), &r_3);

				cc = _subborrow_u64(0, t2_0,
					(t1_0 << s), &t_0);
				cc = _subborrow_u64(cc, t2_1,
					(t1_1 << s) | (t1_0 >> rs), &t_1);
				cc = _subborrow_u64(cc, t2_2,
					(t1_2 << s) | (t1_1 >> rs), &t_2);
				(void)_subborrow_u64(cc, t2_3,
					(t1_3 << s) | (t1_2 >> rs), &t_3);
			}
		} else {
			if (s == 0) {
				cc = _addcarry_u64(0, r2_0, r1_0, &r_0);
				cc = _addcarry_u64(cc, r2_1, r1_1, &r_1);
				cc = _addcarry_u64(cc, r2_2, r1_2, &r_2);
				(void)_addcarry_u64(cc, r2_3, r1_3, &r_3);

				cc = _addcarry_u64(0, t2_0, t1_0, &t_0);
				cc = _addcarry_u64(cc, t2_1, t1_1, &t_1);
				cc = _addcarry_u64(cc, t2_2, t1_2, &t_2);
				(void)_addcarry_u64(cc, t2_3, t1_3, &t_3);

			} else {
				unsigned rs;
				rs = 64 - s;

				cc = _addcarry_u64(0, r2_0,
					(r1_0 << s), &r_0);
				cc = _addcarry_u64(cc, r2_1,
					(r1_1 << s) | (r1_0 >> rs), &r_1);
				cc = _addcarry_u64(cc, r2_2,
					(r1_2 << s) | (r1_1 >> rs), &r_2);
				(void)_addcarry_u64(cc, r2_3,
					(r1_3 << s) | (r1_2 >> rs), &r_3);

				cc = _addcarry_u64(0, t2_0,
					(t1_0 << s), &t_0);
				cc = _addcarry_u64(cc, t2_1,
					(t1_1 << s) | (t1_0 >> rs), &t_1);
				cc = _addcarry_u64(cc, t2_2,
					(t1_2 << s) | (t1_1 >> rs), &t_2);
				(void)_addcarry_u64(cc, t2_3,
					(t1_3 << s) | (t1_2 >> rs), &t_3);
			}
		}

		BITLENGTH_SHRUNK4(bl_r, r_);
		if(bl_r > bl_r1){
			r2_0 = r_0;
			r2_1 = r_1;
			r2_2 = r_2;
			r2_3 = r_3;
			// r2_4 = r_4;
			// r2_5 = r_5;
			// r2_6 = r_6;
			t2_0 = t_0;
			t2_1 = t_1;
			t2_2 = t_2;
			t2_3 = t_3;
			bl_r2 = bl_r;
		}else{
			r2_0 = r1_0;
			r2_1 = r1_1;
			r2_2 = r1_2;
			r2_3 = r1_3;
			// r2_4 = r1_4;
			// r2_5 = r1_5;
			// r2_6 = r1_6;
			r1_0 = r_0;
			r1_1 = r_1;
			r1_2 = r_2;
			r1_3 = r_3;
			// r1_4 = r_4;
			// r1_5 = r_5;
			// r1_6 = r_6;
			t2_0 = t1_0;
			t2_1 = t1_1;
			t2_2 = t1_2;
			t2_3 = t1_3;
			t1_0 = t_0;
			t1_1 = t_1;
			t1_2 = t_2;
			t1_3 = t_3;
			bl_r2 = bl_r1;
			bl_r1 = bl_r;
		}
	}


#undef BITLENGTH
#undef BITLENGTH_SHRUNK6
#undef BITLENGTH_SHRUNK5
#undef BITLENGTH_SHRUNK4
}