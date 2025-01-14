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


int
ED25519_FN(ed25519_sign_open_hEEA) (const unsigned char *m, size_t mlen, const ed25519_public_key pk, const ed25519_signature RS) {
	ge25519 ALIGN(16) R, A, sumBRA;
	hash_512bits hash;
	bignum256modm hram, S1, S2={0}, r, t;
    int r_isneg;
	int t_isneg;

	if ((RS[63] & 224))
		return -1;

	/* hram <-- H(R,A,m) */
	ed25519_hram(hash, RS, pk, m, mlen);
	expand256_modm(hram, hash, 64);

	/* compute r and t s.t. rh = t mod el */
    curve25519_half_size_scalar_vartime_hEEA(r, t, hram, &r_isneg, &t_isneg);
    
    /* unpacking (-R) */
    if (!ge25519_unpack_negative_vartime(&R, RS))
		return -1;

	/*
		[S]B + (-R) + h(-A) =? 0 ==> [rS]B + [r](-R) + [t](-A) =? 0
		------------------------------------------------------------
		r_isneg | t_isneg | [rs]B + [r](-R) + [t](-A) =? 0
		------------------------------------------------------------
		   0    |    0    | [|r| * s]B + [|r|](-R) + [|t|](-A)
		------------------------------------------------------------
		   0    |    1    | [|r| * s]B + [|r|](-R) - [|t|](-A) ==>
		        |         | [|r| * s]B + [|r|](-R) + [|t|](A)
		------------------------------------------------------------
		   1    |    0    | [-|r| * s]B + [-|r|](-R) + [|t|](-A) ==> 
		        |         | [|r| * s] B + [|r|](-R) + [|t|](A)
		------------------------------------------------------------
		   1    |    1    | [-|r| * s]B + [-|r|](-R) + [-|t|](-A) ==> 
		        |         | [|r| * s] B + [|r|](-R) + [|t|](-A)
		------------------------------------------------------------
	*/
	/* unpacking (-A) and adjust the sign based on the previous table */
    if (r_isneg == t_isneg){
		if (!ge25519_unpack_negative_vartime(&A, pk))
			return -1;
	}else{
		if (!ge25519_unpack_positive_vartime(&A, pk))
			return -1;
	}
    
    /* S */
	expand256_modm(S1, RS + 32, 32);

    /* S <-- rS */
    mul256_modm(S1, S1, r);
    
    /* split S to S1 and S2, s.t. S = (S2<<126 | S1), bl(S1) = 126, bl(S2) = 127 */
    S2[0] = (S1[2] >> 14) | ((S1[3] & 0x3FFF) << 42);
    S2[1] = (S1[3] >> 14) | ((S1[4] & 0x3FFF) << 42);
    S2[2] = S1[4] >> 14;
    S1[2] &= 0x3FFF;
    S1[3] = 0;
    S1[4] = 0;

	/* 
	 * [S]B + (-R) + h(-A) =? 0 ==>
	 * [rS]B + [r](-R) + [t](-A) =? 0 ==>
	 * [S1]B + [S2]([2^126]B) + [r](-R) + [t](-A) =? 0 
	*/
	ge25519_quadruple_scalarmult_vartime(&sumBRA, &R, &A, r, t, S1, S2);

    return ge25519_is_neutral_vartime(&sumBRA) ? 0: -1;
}


int
ED25519_FN(ed25519_sign_open_hEEA_samePre) (const unsigned char *m, size_t mlen, const ed25519_public_key pk, const ed25519_signature RS) {
	ge25519 ALIGN(16) R, A, sumBRA;
	hash_512bits hash;
	bignum256modm hram, S1, S2={0}, r, t;
    int r_isneg;
	int t_isneg;

	if ((RS[63] & 224))
		return -1;

	/* hram <-- H(R,A,m) */
	ed25519_hram(hash, RS, pk, m, mlen);
	expand256_modm(hram, hash, 64);

	/* compute r and t s.t. rh = t mod el */
    curve25519_half_size_scalar_vartime_hEEA(r, t, hram, &r_isneg, &t_isneg);
    
    /* unpacking (-R) */
    if (!ge25519_unpack_negative_vartime(&R, RS))
		return -1;

	/*
		[S]B + (-R) + h(-A) =? 0 ==> [rS]B + [r](-R) + [t](-A) =? 0
		------------------------------------------------------------
		r_isneg | t_isneg | [rs]B + [r](-R) + [t](-A) =? 0
		------------------------------------------------------------
		   0    |    0    | [|r| * s]B + [|r|](-R) + [|t|](-A)
		------------------------------------------------------------
		   0    |    1    | [|r| * s]B + [|r|](-R) - [|t|](-A) ==>
		        |         | [|r| * s]B + [|r|](-R) + [|t|](A)
		------------------------------------------------------------
		   1    |    0    | [-|r| * s]B + [-|r|](-R) + [|t|](-A) ==> 
		        |         | [|r| * s] B + [|r|](-R) + [|t|](A)
		------------------------------------------------------------
		   1    |    1    | [-|r| * s]B + [-|r|](-R) + [-|t|](-A) ==> 
		        |         | [|r| * s] B + [|r|](-R) + [|t|](-A)
		------------------------------------------------------------
	*/
	/* unpacking (-A) and adjust the sign based on the previous table */
    if (r_isneg == t_isneg){
		if (!ge25519_unpack_negative_vartime(&A, pk))
			return -1;
	}else{
		if (!ge25519_unpack_positive_vartime(&A, pk))
			return -1;
	}
    
    /* S */
	expand256_modm(S1, RS + 32, 32);

    /* S <-- rS */
    mul256_modm(S1, S1, r);
    
    /* split S to S1 and S2, s.t. S = (S2<<126 | S1), bl(S1) = 126, bl(S2) = 127 */
    S2[0] = (S1[2] >> 14) | ((S1[3] & 0x3FFF) << 42);
    S2[1] = (S1[3] >> 14) | ((S1[4] & 0x3FFF) << 42);
    S2[2] = S1[4] >> 14;
    S1[2] &= 0x3FFF;
    S1[3] = 0;
    S1[4] = 0;

	/* 
	 * [S]B + (-R) + h(-A) =? 0 ==>
	 * [rS]B + [r](-R) + [t](-A) =? 0 ==>
	 * [S1]B + [S2]([2^126]B) + [r](-R) + [t](-A) =? 0 
	*/
	// ge25519_quadruple_scalarmult_vartime(&sumBRA, &R, &A, r, t, S1, S2);
	ge25519_quadruple_scalarmult_samePre_vartime(&sumBRA, &R, &A, r, t, S1, S2);

    return ge25519_is_neutral_vartime(&sumBRA) ? 0: -1;
}



int
ED25519_FN(ed25519_sign_open_hgcd) (const unsigned char *m, size_t mlen, const ed25519_public_key pk, const ed25519_signature RS) {
	ge25519 ALIGN(16) R, A, sumBRA;
	hash_512bits hash;
	bignum256modm hram, S1, S2={0}, r, t;
    int r_isneg;
	int t_isneg;

	if ((RS[63] & 224))
		return -1;

	/* hram <-- H(R,A,m) */
	ed25519_hram(hash, RS, pk, m, mlen);
	expand256_modm(hram, hash, 64);

	/* compute r and t s.t. rh = t mod el */
    // curve25519_half_size_scalar_vartime(r, t, hram, &r_isneg, &t_isneg);
    curve25519_half_size_scalar_vartime_hgcd(r, t, hram, &r_isneg, &t_isneg);
    
    /* unpacking (-R) */
    if (!ge25519_unpack_negative_vartime(&R, RS))
		return -1;

	/*
		[S]B + (-R) + h(-A) =? 0 ==> [rS]B + [r](-R) + [t](-A) =? 0
		------------------------------------------------------------
		r_isneg | t_isneg | [rs]B + [r](-R) + [t](-A) =? 0
		------------------------------------------------------------
		   0    |    0    | [|r| * s]B + [|r|](-R) + [|t|](-A)
		------------------------------------------------------------
		   0    |    1    | [|r| * s]B + [|r|](-R) - [|t|](-A) ==>
		        |         | [|r| * s]B + [|r|](-R) + [|t|](A)
		------------------------------------------------------------
		   1    |    0    | [-|r| * s]B + [-|r|](-R) + [|t|](-A) ==> 
		        |         | [|r| * s] B + [|r|](-R) + [|t|](A)
		------------------------------------------------------------
		   1    |    1    | [-|r| * s]B + [-|r|](-R) + [-|t|](-A) ==> 
		        |         | [|r| * s] B + [|r|](-R) + [|t|](-A)
		------------------------------------------------------------
	*/
	/* unpacking (-A) and adjust the sign based on the previous table */
    if (r_isneg == t_isneg){
		if (!ge25519_unpack_negative_vartime(&A, pk))
			return -1;
	}else{
		if (!ge25519_unpack_positive_vartime(&A, pk))
			return -1;
	}
    
    /* S */
	expand256_modm(S1, RS + 32, 32);

    /* S <-- rS */
    mul256_modm(S1, S1, r);
    
    /* split S to S1 and S2, s.t. S = (S2<<126 | S1), bl(S1) = 126, bl(S2) = 127 */
    S2[0] = (S1[2] >> 14) | ((S1[3] & 0x3FFF) << 42);
    S2[1] = (S1[3] >> 14) | ((S1[4] & 0x3FFF) << 42);
    S2[2] = S1[4] >> 14;
    S1[2] &= 0x3FFF;
    S1[3] = 0;
    S1[4] = 0;

	/* 
	 * [S]B + (-R) + h(-A) =? 0 ==>
	 * [rS]B + [r](-R) + [t](-A) =? 0 ==>
	 * [S1]B + [S2]([2^126]B) + [r](-R) + [t](-A) =? 0 
	*/
	ge25519_quadruple_scalarmult_vartime(&sumBRA, &R, &A, r, t, S1, S2);

    return ge25519_is_neutral_vartime(&sumBRA) ? 0: -1;
}