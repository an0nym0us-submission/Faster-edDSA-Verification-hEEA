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

#include <gmp.h>
#include "../half_size/curve25519/curve25519_hEEA_vartime.h"
#include "../half_size/curve25519/curve25519_hgcd_vartime.h"

#define n_words_bignum256modm bignum256modm_limb_size
#define n_bytes_in_word_bignum256modm sizeof(bignum256modm_element_t)
#define nails_bits_bignum256modm_2_mpz 8*n_bytes_in_word_bignum256modm - bignum256modm_bits_per_limb
#define bignum256modm_2_mpz(out, in) mpz_import(out, n_words_bignum256modm, -1, n_bytes_in_word_bignum256modm, 0, nails_bits_bignum256modm_2_mpz, in);


/* bignum256modm to uint_64 */
void bignum256modm_2_u64(uint64_t *out, const bignum256modm in){
    out[0] = (in[0]      ) | ((in[1] & 0xFF) << 56);
    out[1] = (in[1] >>  8) | ((in[2] & 0xFFFF) << 48);
    out[2] = (in[2] >> 16) | ((in[3] & 0xFFFFFF) << 40);
    out[3] = (in[3] >> 24) | ((in[4] & 0xFFFFFFFF) << 32);
}

/* uint_64 to bignum256modm */
void u64_2_bignum256modm(bignum256modm out, const uint64_t *in){
    out[0] = in[0] & 0xFFFFFFFFFFFFFF;
    out[1] = (in[0] >> 56) | ((in[1] & 0xFFFFFFFFFFFF) << 8);
    out[2] = (in[1] >> 48) | ((in[2] & 0xFFFFFFFFFF) << 16);
    out[3] = (in[2] >> 40) | ((in[3] & 0xFFFFFFFF) << 24);
    out[4] = in[3] >> 32;
}

/* mpz to bignum256modm */
void mpz_2_bignum256modm(bignum256modm out, const mpz_t in) {
    size_t count,i;
    mpz_export((void *) out, &count, -1, n_bytes_in_word_bignum256modm, 0, nails_bits_bignum256modm_2_mpz, in);
    if (count < n_words_bignum256modm){
        for(i=count; i<n_words_bignum256modm; i++){
            out[i] = 0;
        };
    };
};


/* bignum256modm inverse mod el */
void inv256_modm(bignum256modm out, bignum256modm in){
    static mp_limb_t __m[4] = {
    0x5812631a5cf5d3ed,
    0x14def9dea2f79cd6,
    0x0000000000000000,
    0x1000000000000000
    };

	__mpz_struct _m = {
    	._mp_alloc = 4,
    	._mp_size = 4,
    	._mp_d = __m
	};
	mpz_t m = { _m };
    mpz_t x, y;
    mpz_inits(x,y, NULL);
    bignum256modm_2_mpz(x, in);
    mpz_invert(y, x, m);
    mpz_2_bignum256modm(out, y);
};

/* ge25519_unpack_positive_vartime*/
static int
ge25519_unpack_positive_vartime(ge25519 *r, const unsigned char p[32]) {
	static const unsigned char zero[32] = {0};
	static const bignum25519 one = {1};
	unsigned char parity = p[31] >> 7;
	unsigned char check[32];
	bignum25519 t, root, num, den, d3;

	curve25519_expand(r->y, p);
	curve25519_copy(r->z, one);
	curve25519_square(num, r->y); /* x = y^2 */
	curve25519_mul(den, num, ge25519_ecd); /* den = dy^2 */
	curve25519_sub_reduce(num, num, r->z); /* x = y^1 - 1 */
	curve25519_add(den, den, r->z); /* den = dy^2 + 1 */

	/* Computation of sqrt(num/den) */
	/* 1.: computation of num^((p-5)/8)*den^((7p-35)/8) = (num*den^7)^((p-5)/8) */
	curve25519_square(t, den);
	curve25519_mul(d3, t, den);
	curve25519_square(r->x, d3);
	curve25519_mul(r->x, r->x, den);
	curve25519_mul(r->x, r->x, num);
	curve25519_pow_two252m3(r->x, r->x);

	/* 2. computation of r->x = num * den^3 * (num*den^7)^((p-5)/8) */
	curve25519_mul(r->x, r->x, d3);
	curve25519_mul(r->x, r->x, num);

	/* 3. Check if either of the roots works: */
	curve25519_square(t, r->x);
	curve25519_mul(t, t, den);
	curve25519_sub_reduce(root, t, num);
	curve25519_contract(check, root);
	if (!ed25519_verify(check, zero, 32)) {
		curve25519_add_reduce(t, t, num);
		curve25519_contract(check, t);
		if (!ed25519_verify(check, zero, 32))
			return 0;
		curve25519_mul(r->x, r->x, ge25519_sqrtneg1);
	}

	curve25519_contract(check, r->x);
	if ((check[0] & 1) != parity) {
		curve25519_copy(t, r->x);
		curve25519_neg(r->x, t);
	}
	curve25519_mul(r->t, r->x, r->y);
	return 1;
}

/* half_size_scalar */
void curve25519_half_size_scalar_vartime_hEEA(bignum256modm r, bignum256modm t, bignum256modm v, int *r_negative, int *t_negative){
    uint64_t b[4], c0[4] = {0}, c1[4] = {0};
    bignum256modm_2_u64(b, v);
    curve25519_hEEA_vartime(c0, c1, b);
    *t_negative = c0[1] >> 63;
	if (*t_negative){
		c0[0] = ~c0[0];
		c0[1] = ~c0[1];
		if (++c0[0] == 0) {
			c0[1]++;
		}
	}
	*r_negative = c1[1] >> 63;
	if (*r_negative){
		c1[0] = ~c1[0];
		c1[1] = ~c1[1];
		if (++c1[0] == 0) {
			c1[1]++;
		}
	}
    u64_2_bignum256modm(t, c0);
    u64_2_bignum256modm(r, c1);
}
/* half_size_scalar */
void curve25519_half_size_scalar_vartime_hgcd(bignum256modm r, bignum256modm t, bignum256modm v, int *r_negative, int *t_negative){
    uint64_t b[4], c0[4] = {0}, c1[4] = {0};
    bignum256modm_2_u64(b, v);
	int c0_size, c1_size;
    curve25519_hgcd_vartime_enhance2(c0, c1, b, &c0_size, &c1_size);
    *t_negative = 0;
    // *t_negative = c0[1] >> 63;
	// if (*t_negative){
	// 	c0[0] = ~c0[0];
	// 	c0[1] = ~c0[1];
	// 	if (++c0[0] == 0) {
	// 		c0[1]++;
	// 	}
	// }
	*r_negative = c1[1] >> 63;
	if (*r_negative){
		c1[0] = ~c1[0];
		c1[1] = ~c1[1];
		if (++c1[0] == 0) {
			c1[1]++;
		}
	}
    u64_2_bignum256modm(t, c0);
    u64_2_bignum256modm(r, c1);
}


/* 2^126 B*/
static const ge25519_niels ge25519_niels_sliding_multiples2[32] = {
    {{0x0007d60613037524,0x0006d61f784d4a6b,0x0007a642bb8842b7,0x0005fcd646854d91,0x00047204d08d72fd},{0x00042eb30d4b497f,0x0000d7379990e0e4,0x000045bd147be58c,0x0005821bca849a6c,0x00005468d6201405},{0x000565a9f93267de,0x0001b81ab1d1401e,0x0004638a3b3b3f5e,0x0001a9510af16e79,0x0004599ee919b633}},
    {{0x0003e91d6a54c980,0x0004bea21b31f482,0x00003d89195529dc,0x0002dfb56143562d,0x000105ba38985c82},{0x0007016a267dad09,0x00013456fa6c6691,0x00001c884cbd635b,0x0003284bbadd2ca7,0x000759d087ff9e6a},{0x000380dd9d8a5ddb,0x0003320b0498f7f5,0x0002ad92e05bbd89,0x0005fdb1ff38b519,0x0005461548ae00ae}},
    {{0x00069de9527493ab,0x00023843fe9e9a6d,0x0004c1d43378cd4b,0x0007be63da56c8bb,0x00010bef1d8a6961},{0x00007f5d47d18c38,0x00077b1c7b48f0ef,0x000548150d4d51aa,0x0000715463304c16,0x000476fc511210d5},{0x0000d28d62a1b8c1,0x00074d63f2d24c1d,0x00020ad27a48d66b,0x00063f763b89b401,0x0005e9c111f1bbdd}},
    {{0x000317abc651ecc4,0x00045b5b9aa572a5,0x00048880fb38fff1,0x000220c6b87831b7,0x000212ff1fa04498},{0x00010b75a08d0c28,0x000251ebbf08a51d,0x000092dd2fae8118,0x00045b9acc9784d1,0x000489ad71d18519},{0x000195799f52209f,0x0002cf3e4bb3e0c7,0x00041bc78cf28aa1,0x00032742ca068ade,0x00057df53c1085ab}},
    {{0x0002fb343143df91,0x0002e2abb1386e75,0x000191c95796167c,0x0007f2ba7e4434a5,0x000526642f7b8b9b},{0x00076f0a4835ea15,0x00031c9f24e15f4e,0x0005e394cd1f49e9,0x00016e3cfece5b71,0x0003cf0ec6da7418},{0x0004c6ea23d63436,0x0002a302e70fb4be,0x00053c41d1863e44,0x00003f8129248a4f,0x0004e9b5b1c1b001}},
    {{0x00027aca52792fdc,0x00067768afe14104,0x00063d12e5c3e3dc,0x0006935610ebff2b,0x00077781e77abe0f},{0x0003072292c284bf,0x00078bb80408ebf8,0x0007d1de1fe88c02,0x0003c5d2c9afcb90,0x00055162f70acd57},{0x0006e29b5d6c8a2e,0x00074b59260b9bd2,0x0003639e166c1cdc,0x00030c4ad24700b5,0x0000c15656b2383a}},
    {{0x00018f2ede49bb66,0x0000a207b1e13ca7,0x00067119b1b4baeb,0x0001bed53d974613,0x0005612e6ff2d202},{0x0001adb8c2ab09b0,0x0006d7557752d5b5,0x000784f548c0619a,0x0005481758cadfaa,0x0005a5a0fa21d62a},{0x0007be5c45379ca1,0x0007a467c505a3b2,0x000158765ea786ab,0x00013b80a01a33ce,0x00016c56fb9e6b60}},
    {{0x000161e9086365cd,0x000182b07cf04deb,0x00077c3a7a266aa3,0x0007b033899e4a45,0x000241d6697a0c10},{0x00027e63f3be17be,0x0007ef4bf3b88bec,0x00057a7f2508b146,0x0001c911dd650016,0x00060f451fcab5b6},{0x00014cf1a7efff7b,0x000029b1a9689ec3,0x0007b28d8ee6335a,0x000009986887d911,0x0004bd0b8881a482}},
    {{0x0006dfa4449a9164,0x00011f50c12a0bdc,0x00043bab873e5de4,0x000658115a82b5bd,0x000382c794f610b7},{0x000693039b778df1,0x0002d94bd65af452,0x000309b0550a2bc7,0x0007ade29b0e667d,0x0001d6b2f9986fdc},{0x00066dfbf2263b59,0x00075ae83fe1da13,0x000763238b3c7519,0x00070998a0ddc855,0x0002a8dbf0931987}},
    {{0x0003d68616f7433c,0x00020cc54c27826d,0x0004d292c87c1b0f,0x000666989d4f2384,0x0004f9e0c9864502},{0x0007a5498bea9a69,0x00079f439740e455,0x0000a0afa8789917,0x0007e66a50bcad90,0x0005923d05ff17c0},{0x000444d2573f0ee9,0x000557a0e6136632,0x0007670f80cc81be,0x0003dea55014566a,0x0004eaf009ad2613}},
    {{0x000701b2085f0eef,0x00004c7ff66d74d6,0x0001589809bd857b,0x0001cbdb94929e1e,0x000626a05d50d4e6},{0x0007ca37420e18c0,0x00050499b41e1e5b,0x00034bb46a0d56dd,0x000723611896fd4c,0x0005b43e6275b9f5},{0x00023e730e6fc24f,0x00048023f487e0ad,0x0005b0a3794aad5d,0x0000e9be543ad5e9,0x00071288f1b8759b}},
    {{0x00003a2c0d3e47b0,0x00010b99f812d3c6,0x000144de340416e3,0x0005c777232a1095,0x0006ec17664ba6bb},{0x0003a5eb9692e387,0x0005993120b662c1,0x000563083b7cb17f,0x0005336b11d5826a,0x0002c9254222b62e},{0x000607207d1bef26,0x0001614f23f5a3da,0x00065afe553b9092,0x0005f29cf77b44a5,0x000712b8e3986acd}},
    {{0x00052b3cd8a711fc,0x0002cbc15e394341,0x000502bd0bee5478,0x00008956352ad433,0x000356dd9e4ae62f},{0x00032db0ca1b92db,0x000165cfd2c36995,0x000297e7cb7c8ca5,0x0003af4a66d89b2c,0x0006348575c25672},{0x000241c2ef612f4f,0x0002a3c18e5e514e,0x0005b93a900b723e,0x0002fe80e2df8a59,0x0007124b72e418a4}},
    {{0x000217e9928db2f7,0x0001d0266ff07264,0x00038d5add5f17e8,0x0001016c64db3ab4,0x0004a5b586a1ddc4},{0x000712f87fd2ccf4,0x00065f3ebea8edd6,0x000468b29cf92012,0x0001c38877a5c662,0x00074a03e7701f30},{0x00033aa637f2d184,0x0000724b33292843,0x0006cd05a56cf338,0x0006b9b424b7e50b,0x000607dc2e9a157c}},
    {{0x00063d77d0c4ee14,0x000478393888a371,0x0004b1268457ab8b,0x000236b81bcae93a,0x0000564a1c113a97},{0x000659e007d81be4,0x0003b96b2efe1421,0x00034d3b9a730d4a,0x00061b8d78351117,0x00066f0575336154},{0x000388dcf109d53e,0x000155ad4aada861,0x0007cbcafc95d3ec,0x000366d89592f17d,0x00011585acec385f}},
    {{0x0001f19f78c044db,0x00036806aa6c7814,0x0005aab4501a83d3,0x0007f78ea07edd28,0x0000e19c01cc1554},{0x0002d01e6cef9d4e,0x0004f89180615494,0x000416f56ff73f64,0x0002c007d2795965,0x000333c2319a1cdb},{0x000235ae2f695126,0x0003ebd8668ab7ce,0x0007e8039238345f,0x0007966e3f1c87c3,0x0006ed7c6cd428b0}},
    {{0x0007caaf2802991b,0x000327350ca81f2c,0x0004b75e7ba5f70b,0x00055edd1ee93446,0x00027677fa3055f2},{0x00035f61f09b366e,0x00069a10faa736a0,0x0004ac47c40c2897,0x0006d99e1ca16f0f,0x0005a5c9bca68cd1},{0x0005ad9a8b9ad766,0x00067ec631af9520,0x0005924a8a827b19,0x0003287dc2a82c1d,0x00057e763e2e93a9}},
    {{0x0006440add4bc8ab,0x0002b909e072544f,0x00035c5452101a9c,0x0003c9cb114f6cc9,0x00043dc66163e30f},{0x00039f6d612eb4d0,0x00040e64f4cda063,0x0007eb79583e38bf,0x0004467ed3e4cddf,0x00016b1400b9d9c5},{0x0005cf10e4020b7c,0x0007ade5fb94e10d,0x000516c2f28e2b1e,0x0001c9d8fdc8e816,0x0005968804e08936}},
    {{0x00063cbf6e40f5a6,0x00073ca80b2708b7,0x0003e203a0f20e76,0x0001963fa22cc13f,0x0000b936019169e6},{0x00049d308ef60c17,0x00032de95f834da1,0x000198b0b2291739,0x0007fd8601fa5bdd,0x0007258c5e68b5ba},{0x0001531cf5138d59,0x0006107c3f297133,0x00026fd40a8f8bb9,0x00065f5cc960ca9b,0x00055c1d4a3351f2}},
    {{0x0001f0d8da83fe1b,0x0002cfe400dad1bb,0x0001108ce40a675d,0x0005c27b26a0c0d2,0x00053b5625b6820f},{0x000554db2e2d20fa,0x0002dfedec37aae7,0x0003b6cafe9adf02,0x00020e301113913d,0x0000d6a947583924},{0x0004e930171834f7,0x000607b68a94eac1,0x0005ce0a80aad80a,0x0001d78f46322dc5,0x0006044faf707b74}},
    {{0x0000904c4750a99d,0x000505e677290aa7,0x0003385f486e02a1,0x00079d22025820b4,0x00042659d40d3270},{0x0000baf461600752,0x0001c6d1da2eacea,0x000382b0f311d0d4,0x000296b66f101e99,0x00016ce2a475bd5a},{0x000569521c20fa07,0x00027ece81a385dc,0x00024a2a9e1bba2f,0x000721c37c2861e7,0x00055066e8d2d80d}},
    {{0x0003323e7b9c269d,0x0004425734d4e85e,0x0006cbea9e16dfcd,0x00066cdf51505120,0x000518b5023bf335},{0x00013889e8cab503,0x000005fc4ea33e0e,0x0007afeef2672094,0x00005b71b656170f,0x00047c62feef0e23},{0x000538b15b00d1bc,0x0006ff9518e3ea3f,0x000493aaa12b747b,0x0007fcd076fe8235,0x00060683c3e63bd7}},
    {{0x000196c2846a0574,0x0001212f11701ca7,0x00075d10f195742e,0x0002b84a6ae92c03,0x00000512861f11d9},{0x0007621dd41d8fd8,0x00047688a3c72b13,0x0001d2a477db7236,0x0007e8ce57f835e8,0x0001277ac6f337f3},{0x0006d666646682d4,0x000388644cfbfbdc,0x000079bc569c0a8c,0x0001416ae7a020a3,0x0000d7683c6e01cd}},
    {{0x00021334f8ebdfe0,0x0005ad6612be8534,0x0006de86f9b94953,0x000625167d0b452c,0x000638a96f588500},{0x000422d17e5d0176,0x000652c668b742c7,0x0007a1cbaeab72b4,0x0001bbaadf3bbae1,0x00054300374a5e1d},{0x000374cf693957a1,0x00069c1b7cac8b6e,0x0006b5a5a186fb19,0x0005d2387044b03e,0x000340db2ba2ea96}},
    {{0x0007b4ffa8b7c8bf,0x0007dc2a7949daeb,0x000779e7bb926b2f,0x0000210aacf50dfb,0x00027292ced9b817},{0x000306e805f3f6c7,0x0000e8e5009ea9ac,0x00036925c6d82512,0x0004e40f8acc8013,0x000329c9d10a9043},{0x00059805ce47a565,0x00009db4d26aeb10,0x0005e683a53f56ad,0x000117ae6b79d5af,0x000390a93b3622c5}},
    {{0x000650a504820d54,0x0006c24320fe7852,0x0004c43c6cca5f11,0x0001ab587f87768f,0x0000a597ed6425c1},{0x0003eac3c00a848d,0x0007e0c779c28641,0x0004a987c4b30f33,0x00061493ea8ef71b,0x0004400ffc9e28e1},{0x000433e40b727304,0x00015a51d7853867,0x000467f715763f8a,0x0004dcc1d6617f8f,0x00034a1898faadc8}},
    {{0x0001cdb3f141949e,0x00016f93cbc62aa3,0x0002afcb4b38c6ee,0x0006ee99924162ef,0x00027619519526fe},{0x000752e0533be6d8,0x000515931a89aa6d,0x000172e1e58c27b9,0x000041171660f6cc,0x0006d8a6d74898df},{0x0001afbb604976a8,0x00040076ea3c819b,0x0000d76b8f99c1ea,0x00041b31f68fc44a,0x000040116617167a}},
    {{0x0005080d234ee25c,0x00046e1d79d99702,0x0007af905199af6b,0x0002385badbe3a57,0x00062772d2affd09},{0x0000c4f4fc794eb8,0x0001903a8f1bbb20,0x000579d88659a0ce,0x0005c6cbf04ba18c,0x0007e2d5cb6049e0},{0x00017116425389a7,0x0002511de61a73e6,0x000527fb65b877cf,0x00005044d57099e0,0x0004ecf62ab8fb47}},
    {{0x000064bbcb80e261,0x0002ec42dc4ed58c,0x0005711d3cf21183,0x000496826c2416c9,0x00012e8463c36dff},{0x0001c8736a531682,0x000525932b477ace,0x0003c03f2adc1268,0x0005faad75d0ddff,0x00019d994c407845},{0x00016da604a9839c,0x0005a1d81e84e3d1,0x0002c3b74fe50eff,0x0002b88defcec01c,0x000570271f677c38}},
    {{0x0005dc5a6eecbe50,0x0004c96c97db65f4,0x00062079bfc99388,0x0001b706ead2179c,0x0004fb180f4af230},{0x0005a7a4e013d531,0x00008e9730770c96,0x00036ca7e5027e5a,0x0006a199e6524d97,0x00030a43d672f7d1},{0x0000abe1152a5ebf,0x00041890dd0a1257,0x000791d8fd408922,0x00040e35458eef54,0x0005354e47ac0f8b}},
    {{0x000290854a246e5d,0x000156b92e0cc33a,0x00046c23aec7ae80,0x0004f0b86c2d6dbe,0x0005b80e635e8883},{0x00048b301987243c,0x0001374239ad024c,0x0005fb8f64af25aa,0x00048637d226b7e2,0x0002a16b07183262},{0x0005176b03972c7f,0x00040c054dfca2c6,0x0001004115856d50,0x00058c336e4df7f7,0x0000afcba6e89af5}},
    {{0x00053ccd3bf9092f,0x00058bd5e4088660,0x00023146143db175,0x0007b675361d5590,0x0000e8644d95fb23},{0x0005c0e23dde23d2,0x0003152dd99f0b83,0x0002a9b3a7073f1b,0x000640543689625f,0x000388a8e62af42c},{0x0005e7a5f60f452b,0x0002dbd33004ba75,0x0004165797d8a89c,0x00030b6e5272a429,0x00079f06033f05a1}}
};


// /* computes [s1]p1 + [s2]p2 + [s3]basepoint */
// static void 
// ge25519_triple_scalarmult_vartime(ge25519 *r, const ge25519 *p1, const ge25519 *p2, const bignum256modm s1, const bignum256modm s2, const bignum256modm s3) {
// 	signed char slide1[256], slide2[256], slide3[256];
// 	ge25519_pniels pre1[S1_TABLE_SIZE], pre2[S1_TABLE_SIZE];
// 	ge25519 d1, d2;
// 	ge25519_p1p1 t;
// 	int32_t i;

// 	contract256_slidingwindow_modm(slide1, s1, S1_SWINDOWSIZE);
// 	contract256_slidingwindow_modm(slide2, s2, S1_SWINDOWSIZE);
// 	contract256_slidingwindow_modm(slide3, s3, S2_SWINDOWSIZE);
    

// 	ge25519_double(&d1, p1);
// 	ge25519_full_to_pniels(pre1, p1);
// 	ge25519_double(&d2, p2);
// 	ge25519_full_to_pniels(pre2, p2);
// 	for (i = 0; i < S1_TABLE_SIZE - 1; i++){
// 		ge25519_pnielsadd(&pre1[i+1], &d1, &pre1[i]);
// 		ge25519_pnielsadd(&pre2[i+1], &d2, &pre2[i]);
//     }

// 	/* set neutral */
// 	memset(r, 0, sizeof(ge25519));
// 	r->y[0] = 1;
// 	r->z[0] = 1;

// 	i = 255;
// 	while ((i >= 0) && !(slide1[i] | slide2[i] | slide3[i]))
// 		i--;

// 	for (; i >= 0; i--) {
// 		ge25519_double_p1p1(&t, r);
// 		if (slide1[i]) {
// 			ge25519_p1p1_to_full(r, &t);
//             ge25519_pnielsadd_p1p1(&t, r, &pre1[abs(slide1[i]) / 2], (unsigned char)slide1[i] >> 7);
// 		}

// 		if (slide2[i]) {
// 			ge25519_p1p1_to_full(r, &t);
// 			ge25519_pnielsadd_p1p1(&t, r, &pre2[abs(slide2[i]) / 2], (unsigned char)slide2[i] >> 7);
// 		}

// 		if (slide3[i]) {
// 			ge25519_p1p1_to_full(r, &t);
// 			ge25519_nielsadd2_p1p1(&t, r, &ge25519_niels_sliding_multiples[abs(slide3[i]) / 2], (unsigned char)slide3[i] >> 7);
// 		}

// 		ge25519_p1p1_to_partial(r, &t);
// 	}
// }

/* computes [s1]p1 + [s2]p2 + [s3]basepoint + [s4 * 2^126] basepoint*/
static void 
ge25519_quadruple_scalarmult_vartime(ge25519 *r, const ge25519 *p1, const ge25519 *p2, const bignum256modm s1, const bignum256modm s2, const bignum256modm s3, const bignum256modm s4) {
	signed char slide1[256], slide2[256], slide3[256], slide4[256];
	ge25519_pniels pre1[S1_TABLE_SIZE], pre2[S1_TABLE_SIZE];
	ge25519 d1, d2;
	ge25519_p1p1 t;
	int32_t i;

	contract256_slidingwindow_modm(slide1, s1, S1_SWINDOWSIZE);
	contract256_slidingwindow_modm(slide2, s2, S1_SWINDOWSIZE);
	contract256_slidingwindow_modm(slide3, s3, S2_SWINDOWSIZE);
	contract256_slidingwindow_modm(slide4, s4, S2_SWINDOWSIZE);
    
	ge25519_double(&d1, p1);
	ge25519_full_to_pniels(pre1, p1);
	ge25519_double(&d2, p2);
	ge25519_full_to_pniels(pre2, p2);
	for (i = 0; i < S1_TABLE_SIZE - 1; i++){
		ge25519_pnielsadd(&pre1[i+1], &d1, &pre1[i]);
		ge25519_pnielsadd(&pre2[i+1], &d2, &pre2[i]);
    }

	/* set neutral */
	memset(r, 0, sizeof(ge25519));
	r->y[0] = 1;
	r->z[0] = 1;

	i = 255;
	while ((i >= 0) && !(slide1[i] | slide2[i] | slide3[i] | slide4[i]))
		i--;

	for (; i >= 0; i--) {
		ge25519_double_p1p1(&t, r);

		if (slide1[i]) {
			ge25519_p1p1_to_full(r, &t);
            ge25519_pnielsadd_p1p1(&t, r, &pre1[abs(slide1[i]) / 2], (unsigned char)slide1[i] >> 7);
		}

		if (slide2[i]) {
			ge25519_p1p1_to_full(r, &t);
			ge25519_pnielsadd_p1p1(&t, r, &pre2[abs(slide2[i]) / 2], (unsigned char)slide2[i] >> 7);
		}

		if (slide3[i]) {
			ge25519_p1p1_to_full(r, &t);
			ge25519_nielsadd2_p1p1(&t, r, &ge25519_niels_sliding_multiples[abs(slide3[i]) / 2], (unsigned char)slide3[i] >> 7);
		}

		if (slide4[i]) {
			ge25519_p1p1_to_full(r, &t);
			ge25519_nielsadd2_p1p1(&t, r, &ge25519_niels_sliding_multiples2[abs(slide4[i]) / 2], (unsigned char)slide4[i] >> 7);
		}

		ge25519_p1p1_to_partial(r, &t);
	}
}

/* B_par = [2^126] basepoint */
static const ge25519 ge25519_B_par = {
    {0x00022c567d23ea24,0x0001008c10a1cb3c,0x00004fbd2c79d16a,0x0007c2a2c1ffa66d,0x0001f13202c95083},
    {0x0002025c90275f48,0x0007d6ab88ef15a8,0x0007f5ffe8021421,0x0005bf790884f3fe,0x00066336d356c381},
    {0x0000000000000001,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000},
    {0x00018564b1c35b4d,0x0007e436bbed8835,0x0000253b45153643,0x0003c29806a11a66,0x00058f25481ba473}
	};
	
/* computes [s1]p1 + [s2]p2 + [s3]basepoint + [s4 * 2^126] basepoint*/
static void 
ge25519_quadruple_scalarmult_samePre_vartime(ge25519 *r, const ge25519 *p1, const ge25519 *p2, const bignum256modm s1, const bignum256modm s2, const bignum256modm s3, const bignum256modm s4) {
	signed char slide1[256], slide2[256], slide3[256], slide4[256];
	ge25519_pniels pre1[S1_TABLE_SIZE], pre2[S1_TABLE_SIZE], pre4[S1_TABLE_SIZE];
	ge25519 d1, d2, d4;
	ge25519_p1p1 t;
	int32_t i;

	contract256_slidingwindow_modm(slide1, s1, S1_SWINDOWSIZE);
	contract256_slidingwindow_modm(slide2, s2, S1_SWINDOWSIZE);
	contract256_slidingwindow_modm(slide3, s3, S2_SWINDOWSIZE);
	contract256_slidingwindow_modm(slide4, s4, S1_SWINDOWSIZE);
    
	ge25519_double(&d1, p1);
	ge25519_full_to_pniels(pre1, p1);
	ge25519_double(&d2, p2);
	ge25519_full_to_pniels(pre2, p2);
	ge25519_double(&d4, &ge25519_B_par);
	ge25519_full_to_pniels(pre4, &ge25519_B_par);



	for (i = 0; i < S1_TABLE_SIZE - 1; i++){
		ge25519_pnielsadd(&pre1[i+1], &d1, &pre1[i]);
		ge25519_pnielsadd(&pre2[i+1], &d2, &pre2[i]);
		ge25519_pnielsadd(&pre4[i+1], &d4, &pre4[i]);
    }

	/* set neutral */
	memset(r, 0, sizeof(ge25519));
	r->y[0] = 1;
	r->z[0] = 1;

	i = 255;
	while ((i >= 0) && !(slide1[i] | slide2[i] | slide3[i] | slide4[i]))
		i--;

	for (; i >= 0; i--) {
		ge25519_double_p1p1(&t, r);

		if (slide1[i]) {
			ge25519_p1p1_to_full(r, &t);
            ge25519_pnielsadd_p1p1(&t, r, &pre1[abs(slide1[i]) / 2], (unsigned char)slide1[i] >> 7);
		}

		if (slide2[i]) {
			ge25519_p1p1_to_full(r, &t);
			ge25519_pnielsadd_p1p1(&t, r, &pre2[abs(slide2[i]) / 2], (unsigned char)slide2[i] >> 7);
		}

		if (slide3[i]) {
			ge25519_p1p1_to_full(r, &t);
			ge25519_nielsadd2_p1p1(&t, r, &ge25519_niels_sliding_multiples[abs(slide3[i]) / 2], (unsigned char)slide3[i] >> 7);
		}

		// if (slide4[i]) {
		// 	ge25519_p1p1_to_full(r, &t);
		// 	ge25519_nielsadd2_p1p1(&t, r, &ge25519_niels_sliding_multiples2[abs(slide4[i]) / 2], (unsigned char)slide4[i] >> 7);
		// }
		if (slide4[i]) {
			ge25519_p1p1_to_full(r, &t);
			ge25519_pnielsadd_p1p1(&t, r, &pre4[abs(slide4[i]) / 2], (unsigned char)slide4[i] >> 7);
		}

		ge25519_p1p1_to_partial(r, &t);
	}
}