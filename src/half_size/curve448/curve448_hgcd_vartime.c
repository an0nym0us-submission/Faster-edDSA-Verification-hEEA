#include "curve448_hgcd_vartime.h"
#include "gmp-impl_custom.h"
#include <stdio.h>


/* return c0 and c1 s.t. c1 * v = c0 mod el */
void
curve448_hgcd_vartime(
	uint64_t *c0, uint64_t *c1, const uint64_t *v, int *c0_size, int *c1_size)
{

    mp_size_t n = 7; // Number of limbs
    
    /* ap = el */
    mp_limb_t ap[] = {
        0x2378c292ab5844f3,
		0x216cc2728dc58f55,
		0xc44edb49aed63690,
		0xffffffff7cca23e9,
		0xffffffffffffffff,
		0xffffffffffffffff,
		0x3fffffffffffffff
    };
    
    /* bp = v */
    mp_limb_t bp[n];
    bp[0] = v[0];
    bp[1] = v[1];
    bp[2] = v[2];
    bp[3] = v[3];
    bp[4] = v[4];
    bp[5] = v[5];
    bp[6] = v[6];

    
    mp_size_t s_matrix = (n+1)/2 + 1;
    struct hgcd_matrix M;        
    mp_limb_t tmp[4 * n];         
    mp_limb_t tmp_matrix[4 * s_matrix];
    mp_size_t n_reduced;

    mp_limb_t rp[n];
    
    mpn_hgcd_matrix_init(&M, n, tmp_matrix);
   
    /*
     (c;d) = M^{-1}(a;b)
     M = [u00 u01]
         [u10 u11]
     size(u**) ~< N/2
     c =  u11 * a - u01 * b
     d = -u10 * a + u00 * b
     size(c) = size(d) > s
     size(abs(c-d)) = s
     */
    n_reduced = mpn_hgcd(ap, bp, n, &M, tmp);

    mpz_t r2, r1, t2, t1, q;
    mpz_inits(r2, r1, t2, t1, q, NULL);

    int c;
    c = mpn_cmp(ap, bp, n_reduced);
    if(c>0){
        mp_limb_t_2_mpz(r2, ap, n);
        mp_limb_t_2_mpz(r1, bp, n);
        mp_limb_t_2_mpz(t2, M.p[0][1], M.n);
        mpz_neg(t2, t2);
        mp_limb_t_2_mpz(t1, M.p[0][0], M.n);
    }else{
        mp_limb_t_2_mpz(r1, ap, n);
        mp_limb_t_2_mpz(r2, bp, n);
        mp_limb_t_2_mpz(t1, M.p[0][1], M.n);
        mpz_neg(t1, t1);
        mp_limb_t_2_mpz(t2, M.p[0][0], M.n);
    }

    for(;;){
    	if (mpz_sizeinbase(r1, 2) <= 223) {
            for(int i=0; i<r1->_mp_size;i++){
                c0[i] = r1->_mp_d[i];
            }
            *c0_size = r1->_mp_size;

            if (t1->_mp_size < 0 ){
                t1->_mp_size *= -1;
                for(int i=0; i< t1->_mp_size;i++){
                    t1->_mp_d[i] = ~ t1->_mp_d[i];
                }
                mpz_add_ui(t1, t1, 1);
            }
            *c1_size = t1->_mp_size;
            for(int i=0; i< *c1_size;i++){
                c1[i] = t1->_mp_d[i];
            }
			mpz_clears(q, r1, r2, t1, t2, NULL);
			return;
		}

        mpz_tdiv_qr(q, r2, r2, r1);
        mpz_swap(r2, r1);
        mpz_submul(t2, q, t1);
        mpz_swap(t2, t1);
    }
}
