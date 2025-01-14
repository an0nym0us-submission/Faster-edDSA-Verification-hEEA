#include "curve448_hEEA_div_vartime.h"
#include "gmp-impl_custom.h"

/* return c0 and c1 s.t. c1 * v = c0 mod el */
void
curve448_hEEA_div_vartime(
	uint64_t *c0, uint64_t *c1, const uint64_t *v)
{

    mp_size_t n = 7; // Number of limbs

    mpz_t r2, r1, t2, t1, q;
    mpz_inits(r2, r1, t2, t1, q, NULL);
    mpz_set_str(r2, "181709681073901722637330951972001133588410340171829515070372549795146003961539585716195755291692375963310293709091662304773755859649779", 10);
    /* bp = v */
    mp_limb_t bp[n];
    bp[0] = v[0];
    bp[1] = v[1];
    bp[2] = v[2];
    bp[3] = v[3];
    bp[4] = v[4];
    bp[5] = v[5];
    bp[6] = v[6];
    
    mp_limb_t_2_mpz(r1, bp, n);
    mpz_set_ui(t1, 1);

    for(;;){
    	if (mpz_sizeinbase(r1, 2) <= 223) {
            for(int i=0; i<r1->_mp_size;i++){
                c0[i] = r1->_mp_d[i];
            }
            // *c0_size = r1->_mp_size;

            if (t1->_mp_size < 0 ){
                t1->_mp_size *= -1;
                for(int i=0; i< t1->_mp_size;i++){
                    t1->_mp_d[i] = ~ t1->_mp_d[i];
                }
                mpz_add_ui(t1, t1, 1);
            }
            // *c1_size = t1->_mp_size;
            for(int i=0; i< t1->_mp_size;i++){
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
