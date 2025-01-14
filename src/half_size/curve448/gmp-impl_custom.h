/*
 function declaration from internal-used gmp-imp.h
 */
#include "gmp.h"


#define MPN_NORMALIZE(DST, NLIMBS) \
  do {									\
    while ((NLIMBS) > 0)						\
      {									\
	if ((DST)[(NLIMBS) - 1] != 0)					\
	  break;							\
	(NLIMBS)--;							\
      }									\
  } while (0)


#define MP_PTR_SWAP(x, y)						\
  do {									\
    mp_ptr __mp_ptr_swap__tmp = (x);					\
    (x) = (y);								\
    (y) = __mp_ptr_swap__tmp;						\
  } while (0)


struct hgcd_matrix
{
  mp_size_t alloc;		/* for sanity checking only */
  mp_size_t n;
  mp_ptr p[2][2];
};



#define mpn_hgcd_matrix_init __MPN(hgcd_matrix_init)
__GMP_DECLSPEC void mpn_hgcd_matrix_init (struct hgcd_matrix *, mp_size_t, mp_ptr);

#define mpn_hgcd_step __MPN(hgcd_step)
__GMP_DECLSPEC mp_size_t mpn_hgcd_step (mp_size_t, mp_ptr, mp_ptr, mp_size_t, struct hgcd_matrix *, mp_ptr);

#define mpn_hgcd __MPN(hgcd)
__GMP_DECLSPEC mp_size_t mpn_hgcd (mp_ptr, mp_ptr, mp_size_t, struct hgcd_matrix *, mp_ptr);


inline void mp_limb_t_2_mpz(mpz_t out, mp_limb_t* inp, mp_size_t n_limbs){
    // size_t n_limbs = sizeof(inp) / sizeof(mp_limb_t);
// Import multi-limb value into mpz_t
    mpz_import(out, n_limbs, -1, sizeof(mp_limb_t), 0, 0, inp);
    // gmp_printf("%ZX\n", out);
}

// /* custom hgcd */
// inline mp_size_t
// mpn_hgcd_custom (mp_ptr ap, mp_ptr bp, mp_size_t n, mp_size_t s,
// 	  struct hgcd_matrix *M, mp_ptr tp)
// {
//     // mp_size_t s = n/2 + 1;
//     mp_size_t nn;
//     int success = 0;

//     for (;;)
//     {
//       nn = mpn_hgcd_step (n, ap, bp, s, M, tp);
//       if (!nn)
// 	    return success ? n : 0;

//       n = nn;
//       success = 1;
//     }
// }


