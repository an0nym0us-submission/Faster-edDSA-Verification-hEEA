#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <openssl/rand.h>
#include "test-ticks.h"
#include "../src/half_size/curve448/curve448_hEEA_vartime.h"
#include "../src/half_size/curve448/curve448_hEEA_div_vartime.h"
#include "../src/half_size/curve448/curve448_hgcd_vartime.h"
#include "../src/half_size/curve448/curve448_reduce_basis_vartime.h"

#define number_of_samples 10000
#define number_of_rounds 10

struct benchmark_result{
	uint64_t best;
	uint64_t median;
	double average;
};

/* u64 to mpz */
void u64_2_mpz(mpz_t out, const uint64_t* in, size_t l) 
{
	mpz_import(out, l, -1, sizeof(uint64_t), 0, 0, in);
}

/* mpz to u64 */
void mpz_2_u64(uint64_t *out, const mpz_t in, size_t l) {
    size_t count,i;
    mpz_export(out, &count, -1, sizeof(uint64_t), 0, 0, in);
    if (count < l){
        for(i=count; i<l; i++){
            out[i] = 0;
        };
    };
};


void rand_mpz(mpz_t b_mpz, gmp_randstate_t state, const mpz_t L){
	for(;;){
		mpz_urandomb(b_mpz, state, 447);
		if(mpz_cmp_ui(b_mpz, 0))
			break;
	}
	
	mpz_mod(b_mpz, b_mpz, L);
}

/* check if c0 = c1 * b mod el */
int check_correctness(uint64_t *c0, uint64_t *c1, const mpz_t b_mpz, const mpz_t L){
	int c0_isneg, c1_isneg;

	mpz_t c0_mpz, c1_mpz, c1_mul_b;
	mpz_inits(c0_mpz, c1_mpz, c1_mul_b, NULL);

    
	c0_isneg = c0[3] >> 63;
	if (c0_isneg){
		c0[0] = ~c0[0];
		c0[1] = ~c0[1];
		c0[2] = ~c0[2];
		c0[3] = ~c0[3];
	
		if (++c0[0] == 0) {
			if (++c0[1] == 0) {
					if (++c0[2] == 0)
						c0[3]++;
				}	
		}
	}
	u64_2_mpz(c0_mpz, c0, 4);

	c1_isneg = c1[3] >> 63;
	if (c1_isneg){
		c1[0] = ~c1[0];
		c1[1] = ~c1[1];
		c1[2] = ~c1[2];
		c1[3] = ~c1[3];
	
		if (++c1[0] == 0) {
			if (++c1[1] == 0) {
					if (++c1[2] == 0)
						c1[3]++;
				}	
		}
	}
	u64_2_mpz(c1_mpz, c1, 4);

	mpz_mul(c1_mul_b, c1_mpz, b_mpz);
	mpz_mod(c1_mul_b, c1_mul_b, L);

	if (c0_isneg != c1_isneg){
		mpz_sub(c1_mul_b, L, c1_mul_b);
	}

	/* check if c1*b = c0 mod el*/
	return mpz_cmp(c1_mul_b, c0_mpz) ? 0:1;	

}

static int cmp_int64(const void *v1, const void *v2)
{
	int64_t x1, x2;

	x1 = *(const int64_t *)v1;
	x2 = *(const int64_t *)v2;
	if (x1 < x2) {
		return -1;
	} else if (x1 == x2) {
		return 0;
	} else {
		return 1;
	}
}

int test_instance(int test_count){

	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, time(NULL));

	mpz_t L;
	mpz_init(L);
	mpz_set_str(L, "181709681073901722637330951972001133588410340171829515070372549795146003961539585716195755291692375963310293709091662304773755859649779", 10);

	mpz_t b_mpz;
	mpz_init(b_mpz);
	uint64_t b[test_count][7], c0[4], c1[4];
	int c0_size, c1_size;

	uint64_t t_begin;
	
	uint64_t t[test_count];
	uint64_t total_t = 0;
	struct benchmark_result benchmark_reduce_basis; 
	struct benchmark_result benchmark_hEEA; 
	struct benchmark_result benchmark_hEEA_div; 
	struct benchmark_result benchmark_gmp_hgcd; 


	/* pick a random b, s.t. 0 < b < L */
	for(size_t j=0; j<test_count; j++)
	{
		rand_mpz(b_mpz, state, L);
		mpz_2_u64(b[j], b_mpz, 7);
	}

	/* Test hEEA_div */
	total_t = 0;
	for(size_t j=0; j<test_count; j++)
	{
		c0[0] = c0[1] = c0[2] = c0[3] = 0ul;
		c1[0] = c1[1] = c1[2] = c1[3] = 0ul;
    	
		t_begin = get_ticks();
		for(int i=0; i<number_of_rounds;i++)
    		curve448_hEEA_div_vartime(c0, c1, b[j]);

		t[j] = get_ticks() - t_begin;
		total_t += t[j]; 
		u64_2_mpz(b_mpz, b[j], 7);
		if (!check_correctness(c0, c1, b_mpz, L)){
			fprintf(stderr, "ERR: wrong reduction result using `curve25519_hEEA_div_vartime`\n");
			exit(EXIT_FAILURE);
		}
	}
	qsort(t, test_count, sizeof(uint64_t), cmp_int64);
	benchmark_hEEA_div.best =  t[0]/number_of_rounds;
	benchmark_hEEA_div.median = t[test_count/2]/number_of_rounds;
	benchmark_hEEA_div.average =  total_t/(number_of_rounds*(double)test_count);


	/* Test reduce_basis */
	total_t = 0;
	for(size_t j=0; j<test_count; j++)
	{
		c0[0] = c0[1] = c0[2] = c0[3] = 0ul;
		c1[0] = c1[1] = c1[2] = c1[3] = 0ul;

		t_begin = get_ticks();
		for(int i=0; i<number_of_rounds;i++)
    		curve448_reduce_basis_vartime(c0, c1, b[j]);

		t[j] = get_ticks() - t_begin;
		total_t += t[j]; 
		u64_2_mpz(b_mpz, b[j], 7);
		if (!check_correctness(c0, c1, b_mpz, L)){
			fprintf(stderr, "ERR: wrong reduction result using `curve25519_reduce_basis_vartime`\n");
			exit(EXIT_FAILURE);
		}
	}
	qsort(t, test_count, sizeof(uint64_t), cmp_int64);
	benchmark_reduce_basis.best =  t[0]/number_of_rounds;
	benchmark_reduce_basis.median = t[test_count/2]/number_of_rounds;
	benchmark_reduce_basis.average =  total_t/(number_of_rounds*(double)test_count);

	/* Test hEEA */
	total_t = 0;
	for(size_t j=0; j<test_count; j++)
	{
		c0[0] = c0[1] = c0[2] = c0[3] = 0ul;
		c1[0] = c1[1] = c1[2] = c1[3] = 0ul;
    	
		t_begin = get_ticks();
		for(int i=0; i<number_of_rounds;i++)
    		curve448_hEEA_vartime(c0, c1, b[j]);

		t[j] = get_ticks() - t_begin;
		total_t += t[j]; 
		u64_2_mpz(b_mpz, b[j], 7);
		if (!check_correctness(c0, c1, b_mpz, L)){
			fprintf(stderr, "ERR: wrong reduction result using `curve25519_hEEA_vartime`\n");
			exit(EXIT_FAILURE);
		}
		
	}
	qsort(t, test_count, sizeof(uint64_t), cmp_int64);
	benchmark_hEEA.best =  t[0]/number_of_rounds;
	benchmark_hEEA.median = t[test_count/2]/number_of_rounds;
	benchmark_hEEA.average =  total_t/(number_of_rounds*(double)test_count);


	/* Test hgcd */
	total_t = 0;
	for(size_t j=0; j<test_count; j++)
	{
		c0[0] = c0[1] = c0[2] = c0[3] = 0ul;
		c1[0] = c1[1] = c1[2] = c1[3] = 0ul;
    	
		t_begin = get_ticks();
		for(int i=0; i<number_of_rounds;i++)
    		curve448_hgcd_vartime(c0, c1, b[j], &c0_size, &c1_size);

		t[j] = get_ticks() - t_begin;
		total_t += t[j]; 
		u64_2_mpz(b_mpz, b[j], 7);
		if (!check_correctness(c0, c1, b_mpz, L)){
			fprintf(stderr, "ERR: wrong reduction result using `curve25519_hgcd_vartime`\n");
			exit(EXIT_FAILURE);
		}
	}
	qsort(t, test_count, sizeof(uint64_t), cmp_int64);
	benchmark_gmp_hgcd.best =  t[0]/number_of_rounds;
	benchmark_gmp_hgcd.median = t[test_count/2]/number_of_rounds;
	benchmark_gmp_hgcd.average =  total_t/(number_of_rounds*(double)test_count);

	printf("───────────────────────────────────────────────────────────────────────\n");
	printf("Time (ticks)|   hEEA_div   |    hEEA_q    | Speed up     | Improvement\n");
	printf("───────────────────────────────────────────────────────────────────────\n");
	printf("Best        | %-12lu | %-12lu | %-12.4f | %.2f %%\n", benchmark_hEEA_div.best, benchmark_hEEA.best, (double)benchmark_hEEA_div.best/(double)benchmark_hEEA.best,((double)benchmark_hEEA_div.best - (double)benchmark_hEEA.best)/((double)benchmark_hEEA_div.best) * 100);
	printf("Median      | %-12lu | %-12lu | %-12.4f | %.2f %%\n", benchmark_hEEA_div.median, benchmark_hEEA.median, (double)benchmark_hEEA_div.median/(double)benchmark_hEEA.median,((double)benchmark_hEEA_div.median - (double)benchmark_hEEA.median)/((double)benchmark_hEEA_div.median) * 100);
	printf("Average     | %-12.2f | %-12.2f | %-12.4f | %.2f %%\n", benchmark_hEEA_div.average, benchmark_hEEA.average, (double)benchmark_hEEA_div.average/(double)benchmark_hEEA.average,((double)benchmark_hEEA_div.average - (double)benchmark_hEEA.average)/((double)benchmark_hEEA_div.average) * 100);
	printf("───────────────────────────────────────────────────────────────────────\n");

	printf("───────────────────────────────────────────────────────────────────────\n");
	printf("Time (ticks)| reduce_basis |    hEEA_q    | Speed up     | Improvement\n");
	printf("───────────────────────────────────────────────────────────────────────\n");
	printf("Best        | %-12lu | %-12lu | %-12.4f | %.2f %%\n", benchmark_reduce_basis.best, benchmark_hEEA.best, (double)benchmark_reduce_basis.best/(double)benchmark_hEEA.best,((double)benchmark_reduce_basis.best - (double)benchmark_hEEA.best)/((double)benchmark_reduce_basis.best) * 100);
	printf("Median      | %-12lu | %-12lu | %-12.4f | %.2f %%\n", benchmark_reduce_basis.median, benchmark_hEEA.median, (double)benchmark_reduce_basis.median/(double)benchmark_hEEA.median,((double)benchmark_reduce_basis.median - (double)benchmark_hEEA.median)/((double)benchmark_reduce_basis.median) * 100);
	printf("Average     | %-12.2f | %-12.2f | %-12.4f | %.2f %%\n", benchmark_reduce_basis.average, benchmark_hEEA.average, (double)benchmark_reduce_basis.average/(double)benchmark_hEEA.average,((double)benchmark_reduce_basis.average - (double)benchmark_hEEA.average)/((double)benchmark_reduce_basis.average) * 100);
	printf("───────────────────────────────────────────────────────────────────────\n");

	printf("───────────────────────────────────────────────────────────────────────\n");
	printf("Time (ticks)|   GMP_hgcd   |    hEEA_q    | Speed up     | Improvement\n");
	printf("───────────────────────────────────────────────────────────────────────\n");
	printf("Best        | %-12lu | %-12lu | %-12.4f | %.2f %%\n", benchmark_gmp_hgcd.best, benchmark_hEEA.best, (double)benchmark_gmp_hgcd.best/(double)benchmark_hEEA.best,((double)benchmark_gmp_hgcd.best - (double)benchmark_hEEA.best)/((double)benchmark_gmp_hgcd.best) * 100);
	printf("Median      | %-12lu | %-12lu | %-12.4f | %.2f %%\n", benchmark_gmp_hgcd.median, benchmark_hEEA.median, (double)benchmark_gmp_hgcd.median/(double)benchmark_hEEA.median,((double)benchmark_gmp_hgcd.median - (double)benchmark_hEEA.median)/((double)benchmark_gmp_hgcd.median) * 100);
	printf("Average     | %-12.2f | %-12.2f | %-12.4f | %.2f %%\n", benchmark_gmp_hgcd.average, benchmark_hEEA.average, (double)benchmark_gmp_hgcd.average/(double)benchmark_hEEA.average,((double)benchmark_gmp_hgcd.average - (double)benchmark_hEEA.average)/((double)benchmark_gmp_hgcd.average) * 100);
	printf("───────────────────────────────────────────────────────────────────────\n");

    return 0;
}


int main(){
	
	printf("Benchmark of half-size-scalars for Ed448:\n");
	printf("Number of samples = %i \n", number_of_samples);
	printf("Number of rounds = %i \n", number_of_rounds);
	test_instance(number_of_samples);
	printf("Done!\n");
		
}