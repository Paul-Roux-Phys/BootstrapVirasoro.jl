#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpfr.h>

#define PREC 256

void eval_out_of_place(mpfr_t result, mpfr_t *a, int n, mpfr_t q) {
    mpfr_t power, tmp;
    mpfr_init_set_ui(power, 1, MPFR_RNDN); // q^0
    mpfr_init(tmp);
    mpfr_set_ui(result, 0, MPFR_RNDN);

    for (int i = 0; i < n; i++) {
        mpfr_mul(tmp, a[i], power, MPFR_RNDN);
        mpfr_add(result, result, tmp, MPFR_RNDN);
        mpfr_mul(power, power, q, MPFR_RNDN);
    }

    mpfr_clear(power);
    mpfr_clear(tmp);
}

void eval_in_place(mpfr_t result, mpfr_t *a, int n, mpfr_t q) {
    mpfr_set(result, a[n-1], MPFR_RNDN);
    for (int i = n - 2; i >= 0; i--) {
        mpfr_mul(result, result, q, MPFR_RNDN);     // result *= q
        mpfr_add(result, result, a[i], MPFR_RNDN);  // result += a[i]
    }
}

double benchmark(void (*f)(mpfr_t, mpfr_t*, int, mpfr_t), mpfr_t *a, int n, mpfr_t q, int trials) {
    mpfr_t result;
    mpfr_init(result);
    clock_t start = clock();
    for (int i = 0; i < trials; i++) {
        f(result, a, n, q);
    }
    clock_t end = clock();
    mpfr_clear(result);
    return (double)(end - start) / trials / CLOCKS_PER_SEC;
}

int main() {
    const int n = 100000;
    const int trials = 1000;
    mpfr_t *a = malloc(n * sizeof(mpfr_t));
    for (int i = 0; i < n; i++) {
        mpfr_init2(a[i], PREC);
        mpfr_set_d(a[i], 1.0 / (i + 1), MPFR_RNDN); // arbitrary values
    }

    mpfr_t q;
    mpfr_init2(q, PREC);
    mpfr_set_d(q, 0.5, MPFR_RNDN);

    double t_out = benchmark(eval_out_of_place, a, n, q, trials);
    /* double t_in  = benchmark(eval_in_place, a, n, q, trials); */

    printf("Out-of-place time: %.6f s\n", t_out);
    printf("In-place     time: %.6f s\n", t_in);

    for (int i = 0; i < n; i++) mpfr_clear(a[i]);
    free(a);
    mpfr_clear(q);
    return 0;
}
