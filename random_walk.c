#define _POSIX_C_SOURCE 200809L

/*
 * Random non-backtracking 2-isogeny walk over F_{p^2}, using FLINT 3.x.
 *
 * The program selects one of two FLINT backends at runtime:
 *   - fq_nmod when p fits in one machine word;
 *   - fq for arbitrary-size p represented by fmpz.
 *
 * Canonical usage:
 *   ./random_walk (--p p | --pbits n)
 *       [--j re im] [--steps length] [--seed seed]
 *
 * With --pbits n, p is the largest prime below 2^n satisfying p = 3
 * (mod 4). The legacy positional form remains accepted:
 *
 *   ./random_walk p [j_re j_im [length]] [--seed seed]
 *
 * The default starting j-invariant is 1728. This default is permitted only
 * when p = 3 (mod 4), so that j = 1728 is supersingular. The default walk
 * length is 2*ceil(log_2(p)), which equals twice the bit length of an odd
 * prime p.
 *
 * Field elements are represented as re + im*a in
 *
 *   F_{p^2} = F_p[a]/(a^2 - q),
 *
 * where q is a deterministic quadratic nonresidue modulo p. In particular,
 * q = -1 when p = 3 (mod 4), giving the basis a^2 + 1.
 *
 * Compile:
 *   gcc -O3 -march=native -DNDEBUG -std=gnu11 -Wall -Wextra \
 *       random_walk.c \
 *       -lflint -lmpfr -lgmp -o random_walk
 *
 * This is research software. It is not constant-time, and its small
 * SplitMix64 generator is not a cryptographic random-number generator.
 */

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/nmod_poly.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>
#include <flint/fq_poly_factor.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_poly_factor.h>

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef enum {
    BACKEND_FQ_NMOD,
    BACKEND_FQ
} backend_kind_t;

typedef struct {
    backend_kind_t kind;
    fmpz_t characteristic;
    union {
        fq_nmod_ctx_struct word[1];
        fq_ctx_struct big[1];
    } context;
} field_context_t;

typedef struct {
    union {
        struct {
            ulong re, im;
        } word;
        struct {
            fmpz re, im;
        } big;
    } value;
} jkey_t;

enum {
    PHI2_C1488,
    PHI2_C162000,
    PHI2_C40773375,
    PHI2_C8748000000,
    PHI2_C157464000000000,
    PHI2_NCONSTANTS
};

typedef struct {
    union {
        fq_nmod_struct word[PHI2_NCONSTANTS];
        fq_struct big[PHI2_NCONSTANTS];
    } value;
} phi2_constants_t;

typedef struct {
    uint64_t state;
} random_state_t;

static backend_kind_t active_backend;
static ulong field_characteristic_word;
static fmpz_t field_characteristic;

static void die(const char *message)
{
    fprintf(stderr, "error: %s\n", message);
    exit(EXIT_FAILURE);
}

static void usage(FILE *stream, const char *program)
{
    fprintf(stream,
            "usage: %s (--p p | --pbits n) [--j re im]\n"
            "          [--steps length] [--seed seed]\n"
            "legacy: %s p [j_re j_im [length]] [--seed seed]\n",
            program, program);
}

static void parse_characteristic(fmpz_t p, const char *s)
{
    if (fmpz_set_str(p, s, 10) != 0 || fmpz_cmp_ui(p, 3) < 0 ||
        fmpz_is_even(p))
        die("p must be an odd decimal prime");
    if (!fmpz_is_prime(p))
        die("p must be prime");
}

static ulong parse_characteristic_bits(const char *s)
{
    char *end;
    unsigned long bits;

    errno = 0;
    bits = strtoul(s, &end, 10);
    if (errno || *end != '\0' || bits < 2)
        die("pbits must be a decimal integer at least 2");
    return (ulong)bits;
}

static void largest_prime_3mod4_below_power_of_two(
    fmpz_t p, ulong bits)
{
    fmpz_one(p);
    fmpz_mul_2exp(p, p, bits);
    fmpz_sub_ui(p, p, 1);
    while (!fmpz_is_prime(p))
        fmpz_sub_ui(p, p, 4);
}

static size_t parse_steps(const char *s)
{
    char *end;
    unsigned long long value;

    errno = 0;
    value = strtoull(s, &end, 10);
    if (s[0] == '-' || errno || *end != '\0' ||
        value > (unsigned long long)SIZE_MAX)
        die("steps must be a decimal integer in [0,SIZE_MAX]");
    return (size_t)value;
}

static uint64_t parse_seed(const char *s)
{
    char *end;
    unsigned long long value;

    errno = 0;
    value = strtoull(s, &end, 10);
    if (s[0] == '-' || errno || *end != '\0')
        die("seed must be a decimal integer in [0,2^64)");
    return (uint64_t)value;
}

static void parse_field_coefficient(fmpz_t value, const char *s)
{
    if (fmpz_set_str(value, s, 10) != 0 || fmpz_sgn(value) < 0 ||
        fmpz_cmp(value, field_characteristic) >= 0)
        die("j-invariant coefficients must be decimal integers in [0,p)");
}

static void key_init(jkey_t *key)
{
    if (active_backend == BACKEND_FQ) {
        fmpz_init(&key->value.big.re);
        fmpz_init(&key->value.big.im);
    } else {
        key->value.word.re = 0;
        key->value.word.im = 0;
    }
}

static void key_clear(jkey_t *key)
{
    if (active_backend == BACKEND_FQ) {
        fmpz_clear(&key->value.big.re);
        fmpz_clear(&key->value.big.im);
    }
}

static void key_set(jkey_t *result, const jkey_t *source)
{
    if (active_backend == BACKEND_FQ_NMOD) {
        result->value.word = source->value.word;
    } else {
        fmpz_set(&result->value.big.re, &source->value.big.re);
        fmpz_set(&result->value.big.im, &source->value.big.im);
    }
}

static int key_equal(const jkey_t *left, const jkey_t *right)
{
    if (active_backend == BACKEND_FQ_NMOD)
        return left->value.word.re == right->value.word.re &&
               left->value.word.im == right->value.word.im;
    return fmpz_equal(&left->value.big.re, &right->value.big.re) &&
           fmpz_equal(&left->value.big.im, &right->value.big.im);
}

static void key_set_strings(jkey_t *key, const char *real,
                            const char *imaginary)
{
    fmpz_t re, im;

    fmpz_init(re);
    fmpz_init(im);
    parse_field_coefficient(re, real);
    parse_field_coefficient(im, imaginary);
    if (active_backend == BACKEND_FQ_NMOD) {
        key->value.word.re = fmpz_get_ui(re);
        key->value.word.im = fmpz_get_ui(im);
    } else {
        fmpz_set(&key->value.big.re, re);
        fmpz_set(&key->value.big.im, im);
    }
    fmpz_clear(im);
    fmpz_clear(re);
}

static void key_set_default_1728(jkey_t *key)
{
    if (active_backend == BACKEND_FQ_NMOD) {
        key->value.word.re = UWORD(1728) % field_characteristic_word;
        key->value.word.im = 0;
    } else {
        fmpz_set_ui(&key->value.big.re, UWORD(1728));
        fmpz_mod(&key->value.big.re, &key->value.big.re,
                 field_characteristic);
        fmpz_zero(&key->value.big.im);
    }
}

static void fprint_key(FILE *stream, const jkey_t *key)
{
    if (active_backend == BACKEND_FQ_NMOD) {
        fprintf(stream, "(%lu, %lu)",
                (unsigned long)key->value.word.re,
                (unsigned long)key->value.word.im);
    } else {
        fprintf(stream, "(");
        fmpz_fprint(stream, &key->value.big.re);
        fprintf(stream, ", ");
        fmpz_fprint(stream, &key->value.big.im);
        fprintf(stream, ")");
    }
}

static void choose_quadratic_nonresidue(fmpz_t q, const fmpz_t p)
{
    if (fmpz_fdiv_ui(p, UWORD(4)) == UWORD(3)) {
        fmpz_sub_ui(q, p, 1);
        return;
    }
    for (fmpz_set_ui(q, 2); fmpz_cmp(q, p) < 0;
         fmpz_add_ui(q, q, 1))
        if (fmpz_jacobi(q, p) == -1)
            return;
    die("could not find a quadratic nonresidue modulo p");
}

static void field_context_init(field_context_t *field, const fmpz_t p,
                               const fmpz_t quadratic_nonresidue)
{
    field->kind = active_backend;
    fmpz_init_set(field->characteristic, p);

    if (field->kind == BACKEND_FQ_NMOD) {
        nmod_poly_t modulus;
        ulong q = fmpz_get_ui(quadratic_nonresidue);

        field_characteristic_word = fmpz_get_ui(p);
        nmod_poly_init(modulus, field_characteristic_word);
        nmod_poly_set_coeff_ui(modulus, 0,
                               field_characteristic_word - q);
        nmod_poly_set_coeff_ui(modulus, 2, UWORD(1));
        fq_nmod_ctx_init_modulus(field->context.word, modulus, "a");
        nmod_poly_clear(modulus);
    } else {
        fmpz_mod_ctx_t prime_context;
        fmpz_mod_poly_t modulus;
        fmpz_t constant;

        fmpz_mod_ctx_init(prime_context, p);
        fmpz_mod_poly_init(modulus, prime_context);
        fmpz_init(constant);
        fmpz_sub(constant, p, quadratic_nonresidue);
        fmpz_mod_poly_set_coeff_fmpz(modulus, 0, constant,
                                     prime_context);
        fmpz_mod_poly_set_coeff_ui(modulus, 2, UWORD(1),
                                   prime_context);
        fq_ctx_init_modulus(field->context.big, modulus, prime_context,
                            "a");
        fmpz_clear(constant);
        fmpz_mod_poly_clear(modulus, prime_context);
        fmpz_mod_ctx_clear(prime_context);
    }
}

static void field_context_clear(field_context_t *field)
{
    if (field->kind == BACKEND_FQ_NMOD)
        fq_nmod_ctx_clear(field->context.word);
    else
        fq_ctx_clear(field->context.big);
    fmpz_clear(field->characteristic);
}

static void key_to_fq_nmod(fq_nmod_t x, const jkey_t *key,
                           nmod_poly_t representative,
                           const fq_nmod_ctx_t ctx)
{
    nmod_poly_zero(representative);
    nmod_poly_set_coeff_ui(representative, 0, key->value.word.re);
    nmod_poly_set_coeff_ui(representative, 1, key->value.word.im);
    fq_nmod_set_nmod_poly(x, representative, ctx);
}

static void fq_nmod_to_key(jkey_t *key, const fq_nmod_t x,
                           nmod_poly_t representative,
                           const fq_nmod_ctx_t ctx)
{
    fq_nmod_get_nmod_poly(representative, x, ctx);
    key->value.word.re = nmod_poly_get_coeff_ui(representative, 0);
    key->value.word.im = nmod_poly_get_coeff_ui(representative, 1);
}

static void key_to_fq(fq_t x, const jkey_t *key,
                      fmpz_poly_t representative,
                      const fq_ctx_t ctx)
{
    fmpz_poly_zero(representative);
    fmpz_poly_set_coeff_fmpz(representative, 0,
                             &key->value.big.re);
    fmpz_poly_set_coeff_fmpz(representative, 1,
                             &key->value.big.im);
    fq_set_fmpz_poly(x, representative, ctx);
}

static void fq_to_key(jkey_t *key, const fq_t x,
                      fmpz_poly_t representative,
                      const fq_ctx_t ctx)
{
    fq_get_fmpz_poly(representative, x, ctx);
    fmpz_poly_get_coeff_fmpz(&key->value.big.re, representative, 0);
    fmpz_poly_get_coeff_fmpz(&key->value.big.im, representative, 1);
}

static void phi2_constants_init(phi2_constants_t *constants,
                                const field_context_t *field)
{
    static const char *decimal[PHI2_NCONSTANTS] = {
        "1488",
        "162000",
        "40773375",
        "8748000000",
        "157464000000000"
    };
    fmpz_t value;
    size_t i;

    fmpz_init(value);
    for (i = 0; i < PHI2_NCONSTANTS; i++) {
        if (fmpz_set_str(value, decimal[i], 10) != 0)
            die("internal error parsing a Phi_2 coefficient");
        if (field->kind == BACKEND_FQ_NMOD) {
            fq_nmod_init(constants->value.word + i,
                         field->context.word);
            fq_nmod_set_ui(
                constants->value.word + i,
                fmpz_fdiv_ui(value, field_characteristic_word),
                field->context.word);
        } else {
            fq_init(constants->value.big + i, field->context.big);
            fmpz_mod(value, value, field_characteristic);
            fq_set_fmpz(constants->value.big + i, value,
                        field->context.big);
        }
    }
    fmpz_clear(value);
}

static void phi2_constants_clear(phi2_constants_t *constants,
                                 const field_context_t *field)
{
    size_t i;

    if (field->kind == BACKEND_FQ_NMOD) {
        for (i = 0; i < PHI2_NCONSTANTS; i++)
            fq_nmod_clear(constants->value.word + i,
                          field->context.word);
    } else {
        for (i = 0; i < PHI2_NCONSTANTS; i++)
            fq_clear(constants->value.big + i, field->context.big);
    }
}

/*
 * Set f(Y) = Phi_2(j,Y):
 *
 * Y^3 + (-j^2 + 1488j - 162000)Y^2
 *     + (1488j^2 + 40773375j + 8748000000)Y
 *     + j^3 - 162000j^2 + 8748000000j - 157464000000000.
 */
static void phi2_in_y_word(
    fq_nmod_poly_t f, const fq_nmod_t j,
    const phi2_constants_t *constants, const fq_nmod_ctx_t ctx)
{
    const fq_nmod_struct *c = constants->value.word;
    fq_nmod_t j2, j3, c0, c1, c2, temporary;

    fq_nmod_init(j2, ctx);
    fq_nmod_init(j3, ctx);
    fq_nmod_init(c0, ctx);
    fq_nmod_init(c1, ctx);
    fq_nmod_init(c2, ctx);
    fq_nmod_init(temporary, ctx);

    fq_nmod_mul(j2, j, j, ctx);
    fq_nmod_mul(j3, j2, j, ctx);

    fq_nmod_mul(c2, j, c + PHI2_C1488, ctx);
    fq_nmod_sub(c2, c2, j2, ctx);
    fq_nmod_sub(c2, c2, c + PHI2_C162000, ctx);

    fq_nmod_mul(c1, j2, c + PHI2_C1488, ctx);
    fq_nmod_mul(temporary, j, c + PHI2_C40773375, ctx);
    fq_nmod_add(c1, c1, temporary, ctx);
    fq_nmod_add(c1, c1, c + PHI2_C8748000000, ctx);

    fq_nmod_set(c0, j3, ctx);
    fq_nmod_mul(temporary, j2, c + PHI2_C162000, ctx);
    fq_nmod_sub(c0, c0, temporary, ctx);
    fq_nmod_mul(temporary, j, c + PHI2_C8748000000, ctx);
    fq_nmod_add(c0, c0, temporary, ctx);
    fq_nmod_sub(c0, c0, c + PHI2_C157464000000000, ctx);

    fq_nmod_poly_zero(f, ctx);
    fq_nmod_poly_set_coeff(f, 0, c0, ctx);
    fq_nmod_poly_set_coeff(f, 1, c1, ctx);
    fq_nmod_poly_set_coeff(f, 2, c2, ctx);
    fq_nmod_one(temporary, ctx);
    fq_nmod_poly_set_coeff(f, 3, temporary, ctx);

    fq_nmod_clear(temporary, ctx);
    fq_nmod_clear(c2, ctx);
    fq_nmod_clear(c1, ctx);
    fq_nmod_clear(c0, ctx);
    fq_nmod_clear(j3, ctx);
    fq_nmod_clear(j2, ctx);
}

static void phi2_in_y_big(
    fq_poly_t f, const fq_t j,
    const phi2_constants_t *constants, const fq_ctx_t ctx)
{
    const fq_struct *c = constants->value.big;
    fq_t j2, j3, c0, c1, c2, temporary;

    fq_init(j2, ctx);
    fq_init(j3, ctx);
    fq_init(c0, ctx);
    fq_init(c1, ctx);
    fq_init(c2, ctx);
    fq_init(temporary, ctx);

    fq_mul(j2, j, j, ctx);
    fq_mul(j3, j2, j, ctx);

    fq_mul(c2, j, c + PHI2_C1488, ctx);
    fq_sub(c2, c2, j2, ctx);
    fq_sub(c2, c2, c + PHI2_C162000, ctx);

    fq_mul(c1, j2, c + PHI2_C1488, ctx);
    fq_mul(temporary, j, c + PHI2_C40773375, ctx);
    fq_add(c1, c1, temporary, ctx);
    fq_add(c1, c1, c + PHI2_C8748000000, ctx);

    fq_set(c0, j3, ctx);
    fq_mul(temporary, j2, c + PHI2_C162000, ctx);
    fq_sub(c0, c0, temporary, ctx);
    fq_mul(temporary, j, c + PHI2_C8748000000, ctx);
    fq_add(c0, c0, temporary, ctx);
    fq_sub(c0, c0, c + PHI2_C157464000000000, ctx);

    fq_poly_zero(f, ctx);
    fq_poly_set_coeff(f, 0, c0, ctx);
    fq_poly_set_coeff(f, 1, c1, ctx);
    fq_poly_set_coeff(f, 2, c2, ctx);
    fq_one(temporary, ctx);
    fq_poly_set_coeff(f, 3, temporary, ctx);

    fq_clear(temporary, ctx);
    fq_clear(c2, ctx);
    fq_clear(c1, ctx);
    fq_clear(c0, ctx);
    fq_clear(j3, ctx);
    fq_clear(j2, ctx);
}

static size_t phi2_edges_word(
    jkey_t edges[3], const jkey_t *point,
    const phi2_constants_t *constants, const fq_nmod_ctx_t ctx)
{
    fq_nmod_poly_t specialization;
    fq_nmod_poly_factor_t roots;
    fq_nmod_t j, constant, root;
    nmod_poly_t representative;
    size_t count = 0;
    slong factor_index;

    fq_nmod_poly_init(specialization, ctx);
    fq_nmod_poly_factor_init(roots, ctx);
    fq_nmod_init(j, ctx);
    fq_nmod_init(constant, ctx);
    fq_nmod_init(root, ctx);
    nmod_poly_init(representative, field_characteristic_word);

    key_to_fq_nmod(j, point, representative, ctx);
    phi2_in_y_word(specialization, j, constants, ctx);
    fq_nmod_poly_roots(roots, specialization, 1, ctx);

    for (factor_index = 0; factor_index < roots->num; factor_index++) {
        slong copy;

        if (fq_nmod_poly_degree(roots->poly + factor_index, ctx) != 1)
            continue;
        fq_nmod_poly_get_coeff(constant, roots->poly + factor_index,
                               0, ctx);
        fq_nmod_neg(root, constant, ctx);
        for (copy = 0; copy < roots->exp[factor_index]; copy++) {
            if (count == 3)
                die("Phi_2 specialization has more than three roots");
            fq_nmod_to_key(edges + count++, root, representative, ctx);
        }
    }

    nmod_poly_clear(representative);
    fq_nmod_clear(root, ctx);
    fq_nmod_clear(constant, ctx);
    fq_nmod_clear(j, ctx);
    fq_nmod_poly_factor_clear(roots, ctx);
    fq_nmod_poly_clear(specialization, ctx);
    return count;
}

static size_t phi2_edges_big(
    jkey_t edges[3], const jkey_t *point,
    const phi2_constants_t *constants, const fq_ctx_t ctx)
{
    fq_poly_t specialization;
    fq_poly_factor_t roots;
    fq_t j, constant, root;
    fmpz_poly_t representative;
    size_t count = 0;
    slong factor_index;

    fq_poly_init(specialization, ctx);
    fq_poly_factor_init(roots, ctx);
    fq_init(j, ctx);
    fq_init(constant, ctx);
    fq_init(root, ctx);
    fmpz_poly_init(representative);

    key_to_fq(j, point, representative, ctx);
    phi2_in_y_big(specialization, j, constants, ctx);
    fq_poly_roots(roots, specialization, 1, ctx);

    for (factor_index = 0; factor_index < roots->num; factor_index++) {
        slong copy;

        if (fq_poly_degree(roots->poly + factor_index, ctx) != 1)
            continue;
        fq_poly_get_coeff(constant, roots->poly + factor_index, 0, ctx);
        fq_neg(root, constant, ctx);
        for (copy = 0; copy < roots->exp[factor_index]; copy++) {
            if (count == 3)
                die("Phi_2 specialization has more than three roots");
            fq_to_key(edges + count++, root, representative, ctx);
        }
    }

    fmpz_poly_clear(representative);
    fq_clear(root, ctx);
    fq_clear(constant, ctx);
    fq_clear(j, ctx);
    fq_poly_factor_clear(roots, ctx);
    fq_poly_clear(specialization, ctx);
    return count;
}

static size_t phi2_edges(jkey_t edges[3], const jkey_t *point,
                         const phi2_constants_t *constants,
                         const field_context_t *field)
{
    if (field->kind == BACKEND_FQ_NMOD)
        return phi2_edges_word(edges, point, constants,
                               field->context.word);
    return phi2_edges_big(edges, point, constants,
                          field->context.big);
}

static uint64_t mix64(uint64_t x)
{
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

static uint64_t random_u64(random_state_t *random_state)
{
    uint64_t z;

    random_state->state += UINT64_C(0x9e3779b97f4a7c15);
    z = random_state->state;
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

static uint64_t random_state_init(random_state_t *random_state,
                                  int seed_supplied, uint64_t supplied_seed)
{
    uint64_t seed = supplied_seed;

    if (!seed_supplied) {
        FILE *source = fopen("/dev/urandom", "rb");

        seed = 0;
        if (source != NULL) {
            if (fread(&seed, sizeof(seed), 1, source) != 1)
                seed = 0;
            fclose(source);
        }
        if (seed == 0) {
            struct timespec realtime;

            if (clock_gettime(CLOCK_REALTIME, &realtime) != 0)
                die("clock_gettime failed while seeding the walk");
            seed = mix64((uint64_t)realtime.tv_sec ^
                         ((uint64_t)realtime.tv_nsec << 21) ^
                         (uint64_t)(uintptr_t)&seed);
        }
    }
    random_state->state = seed;
    return seed;
}

static size_t random_below(random_state_t *random_state, size_t bound)
{
    uint64_t threshold, value;

    if (bound == 0)
        die("invalid random-choice bound");
    threshold = (UINT64_C(0) - (uint64_t)bound) % (uint64_t)bound;
    do {
        value = random_u64(random_state);
    } while (value < threshold);
    return (size_t)(value % (uint64_t)bound);
}

int main(int argc, char **argv)
{
    const char *positional[4];
    const char *p_string = NULL, *pbits_string = NULL;
    const char *j_re_string = NULL, *j_im_string = NULL;
    const char *steps_string = NULL, *seed_string = NULL;
    size_t n_positional = 0, steps = 0, step, i;
    int p_supplied = 0, pbits_supplied = 0, j_supplied = 0;
    int steps_supplied = 0, seed_supplied = 0, help_requested = 0;
    uint64_t supplied_seed = 0, used_seed;
    ulong pbits = 0;
    fmpz_t quadratic_nonresidue;
    field_context_t field;
    phi2_constants_t constants;
    jkey_t current, previous, edges[3];
    random_state_t random_state;
    int have_previous = 0, completed = 1;
    int argument_index;

    fmpz_init(field_characteristic);
    fmpz_init(quadratic_nonresidue);

    for (argument_index = 1; argument_index < argc; argument_index++) {
        if (strcmp(argv[argument_index], "--p") == 0) {
            if (p_supplied || argument_index + 1 == argc)
                die("--p must occur once and have a value");
            p_string = argv[++argument_index];
            p_supplied = 1;
        } else if (strcmp(argv[argument_index], "--pbits") == 0) {
            if (pbits_supplied || argument_index + 1 == argc)
                die("--pbits must occur once and have a value");
            pbits_string = argv[++argument_index];
            pbits_supplied = 1;
        } else if (strcmp(argv[argument_index], "--j") == 0 ||
                   strcmp(argv[argument_index], "--j1") == 0) {
            if (j_supplied || argument_index + 2 >= argc)
                die("--j must occur once and have two values");
            j_re_string = argv[++argument_index];
            j_im_string = argv[++argument_index];
            j_supplied = 1;
        } else if (strcmp(argv[argument_index], "--steps") == 0) {
            if (steps_supplied || argument_index + 1 == argc)
                die("--steps must occur once and have a value");
            steps_string = argv[++argument_index];
            steps_supplied = 1;
        } else if (strcmp(argv[argument_index], "--seed") == 0) {
            if (seed_supplied || argument_index + 1 == argc)
                die("--seed must occur once and have a value");
            seed_string = argv[++argument_index];
            seed_supplied = 1;
        } else if (strcmp(argv[argument_index], "--help") == 0 ||
                   strcmp(argv[argument_index], "-h") == 0) {
            help_requested = 1;
        } else {
            if (argv[argument_index][0] == '-')
                die("unknown command-line option");
            if (n_positional ==
                sizeof(positional) / sizeof(positional[0]))
                die("too many positional arguments");
            positional[n_positional++] = argv[argument_index];
        }
    }

    if (help_requested) {
        usage(stdout, argv[0]);
        fmpz_clear(quadratic_nonresidue);
        fmpz_clear(field_characteristic);
        return EXIT_SUCCESS;
    }
    if (p_supplied && pbits_supplied)
        die("--p and --pbits are mutually exclusive");

    if (n_positional != 0) {
        if (p_supplied || pbits_supplied || j_supplied ||
            steps_supplied)
            die("do not mix legacy positional input with --p, --pbits, "
                "--j, or --steps");
        if (n_positional != 1 && n_positional != 3 &&
            n_positional != 4) {
            usage(stderr, argv[0]);
            die("the positional form needs 1, 3, or 4 arguments");
        }
        p_string = positional[0];
        p_supplied = 1;
        if (n_positional >= 3) {
            j_re_string = positional[1];
            j_im_string = positional[2];
            j_supplied = 1;
        }
        if (n_positional == 4) {
            steps_string = positional[3];
            steps_supplied = 1;
        }
    }
    if (!p_supplied && !pbits_supplied) {
        usage(stderr, argv[0]);
        die("specify the characteristic with --p, --pbits, or positional p");
    }

    if (pbits_supplied) {
        pbits = parse_characteristic_bits(pbits_string);
        largest_prime_3mod4_below_power_of_two(
            field_characteristic, pbits);
    } else {
        parse_characteristic(field_characteristic, p_string);
    }
    active_backend = fmpz_abs_fits_ui(field_characteristic)
                         ? BACKEND_FQ_NMOD
                         : BACKEND_FQ;
    if (!j_supplied &&
        fmpz_fdiv_ui(field_characteristic, UWORD(4)) != UWORD(3))
        die("the default j = 1728 is supersingular only for p = 3 "
            "(mod 4); supply a supersingular j with --j");

    if (steps_supplied) {
        steps = parse_steps(steps_string);
    } else {
        size_t characteristic_bits =
            (size_t)fmpz_bits(field_characteristic);

        if (characteristic_bits > SIZE_MAX / 2)
            die("the default walk length does not fit in size_t");
        steps = 2 * characteristic_bits;
    }
    if (seed_supplied)
        supplied_seed = parse_seed(seed_string);

    choose_quadratic_nonresidue(quadratic_nonresidue,
                                field_characteristic);
    field_context_init(&field, field_characteristic,
                       quadratic_nonresidue);
    key_init(&current);
    key_init(&previous);
    for (i = 0; i < 3; i++)
        key_init(edges + i);
    if (j_supplied)
        key_set_strings(&current, j_re_string, j_im_string);
    else
        key_set_default_1728(&current);
    phi2_constants_init(&constants, &field);
    used_seed = random_state_init(&random_state, seed_supplied,
                                  supplied_seed);

    printf("Field characteristic p = ");
    fmpz_fprint(stdout, field_characteristic);
    if (pbits_supplied)
        printf(" (largest prime below 2^%lu congruent to 3 modulo 4)",
               (unsigned long)pbits);
    printf("\nArithmetic backend = %s\n",
           active_backend == BACKEND_FQ_NMOD
               ? "fq_nmod (single-word characteristic)"
               : "fq (arbitrary-precision characteristic)");
    if (fmpz_fdiv_ui(field_characteristic, UWORD(4)) == UWORD(3)) {
        printf("Field basis: F_p[a]/(a^2 + 1)\n");
    } else {
        printf("Field basis: F_p[a]/(a^2 - ");
        fmpz_fprint(stdout, quadratic_nonresidue);
        printf(")\n");
    }
    printf("Starting j-invariant = ");
    fprint_key(stdout, &current);
    printf("%s\n", j_supplied ? "" : " (default)");
    printf("Walk length = %zu%s\n", steps,
           steps_supplied ? "" : " (default 2*ceil(log_2(p)))");
    printf("Random seed = %llu\n\n", (unsigned long long)used_seed);

    for (step = 0; step < steps; step++) {
        size_t edge_count =
            phi2_edges(edges, &current, &constants, &field);
        size_t candidates[3], candidate_count = 0;
        int removed_back_edge = 0;
        size_t choice;

        for (i = 0; i < edge_count; i++) {
            if (have_previous && !removed_back_edge &&
                key_equal(edges + i, &previous)) {
                removed_back_edge = 1;
                continue;
            }
            candidates[candidate_count++] = i;
        }
        if (candidate_count == 0) {
            fprintf(stderr,
                    "step %zu: Phi_2(j,Y) has no usable root in F_{p^2}; "
                    "is the starting j-invariant supersingular?\n",
                    step + 1);
            completed = 0;
            break;
        }

        choice = candidates[random_below(&random_state, candidate_count)];
        printf("Step %zu: ", step + 1);
        fprint_key(stdout, &current);
        printf("  --[2]-->  ");
        fprint_key(stdout, edges + choice);
        printf("  (%zu non-backtracking edge%s)\n",
               candidate_count, candidate_count == 1 ? "" : "s");

        key_set(&previous, &current);
        key_set(&current, edges + choice);
        have_previous = 1;
    }

    printf("\nFinal j-invariant = ");
    fprint_key(stdout, &current);
    printf("\n");

    phi2_constants_clear(&constants, &field);
    for (i = 0; i < 3; i++)
        key_clear(edges + i);
    key_clear(&previous);
    key_clear(&current);
    field_context_clear(&field);
    fmpz_clear(quadratic_nonresidue);
    fmpz_clear(field_characteristic);
    flint_cleanup();

    return completed ? EXIT_SUCCESS : 2;
}
