#define _POSIX_C_SOURCE 200809L

/*
 * Interleaved meet-in-the-middle search between supersingular
 * j-invariants. With zero or one input j-invariant, the target is the
 * p-power Frobenius conjugate. Rerandomization defaults to enabled with zero
 * or one supplied j-invariant and disabled with two supplied j-invariants;
 * the default can be overridden explicitly.
 *
 * Field: F_{p^2} = F_p[a]/(a^2 - q), where p is an arbitrary-size input
 * prime and q is a deterministic quadratic nonresidue modulo p. For
 * p = 3 (mod 4), q = -1, preserving the representation
 * F_p[a]/(a^2 + 1). An element is entered and printed as (real, imag),
 * meaning real + imag*a.
 *
 * The program selects one of two FLINT backends at runtime:
 *   - fq_nmod for p values that fit in one machine word;
 *   - fq for larger p values represented by fmpz.
 *
 * Canonical usage:
 *   ./ordered_prime_search_20260725_064545 [N [phi_directory]] [options]
 *
 * Options:
 *   --p p         field characteristic (default: 4294967311)
 *   --pbits n     instead choose the largest prime p < 2^n with p = 3
 *                 (mod 4); mutually exclusive with --p
 *   --B B         smoothness bound, 2 <= B <= 10000
 *                 (default: floor(N^(1/4)))
 *   --N N         middle smooth-path degree bound
 *                 (default: floor(cuberoot(floor(p/2))))
 *   --phi-dir d   modular-polynomial directory (default: mod_pols)
 *   --j1 re im    use the supplied starting j-invariant
 *   --j2 re im    use this target; requires --j1
 *   --rerandomization on|off
 *                 override the mode-dependent rerandomization default
 *   --threads n   use n FLINT worker threads, 1 <= n <= 256
 *   --multipoint-batch n
 *                 cap one multipoint batch at n source points; n = 0
 *                 means no cap (default: no cap)
 *   --force-multipoint
 *                 use a FLINT multipoint tree for every nonempty
 *                 ordered-search batch (default: off)
 *   --rerandomization-jobs n
 *                 test up to n rerandomization endpoints concurrently,
 *                 1 <= n <= 256 (default: 1)
 *   --seed s      seed automatic j-generation with the 64-bit integer s
 *
 * Legacy positional forms remain accepted:
 *   ./ordered_prime_search_20260725_064545 p B N phi_directory [options]
 *   ./ordered_prime_search_20260725_064545 p j1_re j1_im B N \
 *       phi_directory [options]
 *   ./ordered_prime_search_20260725_064545 p j1_re j1_im j2_re j2_im B N \
 *       phi_directory [options]
 *
 * In automatic mode, p must be 3 modulo 4, so j = 1728 is supersingular.
 * Starting from 1728, the program takes a random non-backtracking 2-isogeny
 * walk of length 2*ceil(log2(p)). Phi_2 is reused from the already-loaded
 * smooth modular polynomials.
 *
 * Let ell' be the smallest prime larger than B. For every prime ell <= B,
 * phi_directory must contain Phi_<ell>.txt or phi_j_<ell>.txt. When
 * rerandomization is enabled, a failed initial search also requires ell'.
 * Two sparse formats are accepted:
 *
 *   x_exponent y_exponent integer_coefficient
 *   [x_exponent,y_exponent] integer_coefficient
 *
 * The bracketed format stores one symmetric half; the swapped monomial is
 * inserted automatically. Coefficients may be arbitrary signed integers.
 *
 * Compile with FLINT 3.x:
 *   gcc -O3 -march=native -DNDEBUG -std=gnu11 -Wall -Wextra \
 *       -Wno-unused-parameter \
 *       ordered_prime_search_20260725_064545.c \
 *       -pthread -lflint -lmpfr -lgmp \
 *       -o ordered_prime_search_20260725_064545
 *
 * Performance design:
 *   - Every Phi_ell is read exactly once. Its coefficient polynomials in X
 *     are reduced modulo p and retained as fq_nmod_poly or fq_poly objects.
 *   - The balls are enumerated in canonical prime-factor order. First all
 *     admissible non-backtracking 2-power extensions are generated, then
 *     all 3-power extensions of every existing endpoint, and so on through
 *     the primes ell <= B. This represents each smooth cyclic degree in one
 *     fixed prime order and removes the Dijkstra heap.
 *   - At an ell-stage, endpoints of different current degrees are evaluated
 *     together with the same Phi_ell, subject to the individual condition
 *     degree*ell <= floor(sqrt(N)). Consecutive ell-power layers remain
 *     sequential because one layer supplies the points for the next.
 *   - The two balls alternate after every batch. Newly discovered vertices
 *     are checked immediately against the other ball.
 *   - By default, a complete eligible ell-layer forms one batch. The
 *     optional --multipoint-batch n cap trades batching efficiency for
 *     lower temporary memory and more frequent intersection checks.
 *   - By default, smaller batches use Horner evaluation, since a multipoint
 *     tree costs more than it saves at that size. --force-multipoint
 *     disables this crossover for ordered-search batches.
 *   - FLINT's worker pool evaluates independent coefficient polynomials and
 *     factors independent specializations in parallel. Graph updates remain
 *     serial, so the shortest-path state needs no locks.
 *   - The first intersection is a valid degree-bounded path, but prime-major
 *     order does not promise that it has minimum total degree.
 *   - If j is omitted, it is generated once by a random non-backtracking
 *     2-isogeny walk of length 2*ceil(log2(p)) from j = 1728. A supplied
 *     --seed makes this generation reproducible.
 *   - If rerandomization is enabled after a failed Frobenius-target search,
 *     non-backtracking ell'-walks from j are enumerated breadth first until
 *     an endpoint j' admits a B-smooth path to (j')^p. The conjugate
 *     rerandomization walk is appended in reverse, so a depth-r result has
 *     total degree ell'^(2r) times the middle degree.
 *   - Rerandomization endpoints can be searched in concurrent batches with
 *     --rerandomization-jobs. With one job, the search retains the usual
 *     --threads FLINT parallelism. With multiple jobs, each independent
 *     search is single-threaded, avoiding nested use of FLINT's shared worker
 *     pool while parallelizing the coarser candidate-search workload.
 *   - If rerandomization is explicitly enabled with two supplied
 *     j-invariants, only j1 is rerandomized. A depth-r endpoint j' is
 *     searched against the fixed j2, giving total degree ell'^r times the
 *     middle degree.
 *   - The final report includes wall-clock timings and peak resident memory.
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
#include <flint/thread_support.h>
#include <flint/ulong_extras.h>

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <pthread.h>
#include <stdint.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <time.h>

#define NO_INDEX ((size_t)-1)
/* Retained only by the unused reference Dijkstra implementation below. */
#define MULTIPOINT_BATCH_SIZE 256
#define MULTIPOINT_THRESHOLD 24
#define MULTIPOINT_MIN_DEGREE 12
#define DEFAULT_CHARACTERISTIC "4294967311"
#define DEFAULT_PHI_DIRECTORY "mod_pols"
#define MAX_SMOOTHNESS_BOUND 10000U
#define MAX_RERANDOMIZATION_JOBS 256

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

static backend_kind_t active_backend;
static ulong field_characteristic_word;
static fmpz_t field_characteristic;

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

/* Phi(X,Y) = sum_y coefficient_x[y](X) Y^y. */
typedef struct {
    unsigned ell;
    unsigned max_x_degree;
    size_t n_coefficients;
    union {
        fq_nmod_poly_struct *word;
        fq_poly_struct *big;
    } coefficient_x;
} modular_poly_t;

typedef struct {
    jkey_t key;
    fmpz_t degree;
    size_t parent;
    unsigned parent_ell;
    size_t heap_pos;
    int settled;
} search_node_t;

typedef struct {
    search_node_t *nodes;
    size_t length, alloc;
    size_t *slots; /* node index + 1; zero means empty */
    size_t table_size;
} search_graph_t;

typedef struct {
    size_t *nodes;
    size_t length, alloc;
} min_heap_t;

typedef struct {
    uint64_t specializations;
    uint64_t root_finds;
    uint64_t fast_batches;
    uint64_t horner_batches;
    uint64_t skipped_backtracks;
} search_stats_t;

typedef struct {
    search_graph_t *graph;
    min_heap_t heap;
    size_t *frontier;
    size_t frontier_length, frontier_alloc;
    size_t settled_count;
    search_stats_t stats;
    const char *label;
    double active_seconds;
} ball_enumerator_t;

typedef struct {
    size_t *items;
    size_t length;
    size_t alloc;
} index_vector_t;

typedef struct {
    search_graph_t *graph;
    search_stats_t stats;
    size_t expanded_states;
    const char *label;
    double active_seconds;
} ordered_ball_t;

/* Zero means that a complete eligible prime-power layer is one batch. */
static size_t ordered_multipoint_batch_limit = 0;
/* Nonzero forces every nonempty ordered-search batch through a product tree. */
static int ordered_force_multipoint = 0;

typedef struct {
    const modular_poly_t *phi;
    const fq_nmod_struct *points;
    size_t count;
    fq_nmod_struct *coefficient_values;
    fq_nmod_poly_struct **tree;
    int use_fast;
    const fq_nmod_ctx_struct *ctx;
} word_evaluation_tasks_t;

typedef struct {
    const modular_poly_t *phi;
    size_t count;
    const fq_nmod_struct *coefficient_values;
    fq_nmod_poly_struct *evaluated;
    const fq_nmod_ctx_struct *ctx;
} word_specialization_tasks_t;

typedef struct {
    const fq_nmod_poly_struct *evaluated;
    size_t count;
    size_t task_count;
    size_t root_stride;
    jkey_t *root_keys;
    size_t *root_counts;
    const fq_nmod_ctx_struct *ctx;
} word_root_tasks_t;

typedef struct {
    const modular_poly_t *phi;
    const fq_struct *points;
    size_t count;
    fq_struct *coefficient_values;
    fq_poly_struct **tree;
    int use_fast;
    const fq_ctx_struct *ctx;
} big_evaluation_tasks_t;

typedef struct {
    const modular_poly_t *phi;
    size_t count;
    const fq_struct *coefficient_values;
    fq_poly_struct *evaluated;
    const fq_ctx_struct *ctx;
} big_specialization_tasks_t;

typedef struct {
    const fq_poly_struct *evaluated;
    size_t count;
    size_t task_count;
    size_t root_stride;
    jkey_t *root_keys;
    size_t *root_counts;
    const fq_ctx_struct *ctx;
} big_root_tasks_t;

typedef struct {
    jkey_t key;
    size_t parent;
    size_t depth;
} rerandomization_node_t;

typedef struct {
    rerandomization_node_t *nodes;
    size_t length, alloc;
} rerandomization_tree_t;

typedef struct {
    uint64_t state;
} random_state_t;

static void frobenius_conjugate(jkey_t *result, const jkey_t *j);

static void die(const char *message)
{
    fprintf(stderr, "error: %s\n", message);
    exit(EXIT_FAILURE);
}

static void *xrealloc(void *ptr, size_t count, size_t size)
{
    void *result;

    if (size != 0 && count > SIZE_MAX / size)
        die("allocation size overflow");
    result = realloc(ptr, count * size);
    if (result == NULL && count != 0)
        die("out of memory");
    return result;
}

static double wall_seconds(void)
{
    struct timespec now;

    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0)
        die("clock_gettime failed");
    return (double)now.tv_sec + 1.0e-9 * (double)now.tv_nsec;
}

static double peak_resident_memory_mib(void)
{
    struct rusage usage;

    if (getrusage(RUSAGE_SELF, &usage) != 0)
        return -1.0;
#if defined(__APPLE__)
    return (double)usage.ru_maxrss / (1024.0 * 1024.0);
#else
    return (double)usage.ru_maxrss / 1024.0;
#endif
}

static int key_equal(jkey_t x, jkey_t y)
{
    if (active_backend == BACKEND_FQ_NMOD)
        return x.value.word.re == y.value.word.re &&
               x.value.word.im == y.value.word.im;
    return fmpz_equal(&x.value.big.re, &y.value.big.re) &&
           fmpz_equal(&x.value.big.im, &y.value.big.im);
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

static void key_set_ui(jkey_t *key, ulong real, ulong imaginary)
{
    if (active_backend == BACKEND_FQ_NMOD) {
        key->value.word.re = real;
        key->value.word.im = imaginary;
    } else {
        fmpz_set_ui(&key->value.big.re, real);
        fmpz_set_ui(&key->value.big.im, imaginary);
    }
}

static void key_array_init(jkey_t *keys, size_t count)
{
    size_t i;

    for (i = 0; i < count; i++)
        key_init(keys + i);
}

static void key_array_clear(jkey_t *keys, size_t count)
{
    size_t i;

    for (i = 0; i < count; i++)
        key_clear(keys + i);
}

static uint64_t mix64(uint64_t x)
{
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

static uint64_t key_hash(jkey_t key)
{
    if (active_backend == BACKEND_FQ_NMOD) {
        return mix64((uint64_t)key.value.word.re) ^
               (mix64((uint64_t)key.value.word.im +
                      UINT64_C(0x9e3779b97f4a7c15)) << 1);
    } else {
        uint64_t re_low =
            (uint64_t)fmpz_fdiv_ui(&key.value.big.re, UWORD(4294967291));
        uint64_t im_low =
            (uint64_t)fmpz_fdiv_ui(&key.value.big.im, UWORD(4294967279));
        uint64_t sizes =
            ((uint64_t)fmpz_bits(&key.value.big.re) << 32) ^
            (uint64_t)fmpz_bits(&key.value.big.im);

        return mix64(re_low ^ sizes) ^
               (mix64(im_low + UINT64_C(0x9e3779b97f4a7c15)) << 1);
    }
}

/* representative is reusable scratch storage, avoiding an allocation here. */
static void key_to_fq_nmod(fq_nmod_t x, jkey_t key,
                           nmod_poly_t representative,
                           const fq_nmod_ctx_t ctx)
{
    nmod_poly_zero(representative);
    nmod_poly_set_coeff_ui(representative, 0, key.value.word.re);
    nmod_poly_set_coeff_ui(representative, 1, key.value.word.im);
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

static void key_to_fq(fq_t x, jkey_t key,
                      fmpz_poly_t representative,
                      const fq_ctx_t ctx)
{
    fmpz_poly_zero(representative);
    fmpz_poly_set_coeff_fmpz(representative, 0, &key.value.big.re);
    fmpz_poly_set_coeff_fmpz(representative, 1, &key.value.big.im);
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

static void fprint_key(FILE *stream, jkey_t key)
{
    if (active_backend == BACKEND_FQ_NMOD) {
        fprintf(stream, "(%lu, %lu)", (unsigned long)key.value.word.re,
                (unsigned long)key.value.word.im);
    } else {
        fprintf(stream, "(");
        fmpz_fprint(stream, &key.value.big.re);
        fprintf(stream, ", ");
        fmpz_fprint(stream, &key.value.big.im);
        fprintf(stream, ")");
    }
}

static void print_key(jkey_t key)
{
    fprint_key(stdout, key);
}

/* ----------------------- modular-polynomial input --------------------- */

static void modular_poly_init(modular_poly_t *phi, unsigned ell,
                              const field_context_t *field)
{
    size_t i;

    phi->ell = ell;
    phi->max_x_degree = 0;
    phi->n_coefficients = (size_t)ell + 2;
    if (field->kind == BACKEND_FQ_NMOD) {
        phi->coefficient_x.word =
            xrealloc(NULL, phi->n_coefficients,
                     sizeof(*phi->coefficient_x.word));
        for (i = 0; i < phi->n_coefficients; i++)
            fq_nmod_poly_init(phi->coefficient_x.word + i,
                              field->context.word);
    } else {
        phi->coefficient_x.big =
            xrealloc(NULL, phi->n_coefficients,
                     sizeof(*phi->coefficient_x.big));
        for (i = 0; i < phi->n_coefficients; i++)
            fq_poly_init(phi->coefficient_x.big + i, field->context.big);
    }
}

static void modular_poly_add_term_word(
    modular_poly_t *phi, unsigned x_exp, unsigned y_exp,
    ulong coefficient_mod_p, fq_nmod_t coefficient,
    fq_nmod_t old_coefficient, const fq_nmod_ctx_t ctx)
{
    if (x_exp > phi->ell + 1 || y_exp > phi->ell + 1)
        die("modular-polynomial exponent exceeds ell+1");
    if (coefficient_mod_p == 0)
        return;

    if (x_exp > phi->max_x_degree)
        phi->max_x_degree = x_exp;

    fq_nmod_set_ui(coefficient, coefficient_mod_p, ctx);
    fq_nmod_poly_get_coeff(old_coefficient,
                           phi->coefficient_x.word + y_exp,
                           x_exp, ctx);
    fq_nmod_add(old_coefficient, old_coefficient, coefficient, ctx);
    fq_nmod_poly_set_coeff(phi->coefficient_x.word + y_exp, x_exp,
                           old_coefficient, ctx);
}

static void modular_poly_add_term_big(
    modular_poly_t *phi, unsigned x_exp, unsigned y_exp,
    const fmpz_t coefficient_mod_p, fq_t coefficient,
    fq_t old_coefficient, const fq_ctx_t ctx)
{
    if (x_exp > phi->ell + 1 || y_exp > phi->ell + 1)
        die("modular-polynomial exponent exceeds ell+1");
    if (fmpz_is_zero(coefficient_mod_p))
        return;

    if (x_exp > phi->max_x_degree)
        phi->max_x_degree = x_exp;

    fq_set_fmpz(coefficient, coefficient_mod_p, ctx);
    fq_poly_get_coeff(old_coefficient, phi->coefficient_x.big + y_exp,
                      x_exp, ctx);
    fq_add(old_coefficient, old_coefficient, coefficient, ctx);
    fq_poly_set_coeff(phi->coefficient_x.big + y_exp, x_exp,
                      old_coefficient, ctx);
}

static void modular_poly_load(modular_poly_t *phi, const char *directory,
                              const field_context_t *field)
{
    char path[4096];
    char *line = NULL;
    size_t capacity = 0;
    ssize_t line_length;
    FILE *input;
    unsigned long line_number = 0;
    size_t term_count = 0;
    fmpz_t parsed_coefficient, coefficient_mod_p_big;
    fq_nmod_t coefficient_word, old_coefficient_word;
    fq_t coefficient_big, old_coefficient_big;

    if (snprintf(path, sizeof(path), "%s/Phi_%u.txt", directory, phi->ell)
        >= (int)sizeof(path))
        die("modular-polynomial path is too long");
    input = fopen(path, "r");
    if (input == NULL) {
        if (snprintf(path, sizeof(path), "%s/phi_j_%u.txt", directory,
                     phi->ell) >= (int)sizeof(path))
            die("modular-polynomial path is too long");
        input = fopen(path, "r");
        if (input == NULL) {
            fprintf(stderr, "error: cannot open Phi_%u.txt or phi_j_%u.txt "
                    "in %s: %s\n", phi->ell, phi->ell, directory,
                    strerror(errno));
            exit(EXIT_FAILURE);
        }
    }

    fmpz_init(parsed_coefficient);
    fmpz_init(coefficient_mod_p_big);
    if (field->kind == BACKEND_FQ_NMOD) {
        fq_nmod_init(coefficient_word, field->context.word);
        fq_nmod_init(old_coefficient_word, field->context.word);
    } else {
        fq_init(coefficient_big, field->context.big);
        fq_init(old_coefficient_big, field->context.big);
    }

    while ((line_length = getline(&line, &capacity, input)) >= 0) {
        char *p = line, *end;
        unsigned long x_exp, y_exp;
        int symmetric_half = 0;

        (void)line_length;
        line_number++;
        while (isspace((unsigned char)*p))
            p++;
        if (*p == '\0' || *p == '#')
            continue;

        errno = 0;
        if (*p == '[') {
            symmetric_half = 1;
            p++;
        }
        x_exp = strtoul(p, &end, 10);
        if (errno || end == p)
            goto malformed;
        p = end;
        if (symmetric_half) {
            while (isspace((unsigned char)*p))
                p++;
            if (*p++ != ',')
                goto malformed;
        }
        y_exp = strtoul(p, &end, 10);
        if (errno || end == p)
            goto malformed;
        p = end;
        if (symmetric_half) {
            while (isspace((unsigned char)*p))
                p++;
            if (*p++ != ']')
                goto malformed;
        }
        while (isspace((unsigned char)*p))
            p++;
        if (*p == '\0')
            goto malformed;
        end = p + strlen(p);
        while (end > p && isspace((unsigned char)end[-1]))
            end--;
        *end = '\0';
        if (x_exp > UINT_MAX || y_exp > UINT_MAX)
            goto malformed;
        if (fmpz_set_str(parsed_coefficient, p, 10) != 0)
            goto malformed;

        if (field->kind == BACKEND_FQ_NMOD) {
            ulong coefficient_mod_p =
                fmpz_fdiv_ui(parsed_coefficient,
                             field_characteristic_word);

            modular_poly_add_term_word(
                phi, (unsigned)x_exp, (unsigned)y_exp,
                coefficient_mod_p, coefficient_word,
                old_coefficient_word, field->context.word);
        } else {
            fmpz_mod(coefficient_mod_p_big, parsed_coefficient,
                     field_characteristic);
            modular_poly_add_term_big(
                phi, (unsigned)x_exp, (unsigned)y_exp,
                coefficient_mod_p_big, coefficient_big,
                old_coefficient_big, field->context.big);
        }
        term_count++;
        if (symmetric_half && x_exp != y_exp) {
            if (field->kind == BACKEND_FQ_NMOD) {
                ulong coefficient_mod_p =
                    fmpz_fdiv_ui(parsed_coefficient,
                                 field_characteristic_word);

                modular_poly_add_term_word(
                    phi, (unsigned)y_exp, (unsigned)x_exp,
                    coefficient_mod_p, coefficient_word,
                    old_coefficient_word, field->context.word);
            } else {
                modular_poly_add_term_big(
                    phi, (unsigned)y_exp, (unsigned)x_exp,
                    coefficient_mod_p_big, coefficient_big,
                    old_coefficient_big, field->context.big);
            }
            term_count++;
        }
        continue;

malformed:
        fprintf(stderr, "error: malformed line %lu in %s\n", line_number,
                path);
        exit(EXIT_FAILURE);
    }

    free(line);
    fclose(input);
    if (field->kind == BACKEND_FQ_NMOD) {
        fq_nmod_clear(coefficient_word, field->context.word);
        fq_nmod_clear(old_coefficient_word, field->context.word);
    } else {
        fq_clear(coefficient_big, field->context.big);
        fq_clear(old_coefficient_big, field->context.big);
    }
    fmpz_clear(coefficient_mod_p_big);
    fmpz_clear(parsed_coefficient);
    if (term_count == 0)
        die("empty modular-polynomial file");
}

static void modular_poly_clear(modular_poly_t *phi,
                               const field_context_t *field)
{
    size_t i;

    if (field->kind == BACKEND_FQ_NMOD) {
        for (i = 0; i < phi->n_coefficients; i++)
            fq_nmod_poly_clear(phi->coefficient_x.word + i,
                               field->context.word);
        free(phi->coefficient_x.word);
    } else {
        for (i = 0; i < phi->n_coefficients; i++)
            fq_poly_clear(phi->coefficient_x.big + i,
                          field->context.big);
        free(phi->coefficient_x.big);
    }
}

/* ------------------------ parallel field work ------------------------- */

static int compare_keys(const void *left, const void *right)
{
    const jkey_t *a = left, *b = right;

    if (active_backend == BACKEND_FQ_NMOD) {
        if (a->value.word.re != b->value.word.re)
            return a->value.word.re < b->value.word.re ? -1 : 1;
        if (a->value.word.im != b->value.word.im)
            return a->value.word.im < b->value.word.im ? -1 : 1;
        return 0;
    } else {
        int comparison =
            fmpz_cmp(&a->value.big.re, &b->value.big.re);

        if (comparison != 0)
            return comparison;
        return fmpz_cmp(&a->value.big.im, &b->value.big.im);
    }
}

static void evaluate_coefficient_task_word(slong y, void *argument)
{
    word_evaluation_tasks_t *tasks = argument;
    const fq_nmod_poly_struct *coefficient =
        tasks->phi->coefficient_x.word + y;
    fq_nmod_struct *values = tasks->coefficient_values + (size_t)y *
                            tasks->count;

    if (fq_nmod_poly_is_zero(coefficient, tasks->ctx))
        return;
    if (tasks->use_fast) {
        _fq_nmod_poly_evaluate_fq_nmod_vec_fast_precomp(
            values, coefficient->coeffs, coefficient->length, tasks->tree,
            (slong)tasks->count, tasks->ctx);
    } else {
        _fq_nmod_poly_evaluate_fq_nmod_vec_iter(
            values, coefficient->coeffs, coefficient->length, tasks->points,
            (slong)tasks->count, tasks->ctx);
    }
}

static void build_specialization_task_word(slong point_index, void *argument)
{
    word_specialization_tasks_t *tasks = argument;
    fq_nmod_poly_struct *specialization = tasks->evaluated + point_index;
    size_t y;

    fq_nmod_poly_zero(specialization, tasks->ctx);
    for (y = 0; y < tasks->phi->n_coefficients; y++) {
        const fq_nmod_poly_struct *coefficient =
            tasks->phi->coefficient_x.word + y;

        if (!fq_nmod_poly_is_zero(coefficient, tasks->ctx)) {
            fq_nmod_poly_set_coeff(
                specialization, (slong)y,
                tasks->coefficient_values + y * tasks->count +
                    (size_t)point_index,
                tasks->ctx);
        }
    }
}

static void find_roots_task_word(slong task_index, void *argument)
{
    word_root_tasks_t *tasks = argument;
    fq_nmod_poly_factor_t roots;
    fq_nmod_t constant, root;
    nmod_poly_t representative;
    size_t first = tasks->count * (size_t)task_index / tasks->task_count;
    size_t last = tasks->count * ((size_t)task_index + 1) /
                  tasks->task_count;
    size_t i;

    fq_nmod_poly_factor_init(roots, tasks->ctx);
    fq_nmod_init(constant, tasks->ctx);
    fq_nmod_init(root, tasks->ctx);
    nmod_poly_init(representative, field_characteristic_word);

    for (i = first; i < last; i++) {
        slong factor_index;
        size_t count = 0;
        jkey_t *point_roots = tasks->root_keys + i * tasks->root_stride;

        fq_nmod_poly_roots(roots, tasks->evaluated + i, 0, tasks->ctx);
        for (factor_index = 0; factor_index < roots->num; factor_index++) {
            if (fq_nmod_poly_degree(roots->poly + factor_index,
                                    tasks->ctx) != 1)
                continue;
            if (count == tasks->root_stride)
                break;
            fq_nmod_poly_get_coeff(constant, roots->poly + factor_index, 0,
                                   tasks->ctx);
            fq_nmod_neg(root, constant, tasks->ctx);
            fq_nmod_to_key(point_roots + count++, root, representative,
                           tasks->ctx);
        }
        qsort(point_roots, count, sizeof(*point_roots), compare_keys);
        tasks->root_counts[i] = count;
    }

    nmod_poly_clear(representative);
    fq_nmod_clear(root, tasks->ctx);
    fq_nmod_clear(constant, tasks->ctx);
    fq_nmod_poly_factor_clear(roots, tasks->ctx);
}

/*
 * Specialize one modular polynomial at one j-invariant and return its
 * F_{p^2}-rational roots with multiplicity. Multiplicity matters here:
 * different ell-isogenies can have the same target j-invariant, and the
 * non-backtracking walk must remove exactly one dual edge.
 *
 * This routine handles one walk endpoint at a time, so a multipoint tree
 * would cost more than direct evaluation.
 */
static size_t modular_poly_edges_at_word(
    jkey_t **edges_out, jkey_t point, const modular_poly_t *phi,
    const fq_nmod_ctx_t ctx)
{
    fq_nmod_struct *coefficient_values;
    fq_nmod_struct evaluation_point;
    fq_nmod_poly_struct specialization;
    fq_nmod_poly_factor_t roots;
    fq_nmod_t constant, root;
    word_evaluation_tasks_t evaluation_tasks;
    word_specialization_tasks_t specialization_tasks;
    nmod_poly_t representative;
    size_t y, count = 0;
    slong factor_index;

    coefficient_values =
        xrealloc(NULL, phi->n_coefficients, sizeof(*coefficient_values));
    *edges_out = xrealloc(NULL, (size_t)phi->ell + 1,
                          sizeof(**edges_out));
    key_array_init(*edges_out, (size_t)phi->ell + 1);

    fq_nmod_init(&evaluation_point, ctx);
    fq_nmod_poly_init(&specialization, ctx);
    fq_nmod_poly_factor_init(roots, ctx);
    fq_nmod_init(constant, ctx);
    fq_nmod_init(root, ctx);
    nmod_poly_init(representative, field_characteristic_word);
    key_to_fq_nmod(&evaluation_point, point, representative, ctx);
    for (y = 0; y < phi->n_coefficients; y++)
        fq_nmod_init(coefficient_values + y, ctx);

    evaluation_tasks.phi = phi;
    evaluation_tasks.points = &evaluation_point;
    evaluation_tasks.count = 1;
    evaluation_tasks.coefficient_values = coefficient_values;
    evaluation_tasks.tree = NULL;
    evaluation_tasks.use_fast = 0;
    evaluation_tasks.ctx = ctx;
    for (y = 0; y < phi->n_coefficients; y++)
        evaluate_coefficient_task_word((slong)y, &evaluation_tasks);

    specialization_tasks.phi = phi;
    specialization_tasks.count = 1;
    specialization_tasks.coefficient_values = coefficient_values;
    specialization_tasks.evaluated = &specialization;
    specialization_tasks.ctx = ctx;
    build_specialization_task_word(0, &specialization_tasks);

    fq_nmod_poly_roots(roots, &specialization, 1, ctx);
    for (factor_index = 0; factor_index < roots->num; factor_index++) {
        slong multiplicity, copy;
        if (fq_nmod_poly_degree(roots->poly + factor_index, ctx) != 1)
            continue;
        fq_nmod_poly_get_coeff(constant, roots->poly + factor_index, 0,
                               ctx);
        fq_nmod_neg(root, constant, ctx);
        multiplicity = roots->exp[factor_index];
        for (copy = 0; copy < multiplicity; copy++) {
            if (count == (size_t)phi->ell + 1)
                die("too many roots in a modular-polynomial specialization");
            fq_nmod_to_key(*edges_out + count++, root, representative, ctx);
        }
    }
    qsort(*edges_out, count, sizeof(**edges_out), compare_keys);

    for (y = 0; y < phi->n_coefficients; y++)
        fq_nmod_clear(coefficient_values + y, ctx);
    nmod_poly_clear(representative);
    fq_nmod_clear(root, ctx);
    fq_nmod_clear(constant, ctx);
    fq_nmod_poly_factor_clear(roots, ctx);
    fq_nmod_poly_clear(&specialization, ctx);
    fq_nmod_clear(&evaluation_point, ctx);
    free(coefficient_values);
    return count;
}

static void evaluate_coefficient_task_big(slong y, void *argument)
{
    big_evaluation_tasks_t *tasks = argument;
    const fq_poly_struct *coefficient =
        tasks->phi->coefficient_x.big + y;
    fq_struct *values =
        tasks->coefficient_values + (size_t)y * tasks->count;

    if (fq_poly_is_zero(coefficient, tasks->ctx))
        return;
    if (tasks->use_fast) {
        _fq_poly_evaluate_fq_vec_fast_precomp(
            values, coefficient->coeffs, coefficient->length, tasks->tree,
            (slong)tasks->count, tasks->ctx);
    } else {
        _fq_poly_evaluate_fq_vec_iter(
            values, coefficient->coeffs, coefficient->length, tasks->points,
            (slong)tasks->count, tasks->ctx);
    }
}

static void build_specialization_task_big(slong point_index, void *argument)
{
    big_specialization_tasks_t *tasks = argument;
    fq_poly_struct *specialization = tasks->evaluated + point_index;
    size_t y;

    fq_poly_zero(specialization, tasks->ctx);
    for (y = 0; y < tasks->phi->n_coefficients; y++) {
        const fq_poly_struct *coefficient =
            tasks->phi->coefficient_x.big + y;

        if (!fq_poly_is_zero(coefficient, tasks->ctx)) {
            fq_poly_set_coeff(
                specialization, (slong)y,
                tasks->coefficient_values + y * tasks->count +
                    (size_t)point_index,
                tasks->ctx);
        }
    }
}

static void find_roots_task_big(slong task_index, void *argument)
{
    big_root_tasks_t *tasks = argument;
    fq_poly_factor_t roots;
    fq_t constant, root;
    fmpz_poly_t representative;
    size_t first = tasks->count * (size_t)task_index / tasks->task_count;
    size_t last = tasks->count * ((size_t)task_index + 1) /
                  tasks->task_count;
    size_t i;

    fq_poly_factor_init(roots, tasks->ctx);
    fq_init(constant, tasks->ctx);
    fq_init(root, tasks->ctx);
    fmpz_poly_init(representative);

    for (i = first; i < last; i++) {
        slong factor_index;
        size_t count = 0;
        jkey_t *point_roots = tasks->root_keys + i * tasks->root_stride;

        fq_poly_roots(roots, tasks->evaluated + i, 0, tasks->ctx);
        for (factor_index = 0; factor_index < roots->num; factor_index++) {
            if (fq_poly_degree(roots->poly + factor_index,
                               tasks->ctx) != 1)
                continue;
            if (count == tasks->root_stride)
                break;
            fq_poly_get_coeff(constant, roots->poly + factor_index, 0,
                              tasks->ctx);
            fq_neg(root, constant, tasks->ctx);
            fq_to_key(point_roots + count++, root, representative,
                      tasks->ctx);
        }
        qsort(point_roots, count, sizeof(*point_roots), compare_keys);
        tasks->root_counts[i] = count;
    }

    fmpz_poly_clear(representative);
    fq_clear(root, tasks->ctx);
    fq_clear(constant, tasks->ctx);
    fq_poly_factor_clear(roots, tasks->ctx);
}

static size_t modular_poly_edges_at_big(
    jkey_t **edges_out, jkey_t point, const modular_poly_t *phi,
    const fq_ctx_t ctx)
{
    fq_struct *coefficient_values;
    fq_struct evaluation_point;
    fq_poly_struct specialization;
    fq_poly_factor_t roots;
    fq_t constant, root;
    big_evaluation_tasks_t evaluation_tasks;
    big_specialization_tasks_t specialization_tasks;
    fmpz_poly_t representative;
    size_t y, count = 0;
    slong factor_index;

    coefficient_values =
        xrealloc(NULL, phi->n_coefficients, sizeof(*coefficient_values));
    *edges_out = xrealloc(NULL, (size_t)phi->ell + 1,
                          sizeof(**edges_out));
    key_array_init(*edges_out, (size_t)phi->ell + 1);

    fq_init(&evaluation_point, ctx);
    fq_poly_init(&specialization, ctx);
    fq_poly_factor_init(roots, ctx);
    fq_init(constant, ctx);
    fq_init(root, ctx);
    fmpz_poly_init(representative);
    key_to_fq(&evaluation_point, point, representative, ctx);
    for (y = 0; y < phi->n_coefficients; y++)
        fq_init(coefficient_values + y, ctx);

    evaluation_tasks.phi = phi;
    evaluation_tasks.points = &evaluation_point;
    evaluation_tasks.count = 1;
    evaluation_tasks.coefficient_values = coefficient_values;
    evaluation_tasks.tree = NULL;
    evaluation_tasks.use_fast = 0;
    evaluation_tasks.ctx = ctx;
    for (y = 0; y < phi->n_coefficients; y++)
        evaluate_coefficient_task_big((slong)y, &evaluation_tasks);

    specialization_tasks.phi = phi;
    specialization_tasks.count = 1;
    specialization_tasks.coefficient_values = coefficient_values;
    specialization_tasks.evaluated = &specialization;
    specialization_tasks.ctx = ctx;
    build_specialization_task_big(0, &specialization_tasks);

    fq_poly_roots(roots, &specialization, 1, ctx);
    for (factor_index = 0; factor_index < roots->num; factor_index++) {
        slong multiplicity, copy;

        if (fq_poly_degree(roots->poly + factor_index, ctx) != 1)
            continue;
        fq_poly_get_coeff(constant, roots->poly + factor_index, 0, ctx);
        fq_neg(root, constant, ctx);
        multiplicity = roots->exp[factor_index];
        for (copy = 0; copy < multiplicity; copy++) {
            if (count == (size_t)phi->ell + 1)
                die("too many roots in a modular-polynomial specialization");
            fq_to_key(*edges_out + count++, root, representative, ctx);
        }
    }
    qsort(*edges_out, count, sizeof(**edges_out), compare_keys);

    for (y = 0; y < phi->n_coefficients; y++)
        fq_clear(coefficient_values + y, ctx);
    fmpz_poly_clear(representative);
    fq_clear(root, ctx);
    fq_clear(constant, ctx);
    fq_poly_factor_clear(roots, ctx);
    fq_poly_clear(&specialization, ctx);
    fq_clear(&evaluation_point, ctx);
    free(coefficient_values);
    return count;
}

static size_t modular_poly_edges_at(
    jkey_t **edges_out, jkey_t point, const modular_poly_t *phi,
    const field_context_t *field)
{
    if (field->kind == BACKEND_FQ_NMOD)
        return modular_poly_edges_at_word(edges_out, point, phi,
                                          field->context.word);
    return modular_poly_edges_at_big(edges_out, point, phi,
                                     field->context.big);
}

/* -------------------------- graph hash table -------------------------- */

static void graph_rehash(search_graph_t *graph, size_t new_size)
{
    size_t i;

    free(graph->slots);
    graph->slots = calloc(new_size, sizeof(*graph->slots));
    if (graph->slots == NULL)
        die("out of memory");
    graph->table_size = new_size;

    for (i = 0; i < graph->length; i++) {
        size_t slot = (size_t)key_hash(graph->nodes[i].key) & (new_size - 1);

        while (graph->slots[slot] != 0)
            slot = (slot + 1) & (new_size - 1);
        graph->slots[slot] = i + 1;
    }
}

static void graph_init(search_graph_t *graph)
{
    graph->nodes = NULL;
    graph->length = graph->alloc = 0;
    graph->slots = NULL;
    graph->table_size = 0;
    graph_rehash(graph, 1024);
}

static size_t graph_find(const search_graph_t *graph, jkey_t key)
{
    size_t slot = (size_t)key_hash(key) & (graph->table_size - 1);

    while (graph->slots[slot] != 0) {
        size_t index = graph->slots[slot] - 1;

        if (key_equal(graph->nodes[index].key, key))
            return index;
        slot = (slot + 1) & (graph->table_size - 1);
    }
    return NO_INDEX;
}

static size_t graph_add(search_graph_t *graph, jkey_t key,
                        const fmpz_t degree, size_t parent,
                        unsigned parent_ell)
{
    size_t index, slot;
    search_node_t *node;

    if ((graph->length + 1) * 10 >= graph->table_size * 7)
        graph_rehash(graph, 2 * graph->table_size);
    if (graph->length == graph->alloc) {
        graph->alloc = graph->alloc ? 2 * graph->alloc : 1024;
        graph->nodes = xrealloc(graph->nodes, graph->alloc,
                                sizeof(*graph->nodes));
    }
    index = graph->length++;
    node = graph->nodes + index;
    key_init(&node->key);
    key_set(&node->key, &key);
    fmpz_init_set(node->degree, degree);
    node->parent = parent;
    node->parent_ell = parent_ell;
    node->heap_pos = NO_INDEX;
    node->settled = 0;

    slot = (size_t)key_hash(key) & (graph->table_size - 1);
    while (graph->slots[slot] != 0)
        slot = (slot + 1) & (graph->table_size - 1);
    graph->slots[slot] = index + 1;
    return index;
}

static void graph_clear(search_graph_t *graph)
{
    size_t i;

    for (i = 0; i < graph->length; i++) {
        key_clear(&graph->nodes[i].key);
        fmpz_clear(graph->nodes[i].degree);
    }
    free(graph->nodes);
    free(graph->slots);
}

/*
 * The former Dijkstra implementation is retained below as disabled reference
 * code to keep this derivative easy to compare with find_smooth_path.c.  The
 * executable path in this file uses only the canonical prime-order search.
 */
#if 0
/* -------------------------- indexed min-heap -------------------------- */

static int heap_node_less(const min_heap_t *heap, size_t a, size_t b,
                          const search_graph_t *graph)
{
    size_t node_a = heap->nodes[a], node_b = heap->nodes[b];
    int comparison = fmpz_cmp(graph->nodes[node_a].degree,
                              graph->nodes[node_b].degree);

    return comparison < 0 || (comparison == 0 && node_a < node_b);
}

static void heap_swap(min_heap_t *heap, size_t a, size_t b,
                      search_graph_t *graph)
{
    size_t temporary = heap->nodes[a];

    heap->nodes[a] = heap->nodes[b];
    heap->nodes[b] = temporary;
    graph->nodes[heap->nodes[a]].heap_pos = a;
    graph->nodes[heap->nodes[b]].heap_pos = b;
}

static void heap_insert(min_heap_t *heap, size_t node, search_graph_t *graph)
{
    size_t i, parent;

    if (heap->length == heap->alloc) {
        heap->alloc = heap->alloc ? 2 * heap->alloc : 1024;
        heap->nodes = xrealloc(heap->nodes, heap->alloc,
                               sizeof(*heap->nodes));
    }
    i = heap->length++;
    heap->nodes[i] = node;
    graph->nodes[node].heap_pos = i;
    while (i > 0) {
        parent = (i - 1) / 2;
        if (!heap_node_less(heap, i, parent, graph))
            break;
        heap_swap(heap, i, parent, graph);
        i = parent;
    }
}

static void heap_decrease_key(min_heap_t *heap, size_t node,
                              search_graph_t *graph)
{
    size_t i = graph->nodes[node].heap_pos;

    while (i > 0) {
        size_t parent = (i - 1) / 2;

        if (!heap_node_less(heap, i, parent, graph))
            break;
        heap_swap(heap, i, parent, graph);
        i = parent;
    }
}

static size_t heap_peek(const min_heap_t *heap)
{
    return heap->length ? heap->nodes[0] : NO_INDEX;
}

static size_t heap_pop(min_heap_t *heap, search_graph_t *graph)
{
    size_t result, i, left, right, smallest;

    if (heap->length == 0)
        return NO_INDEX;
    result = heap->nodes[0];
    graph->nodes[result].heap_pos = NO_INDEX;
    heap->length--;
    if (heap->length == 0)
        return result;

    heap->nodes[0] = heap->nodes[heap->length];
    graph->nodes[heap->nodes[0]].heap_pos = 0;
    i = 0;
    for (;;) {
        left = 2 * i + 1;
        right = left + 1;
        smallest = i;
        if (left < heap->length && heap_node_less(heap, left, smallest, graph))
            smallest = left;
        if (right < heap->length &&
            heap_node_less(heap, right, smallest, graph))
            smallest = right;
        if (smallest == i)
            break;
        heap_swap(heap, i, smallest, graph);
        i = smallest;
    }
    return result;
}

static void heap_clear(min_heap_t *heap)
{
    free(heap->nodes);
}

/* --------------------------- bounded search --------------------------- */

static void relax_node(search_graph_t *graph, min_heap_t *heap, jkey_t key,
                       const fmpz_t degree, size_t parent, unsigned ell)
{
    size_t index = graph_find(graph, key);

    if (index == NO_INDEX) {
        index = graph_add(graph, key, degree, parent, ell);
        heap_insert(heap, index, graph);
    } else if (!graph->nodes[index].settled &&
               fmpz_cmp(degree, graph->nodes[index].degree) < 0) {
        fmpz_set(graph->nodes[index].degree, degree);
        graph->nodes[index].parent = parent;
        graph->nodes[index].parent_ell = ell;
        heap_decrease_key(heap, index, graph);
    }
}

static int process_frontier_batch_word(
    search_graph_t *graph, min_heap_t *heap, const size_t *frontier,
    size_t count, const fmpz_t layer_degree, const fmpz_t radius,
    const modular_poly_t *polynomials, size_t n_polynomials,
    const fq_nmod_ctx_t ctx, unsigned n_threads, search_stats_t *stats,
    const search_graph_t *other, size_t *this_meet, size_t *other_meet,
    const atomic_int *cancel_requested, int *cancelled_out)
{
    fq_nmod_struct *points, *coefficient_values;
    fq_nmod_poly_struct *evaluated;
    fq_nmod_poly_struct **tree = NULL;
    jkey_t *root_keys;
    size_t *root_counts;
    nmod_poly_t representative;
    fmpz_t next_degree;
    size_t i, k, maximum_coefficients = 0, maximum_root_stride = 0;
    int build_tree = 0, found = 0;

    if (cancel_requested != NULL &&
        atomic_load_explicit(cancel_requested, memory_order_relaxed)) {
        if (cancelled_out != NULL)
            *cancelled_out = 1;
        return 0;
    }

    /* Avoid all field allocation for terminal Dijkstra layers. */
    fmpz_init(next_degree);
    fmpz_mul_ui(next_degree, layer_degree, polynomials[0].ell);
    if (fmpz_cmp(next_degree, radius) > 0) {
        fmpz_clear(next_degree);
        return 0;
    }

    /* A product tree only pays off when an admissible coefficient has
       enough degree; this keeps small Phi_ell evaluations on Horner. */
    for (k = 0; k < n_polynomials; k++) {
        const modular_poly_t *phi = polynomials + k;

        fmpz_mul_ui(next_degree, layer_degree, phi->ell);
        if (fmpz_cmp(next_degree, radius) > 0)
            continue;
        if (phi->n_coefficients > maximum_coefficients)
            maximum_coefficients = phi->n_coefficients;
        if ((size_t)phi->ell + 1 > maximum_root_stride)
            maximum_root_stride = (size_t)phi->ell + 1;
        if (count >= MULTIPOINT_THRESHOLD &&
            phi->max_x_degree >= MULTIPOINT_MIN_DEGREE)
            build_tree = 1;
    }

    points = xrealloc(NULL, count, sizeof(*points));
    coefficient_values = xrealloc(NULL, count * maximum_coefficients,
                                   sizeof(*coefficient_values));
    evaluated = xrealloc(NULL, count, sizeof(*evaluated));
    root_keys = xrealloc(NULL, count * maximum_root_stride,
                         sizeof(*root_keys));
    root_counts = xrealloc(NULL, count, sizeof(*root_counts));
    nmod_poly_init(representative, field_characteristic_word);
    key_array_init(root_keys, count * maximum_root_stride);
    for (i = 0; i < count; i++) {
        fq_nmod_init(points + i, ctx);
        fq_nmod_poly_init(evaluated + i, ctx);
        key_to_fq_nmod(points + i, graph->nodes[frontier[i]].key,
                       representative, ctx);
    }
    for (i = 0; i < count * maximum_coefficients; i++)
        fq_nmod_init(coefficient_values + i, ctx);
    if (build_tree) {
        tree = _fq_nmod_poly_tree_alloc((slong)count, ctx);
        _fq_nmod_poly_tree_build(tree, points, (slong)count, ctx);
    }

    for (k = 0; k < n_polynomials; k++) {
        const modular_poly_t *phi = polynomials + k;
        word_evaluation_tasks_t evaluation_tasks;
        word_specialization_tasks_t specialization_tasks;
        word_root_tasks_t root_tasks;
        size_t first_new, root_task_count;
        int use_fast;

        if (cancel_requested != NULL &&
            atomic_load_explicit(cancel_requested,
                                 memory_order_relaxed)) {
            if (cancelled_out != NULL)
                *cancelled_out = 1;
            break;
        }
        fmpz_mul_ui(next_degree, layer_degree, phi->ell);
        if (fmpz_cmp(next_degree, radius) > 0)
            continue;

        use_fast = tree != NULL &&
                   phi->max_x_degree >= MULTIPOINT_MIN_DEGREE;
        evaluation_tasks.phi = phi;
        evaluation_tasks.points = points;
        evaluation_tasks.count = count;
        evaluation_tasks.coefficient_values = coefficient_values;
        evaluation_tasks.tree = tree;
        evaluation_tasks.use_fast = use_fast;
        evaluation_tasks.ctx = ctx;
        flint_parallel_do(evaluate_coefficient_task_word, &evaluation_tasks,
                          (slong)phi->n_coefficients, (int)n_threads,
                          FLINT_PARALLEL_DYNAMIC);
        if (use_fast)
            stats->fast_batches++;
        else
            stats->horner_batches++;

        specialization_tasks.phi = phi;
        specialization_tasks.count = count;
        specialization_tasks.coefficient_values = coefficient_values;
        specialization_tasks.evaluated = evaluated;
        specialization_tasks.ctx = ctx;
        flint_parallel_do(build_specialization_task_word,
                          &specialization_tasks,
                          (slong)count, (int)n_threads,
                          FLINT_PARALLEL_UNIFORM);

        stats->specializations += (uint64_t)count;
        stats->root_finds += (uint64_t)count;

        root_task_count = count;
        if (root_task_count > (size_t)n_threads * 4)
            root_task_count = (size_t)n_threads * 4;
        root_tasks.evaluated = evaluated;
        root_tasks.count = count;
        root_tasks.task_count = root_task_count;
        root_tasks.root_stride = (size_t)phi->ell + 1;
        root_tasks.root_keys = root_keys;
        root_tasks.root_counts = root_counts;
        root_tasks.ctx = ctx;
        flint_parallel_do(find_roots_task_word, &root_tasks,
                          (slong)root_task_count, (int)n_threads,
                          FLINT_PARALLEL_DYNAMIC);

        first_new = graph->length;
        for (i = 0; i < count; i++) {
            size_t root_index;

            for (root_index = 0; root_index < root_counts[i]; root_index++) {
                jkey_t root_key =
                    root_keys[i * root_tasks.root_stride + root_index];
                const search_node_t *current = graph->nodes + frontier[i];

                if (current->parent != NO_INDEX &&
                    current->parent_ell == phi->ell &&
                    key_equal(root_key, graph->nodes[current->parent].key)) {
                    stats->skipped_backtracks++;
                    continue;
                }
                relax_node(graph, heap, root_key, next_degree, frontier[i],
                           phi->ell);
            }
        }

        for (i = first_new; i < graph->length; i++) {
            size_t match = graph_find(other, graph->nodes[i].key);

            if (match != NO_INDEX) {
                *this_meet = i;
                *other_meet = match;
                found = 1;
                break;
            }
        }
        if (found)
            break;
    }

    fmpz_clear(next_degree);
    if (tree != NULL)
        _fq_nmod_poly_tree_free(tree, (slong)count, ctx);
    for (i = 0; i < count * maximum_coefficients; i++)
        fq_nmod_clear(coefficient_values + i, ctx);
    for (i = 0; i < count; i++) {
        fq_nmod_poly_clear(evaluated + i, ctx);
        fq_nmod_clear(points + i, ctx);
    }
    nmod_poly_clear(representative);
    key_array_clear(root_keys, count * maximum_root_stride);
    free(root_counts);
    free(root_keys);
    free(evaluated);
    free(coefficient_values);
    free(points);
    return found;
}

static int process_frontier_batch_big(
    search_graph_t *graph, min_heap_t *heap, const size_t *frontier,
    size_t count, const fmpz_t layer_degree, const fmpz_t radius,
    const modular_poly_t *polynomials, size_t n_polynomials,
    const fq_ctx_t ctx, unsigned n_threads, search_stats_t *stats,
    const search_graph_t *other, size_t *this_meet, size_t *other_meet,
    const atomic_int *cancel_requested, int *cancelled_out)
{
    fq_struct *points, *coefficient_values;
    fq_poly_struct *evaluated;
    fq_poly_struct **tree = NULL;
    jkey_t *root_keys;
    size_t *root_counts;
    fmpz_poly_t representative;
    fmpz_t next_degree;
    size_t i, k, maximum_coefficients = 0, maximum_root_stride = 0;
    int build_tree = 0, found = 0;

    if (cancel_requested != NULL &&
        atomic_load_explicit(cancel_requested, memory_order_relaxed)) {
        if (cancelled_out != NULL)
            *cancelled_out = 1;
        return 0;
    }

    fmpz_init(next_degree);
    fmpz_mul_ui(next_degree, layer_degree, polynomials[0].ell);
    if (fmpz_cmp(next_degree, radius) > 0) {
        fmpz_clear(next_degree);
        return 0;
    }

    for (k = 0; k < n_polynomials; k++) {
        const modular_poly_t *phi = polynomials + k;

        fmpz_mul_ui(next_degree, layer_degree, phi->ell);
        if (fmpz_cmp(next_degree, radius) > 0)
            continue;
        if (phi->n_coefficients > maximum_coefficients)
            maximum_coefficients = phi->n_coefficients;
        if ((size_t)phi->ell + 1 > maximum_root_stride)
            maximum_root_stride = (size_t)phi->ell + 1;
        if (count >= MULTIPOINT_THRESHOLD &&
            phi->max_x_degree >= MULTIPOINT_MIN_DEGREE)
            build_tree = 1;
    }

    points = xrealloc(NULL, count, sizeof(*points));
    coefficient_values = xrealloc(NULL, count * maximum_coefficients,
                                   sizeof(*coefficient_values));
    evaluated = xrealloc(NULL, count, sizeof(*evaluated));
    root_keys = xrealloc(NULL, count * maximum_root_stride,
                         sizeof(*root_keys));
    root_counts = xrealloc(NULL, count, sizeof(*root_counts));
    fmpz_poly_init(representative);
    key_array_init(root_keys, count * maximum_root_stride);
    for (i = 0; i < count; i++) {
        fq_init(points + i, ctx);
        fq_poly_init(evaluated + i, ctx);
        key_to_fq(points + i, graph->nodes[frontier[i]].key,
                  representative, ctx);
    }
    for (i = 0; i < count * maximum_coefficients; i++)
        fq_init(coefficient_values + i, ctx);
    if (build_tree) {
        tree = _fq_poly_tree_alloc((slong)count, ctx);
        _fq_poly_tree_build(tree, points, (slong)count, ctx);
    }

    for (k = 0; k < n_polynomials; k++) {
        const modular_poly_t *phi = polynomials + k;
        big_evaluation_tasks_t evaluation_tasks;
        big_specialization_tasks_t specialization_tasks;
        big_root_tasks_t root_tasks;
        size_t first_new, root_task_count;
        int use_fast;

        if (cancel_requested != NULL &&
            atomic_load_explicit(cancel_requested,
                                 memory_order_relaxed)) {
            if (cancelled_out != NULL)
                *cancelled_out = 1;
            break;
        }
        fmpz_mul_ui(next_degree, layer_degree, phi->ell);
        if (fmpz_cmp(next_degree, radius) > 0)
            continue;

        use_fast = tree != NULL &&
                   phi->max_x_degree >= MULTIPOINT_MIN_DEGREE;
        evaluation_tasks.phi = phi;
        evaluation_tasks.points = points;
        evaluation_tasks.count = count;
        evaluation_tasks.coefficient_values = coefficient_values;
        evaluation_tasks.tree = tree;
        evaluation_tasks.use_fast = use_fast;
        evaluation_tasks.ctx = ctx;
        flint_parallel_do(evaluate_coefficient_task_big, &evaluation_tasks,
                          (slong)phi->n_coefficients, (int)n_threads,
                          FLINT_PARALLEL_DYNAMIC);
        if (use_fast)
            stats->fast_batches++;
        else
            stats->horner_batches++;

        specialization_tasks.phi = phi;
        specialization_tasks.count = count;
        specialization_tasks.coefficient_values = coefficient_values;
        specialization_tasks.evaluated = evaluated;
        specialization_tasks.ctx = ctx;
        flint_parallel_do(build_specialization_task_big,
                          &specialization_tasks, (slong)count,
                          (int)n_threads, FLINT_PARALLEL_UNIFORM);

        stats->specializations += (uint64_t)count;
        stats->root_finds += (uint64_t)count;

        root_task_count = count;
        if (root_task_count > (size_t)n_threads * 4)
            root_task_count = (size_t)n_threads * 4;
        root_tasks.evaluated = evaluated;
        root_tasks.count = count;
        root_tasks.task_count = root_task_count;
        root_tasks.root_stride = (size_t)phi->ell + 1;
        root_tasks.root_keys = root_keys;
        root_tasks.root_counts = root_counts;
        root_tasks.ctx = ctx;
        flint_parallel_do(find_roots_task_big, &root_tasks,
                          (slong)root_task_count, (int)n_threads,
                          FLINT_PARALLEL_DYNAMIC);

        first_new = graph->length;
        for (i = 0; i < count; i++) {
            size_t root_index;

            for (root_index = 0; root_index < root_counts[i]; root_index++) {
                jkey_t root_key =
                    root_keys[i * root_tasks.root_stride + root_index];
                const search_node_t *current =
                    graph->nodes + frontier[i];

                if (current->parent != NO_INDEX &&
                    current->parent_ell == phi->ell &&
                    key_equal(root_key,
                              graph->nodes[current->parent].key)) {
                    stats->skipped_backtracks++;
                    continue;
                }
                relax_node(graph, heap, root_key, next_degree, frontier[i],
                           phi->ell);
            }
        }

        for (i = first_new; i < graph->length; i++) {
            size_t match = graph_find(other, graph->nodes[i].key);

            if (match != NO_INDEX) {
                *this_meet = i;
                *other_meet = match;
                found = 1;
                break;
            }
        }
        if (found)
            break;
    }

    fmpz_clear(next_degree);
    if (tree != NULL)
        _fq_poly_tree_free(tree, (slong)count, ctx);
    for (i = 0; i < count * maximum_coefficients; i++)
        fq_clear(coefficient_values + i, ctx);
    for (i = 0; i < count; i++) {
        fq_poly_clear(evaluated + i, ctx);
        fq_clear(points + i, ctx);
    }
    fmpz_poly_clear(representative);
    key_array_clear(root_keys, count * maximum_root_stride);
    free(root_counts);
    free(root_keys);
    free(evaluated);
    free(coefficient_values);
    free(points);
    return found;
}

static int process_frontier_batch(
    search_graph_t *graph, min_heap_t *heap, const size_t *frontier,
    size_t count, const fmpz_t layer_degree, const fmpz_t radius,
    const modular_poly_t *polynomials, size_t n_polynomials,
    const field_context_t *field, unsigned n_threads,
    search_stats_t *stats, const search_graph_t *other,
    size_t *this_meet, size_t *other_meet,
    const atomic_int *cancel_requested, int *cancelled_out)
{
    if (field->kind == BACKEND_FQ_NMOD)
        return process_frontier_batch_word(
            graph, heap, frontier, count, layer_degree, radius,
            polynomials, n_polynomials, field->context.word, n_threads,
            stats, other, this_meet, other_meet, cancel_requested,
            cancelled_out);
    return process_frontier_batch_big(
        graph, heap, frontier, count, layer_degree, radius,
        polynomials, n_polynomials, field->context.big, n_threads,
        stats, other, this_meet, other_meet, cancel_requested,
        cancelled_out);
}

static void ball_enumerator_init(ball_enumerator_t *enumerator,
                                 search_graph_t *graph, jkey_t start,
                                 const char *label)
{
    fmpz_t one;
    size_t root;

    memset(enumerator, 0, sizeof(*enumerator));
    enumerator->graph = graph;
    enumerator->label = label;
    graph_init(graph);
    fmpz_init_set_ui(one, 1);
    root = graph_add(graph, start, one, NO_INDEX, 0);
    heap_insert(&enumerator->heap, root, graph);
    fmpz_clear(one);
}

static int ball_enumerator_has_work(const ball_enumerator_t *enumerator)
{
    return enumerator->heap.length != 0;
}

static const fmpz *ball_enumerator_next_degree(
    const ball_enumerator_t *enumerator)
{
    size_t index = heap_peek(&enumerator->heap);

    if (index == NO_INDEX)
        return NULL;
    return enumerator->graph->nodes[index].degree;
}

/*
 * Advance one complete multiplicative-degree layer. Collision tests happen
 * after every batch, so a successful search can return before either ball
 * has been fully enumerated. The returned parent chains are valid even for a
 * newly discovered (not yet settled) meeting vertex.
 */
static int ball_enumerator_expand_layer(
    ball_enumerator_t *enumerator, const search_graph_t *other,
    const fmpz_t radius, const modular_poly_t *polynomials,
    size_t n_polynomials, const field_context_t *field,
    unsigned n_threads,
    size_t *this_meet, size_t *other_meet,
    const atomic_int *cancel_requested, int *cancelled_out)
{
    search_graph_t *graph = enumerator->graph;
    fmpz_t layer_degree;
    size_t current_index, offset;

    if (!ball_enumerator_has_work(enumerator))
        return 0;
    if (cancel_requested != NULL &&
        atomic_load_explicit(cancel_requested, memory_order_relaxed)) {
        if (cancelled_out != NULL)
            *cancelled_out = 1;
        return 0;
    }

    fmpz_init(layer_degree);
    current_index = heap_pop(&enumerator->heap, graph);
    fmpz_set(layer_degree, graph->nodes[current_index].degree);
    enumerator->frontier_length = 0;

    for (;;) {
        if (enumerator->frontier_length == enumerator->frontier_alloc) {
            enumerator->frontier_alloc =
                enumerator->frontier_alloc
                    ? 2 * enumerator->frontier_alloc
                    : MULTIPOINT_BATCH_SIZE;
            enumerator->frontier =
                xrealloc(enumerator->frontier,
                         enumerator->frontier_alloc,
                         sizeof(*enumerator->frontier));
        }
        graph->nodes[current_index].settled = 1;
        enumerator->settled_count++;
        enumerator->frontier[enumerator->frontier_length++] =
            current_index;

        current_index = heap_peek(&enumerator->heap);
        if (current_index == NO_INDEX ||
            fmpz_cmp(graph->nodes[current_index].degree,
                     layer_degree) != 0)
            break;
        current_index = heap_pop(&enumerator->heap, graph);
    }

    for (offset = 0; offset < enumerator->frontier_length;
         offset += MULTIPOINT_BATCH_SIZE) {
        size_t count = enumerator->frontier_length - offset;

        if (cancel_requested != NULL &&
            atomic_load_explicit(cancel_requested,
                                 memory_order_relaxed)) {
            if (cancelled_out != NULL)
                *cancelled_out = 1;
            fmpz_clear(layer_degree);
            return 0;
        }
        if (count > MULTIPOINT_BATCH_SIZE)
            count = MULTIPOINT_BATCH_SIZE;
        if (process_frontier_batch(
                graph, &enumerator->heap,
                enumerator->frontier + offset, count, layer_degree,
                radius, polynomials, n_polynomials, field, n_threads,
                &enumerator->stats, other, this_meet, other_meet,
                cancel_requested, cancelled_out)) {
            fmpz_clear(layer_degree);
            return 1;
        }
    }

    fmpz_clear(layer_degree);
    return 0;
}

static void ball_enumerator_report(const ball_enumerator_t *enumerator,
                                   int stopped_early)
{
    fprintf(stderr,
            "%s ball: %lu j-invariants (%lu settled) in %.3f s%s "
            "(%llu specializations, %llu fast-eval batches, "
            "%llu Horner batches, %llu backtracks skipped)\n",
            enumerator->label,
            (unsigned long)enumerator->graph->length,
            (unsigned long)enumerator->settled_count,
            enumerator->active_seconds,
            stopped_early ? ", stopped at first intersection" : "",
            (unsigned long long)enumerator->stats.specializations,
            (unsigned long long)enumerator->stats.fast_batches,
            (unsigned long long)enumerator->stats.horner_batches,
            (unsigned long long)enumerator->stats.skipped_backtracks);
}

static void ball_enumerator_clear(ball_enumerator_t *enumerator)
{
    free(enumerator->frontier);
    heap_clear(&enumerator->heap);
    memset(enumerator, 0, sizeof(*enumerator));
}

/* ---------------------- intersection and output ----------------------- */

static int search_between(search_graph_t *left, search_graph_t *right,
                          jkey_t start, jkey_t target,
                          const fmpz_t degree_bound,
                          const modular_poly_t *polynomials,
                          size_t n_polynomials,
                          const field_context_t *field,
                          unsigned n_threads,
                          const char *left_label,
                          const char *right_label,
                          size_t *left_meet, size_t *right_meet,
                          fmpz_t path_degree, int report,
                          const atomic_int *cancel_requested,
                          int *cancelled_out)
{
    ball_enumerator_t left_enumerator, right_enumerator;
    fmpz_t radius;
    int found = 0, expand_left, alternate_tie = 0;

    fmpz_init(radius);
    fmpz_sqrt(radius, degree_bound);
    if (cancelled_out != NULL)
        *cancelled_out = 0;
    if (report) {
        fprintf(stderr, "Search radius floor(sqrt(bound)) = ");
        fmpz_fprint(stderr, radius);
        fprintf(stderr, "\n");
    }

    ball_enumerator_init(&left_enumerator, left, start, left_label);
    ball_enumerator_init(&right_enumerator, right, target, right_label);

    if (key_equal(start, target)) {
        *left_meet = 0;
        *right_meet = 0;
        fmpz_one(path_degree);
        found = 1;
    }

    while (!found &&
           (ball_enumerator_has_work(&left_enumerator) ||
            ball_enumerator_has_work(&right_enumerator))) {
        double layer_started;

        if (cancel_requested != NULL &&
            atomic_load_explicit(cancel_requested,
                                 memory_order_relaxed)) {
            if (cancelled_out != NULL)
                *cancelled_out = 1;
            break;
        }
        if (!ball_enumerator_has_work(&right_enumerator)) {
            expand_left = 1;
        } else if (!ball_enumerator_has_work(&left_enumerator)) {
            expand_left = 0;
        } else {
            int comparison = fmpz_cmp(
                ball_enumerator_next_degree(&left_enumerator),
                ball_enumerator_next_degree(&right_enumerator));

            if (comparison < 0)
                expand_left = 1;
            else if (comparison > 0)
                expand_left = 0;
            else {
                expand_left = !alternate_tie;
                alternate_tie = !alternate_tie;
            }
        }

        layer_started = wall_seconds();
        if (expand_left) {
            found = ball_enumerator_expand_layer(
                &left_enumerator, right, radius, polynomials,
                n_polynomials, field, n_threads,
                left_meet, right_meet, cancel_requested,
                cancelled_out);
            left_enumerator.active_seconds +=
                wall_seconds() - layer_started;
        } else {
            found = ball_enumerator_expand_layer(
                &right_enumerator, left, radius, polynomials,
                n_polynomials, field, n_threads,
                right_meet, left_meet, cancel_requested,
                cancelled_out);
            right_enumerator.active_seconds +=
                wall_seconds() - layer_started;
        }
        if (cancelled_out != NULL && *cancelled_out)
            break;
    }

    if (found) {
        fmpz_mul(path_degree, left->nodes[*left_meet].degree,
                 right->nodes[*right_meet].degree);
        if (fmpz_cmp(path_degree, degree_bound) > 0)
            die("internal error: meeting path exceeds degree bound");
    }

    if (report) {
        ball_enumerator_report(&left_enumerator, found);
        ball_enumerator_report(&right_enumerator, found);
    }
    ball_enumerator_clear(&left_enumerator);
    ball_enumerator_clear(&right_enumerator);
    fmpz_clear(radius);
    return found;
}
#endif

/* ------------------- canonical prime-order search -------------------- */

static void index_vector_append(index_vector_t *vector, size_t value)
{
    if (vector->length == vector->alloc) {
        vector->alloc = vector->alloc ? 2 * vector->alloc : 1024;
        vector->items = xrealloc(vector->items, vector->alloc,
                                 sizeof(*vector->items));
    }
    vector->items[vector->length++] = value;
}

static void index_vector_clear(index_vector_t *vector)
{
    free(vector->items);
    memset(vector, 0, sizeof(*vector));
}

static int compare_size_t_values(const void *left, const void *right)
{
    size_t a = *(const size_t *)left;
    size_t b = *(const size_t *)right;

    return a < b ? -1 : a > b;
}

static void index_vector_sort_unique(index_vector_t *vector)
{
    size_t source, target;

    if (vector->length < 2)
        return;
    qsort(vector->items, vector->length, sizeof(*vector->items),
          compare_size_t_values);
    target = 1;
    for (source = 1; source < vector->length; source++)
        if (vector->items[source] != vector->items[target - 1])
            vector->items[target++] = vector->items[source];
    vector->length = target;
}

static void index_vector_fill_graph(index_vector_t *vector,
                                    const search_graph_t *graph)
{
    size_t i;

    vector->length = 0;
    if (vector->alloc < graph->length) {
        vector->alloc = graph->length;
        vector->items = xrealloc(vector->items, vector->alloc,
                                 sizeof(*vector->items));
    }
    for (i = 0; i < graph->length; i++)
        vector->items[vector->length++] = i;
}

static void index_vector_filter_eligible(index_vector_t *vector,
                                         const search_graph_t *graph,
                                         unsigned ell,
                                         const fmpz_t radius)
{
    fmpz_t next_degree;
    size_t source, target = 0;

    fmpz_init(next_degree);
    for (source = 0; source < vector->length; source++) {
        size_t index = vector->items[source];

        fmpz_mul_ui(next_degree, graph->nodes[index].degree, ell);
        if (fmpz_cmp(next_degree, radius) <= 0)
            vector->items[target++] = index;
    }
    vector->length = target;
    fmpz_clear(next_degree);
}

/*
 * Insert a newly reached endpoint, or retain the smaller canonical degree
 * when the endpoint was already known. Returning nonzero asks the current
 * ell-stage to propagate this new or improved state through another
 * ell-isogeny layer.
 */
static int ordered_relax_node(search_graph_t *graph, jkey_t key,
                              const fmpz_t degree, size_t parent,
                              unsigned ell, size_t *index_out)
{
    size_t index = graph_find(graph, key);

    if (index == NO_INDEX) {
        index = graph_add(graph, key, degree, parent, ell);
        *index_out = index;
        return 1;
    }
    *index_out = index;
    if (fmpz_cmp(degree, graph->nodes[index].degree) < 0) {
        fmpz_set(graph->nodes[index].degree, degree);
        graph->nodes[index].parent = parent;
        graph->nodes[index].parent_ell = ell;
        return 1;
    }
    return 0;
}

static void graph_chain_degree(fmpz_t degree, const search_graph_t *graph,
                               size_t endpoint)
{
    size_t current = endpoint, traversed = 0;

    fmpz_one(degree);
    while (graph->nodes[current].parent != NO_INDEX) {
        fmpz_mul_ui(degree, degree, graph->nodes[current].parent_ell);
        current = graph->nodes[current].parent;
        if (++traversed > graph->length)
            die("internal error: cycle in stored parent chain");
    }
}

static int process_ordered_batch_word(
    ordered_ball_t *ball, const size_t *frontier, size_t count,
    const fmpz_t radius, const modular_poly_t *phi,
    const fq_nmod_ctx_t ctx, unsigned n_threads,
    const search_graph_t *other, index_vector_t *next_frontier,
    size_t *this_meet, size_t *other_meet,
    const atomic_int *cancel_requested, int *cancelled_out)
{
    search_graph_t *graph = ball->graph;
    fq_nmod_struct *points, *coefficient_values;
    fq_nmod_poly_struct *evaluated;
    fq_nmod_poly_struct **tree = NULL;
    jkey_t *root_keys;
    size_t *root_counts;
    nmod_poly_t representative;
    fmpz_t next_degree;
    size_t i, root_task_count;
    size_t coefficient_count = phi->n_coefficients;
    size_t root_stride = (size_t)phi->ell + 1;
    int use_fast, found = 0;
    word_evaluation_tasks_t evaluation_tasks;
    word_specialization_tasks_t specialization_tasks;
    word_root_tasks_t root_tasks;

    if (count == 0)
        return 0;
    if (count > (size_t)LONG_MAX)
        die("multipoint batch is too large for FLINT's slong interface");
    if (coefficient_count != 0 &&
        count > SIZE_MAX / coefficient_count)
        die("multipoint coefficient buffer size overflow");
    if (root_stride != 0 && count > SIZE_MAX / root_stride)
        die("multipoint root buffer size overflow");
    if (cancel_requested != NULL &&
        atomic_load_explicit(cancel_requested, memory_order_relaxed)) {
        if (cancelled_out != NULL)
            *cancelled_out = 1;
        return 0;
    }

    points = xrealloc(NULL, count, sizeof(*points));
    coefficient_values =
        xrealloc(NULL, count * coefficient_count,
                 sizeof(*coefficient_values));
    evaluated = xrealloc(NULL, count, sizeof(*evaluated));
    root_keys = xrealloc(NULL, count * root_stride, sizeof(*root_keys));
    root_counts = xrealloc(NULL, count, sizeof(*root_counts));
    nmod_poly_init(representative, field_characteristic_word);
    key_array_init(root_keys, count * root_stride);
    for (i = 0; i < count; i++) {
        fq_nmod_init(points + i, ctx);
        fq_nmod_poly_init(evaluated + i, ctx);
        key_to_fq_nmod(points + i, graph->nodes[frontier[i]].key,
                       representative, ctx);
    }
    for (i = 0; i < count * coefficient_count; i++)
        fq_nmod_init(coefficient_values + i, ctx);

    use_fast = ordered_force_multipoint ||
               (count >= MULTIPOINT_THRESHOLD &&
                phi->max_x_degree >= MULTIPOINT_MIN_DEGREE);
    if (use_fast) {
        tree = _fq_nmod_poly_tree_alloc((slong)count, ctx);
        _fq_nmod_poly_tree_build(tree, points, (slong)count, ctx);
    }

    evaluation_tasks.phi = phi;
    evaluation_tasks.points = points;
    evaluation_tasks.count = count;
    evaluation_tasks.coefficient_values = coefficient_values;
    evaluation_tasks.tree = tree;
    evaluation_tasks.use_fast = use_fast;
    evaluation_tasks.ctx = ctx;
    flint_parallel_do(evaluate_coefficient_task_word, &evaluation_tasks,
                      (slong)coefficient_count, (int)n_threads,
                      FLINT_PARALLEL_DYNAMIC);
    if (use_fast)
        ball->stats.fast_batches++;
    else
        ball->stats.horner_batches++;

    specialization_tasks.phi = phi;
    specialization_tasks.count = count;
    specialization_tasks.coefficient_values = coefficient_values;
    specialization_tasks.evaluated = evaluated;
    specialization_tasks.ctx = ctx;
    flint_parallel_do(build_specialization_task_word,
                      &specialization_tasks, (slong)count,
                      (int)n_threads, FLINT_PARALLEL_UNIFORM);

    ball->stats.specializations += (uint64_t)count;
    ball->stats.root_finds += (uint64_t)count;
    ball->expanded_states += count;

    root_task_count = count;
    if (root_task_count > (size_t)n_threads * 4)
        root_task_count = (size_t)n_threads * 4;
    root_tasks.evaluated = evaluated;
    root_tasks.count = count;
    root_tasks.task_count = root_task_count;
    root_tasks.root_stride = root_stride;
    root_tasks.root_keys = root_keys;
    root_tasks.root_counts = root_counts;
    root_tasks.ctx = ctx;
    flint_parallel_do(find_roots_task_word, &root_tasks,
                      (slong)root_task_count, (int)n_threads,
                      FLINT_PARALLEL_DYNAMIC);

    fmpz_init(next_degree);
    for (i = 0; i < count && !found; i++) {
        size_t root_index;

        fmpz_mul_ui(next_degree, graph->nodes[frontier[i]].degree,
                    phi->ell);
        if (fmpz_cmp(next_degree, radius) > 0)
            continue;
        for (root_index = 0; root_index < root_counts[i]; root_index++) {
            jkey_t root_key = root_keys[i * root_stride + root_index];
            const search_node_t *current =
                graph->nodes + frontier[i];
            size_t index, match;

            if (current->parent != NO_INDEX &&
                current->parent_ell == phi->ell &&
                key_equal(root_key,
                          graph->nodes[current->parent].key)) {
                ball->stats.skipped_backtracks++;
                continue;
            }
            if (!ordered_relax_node(graph, root_key, next_degree,
                                    frontier[i], phi->ell, &index))
                continue;
            index_vector_append(next_frontier, index);
            match = graph_find(other, graph->nodes[index].key);
            if (match != NO_INDEX) {
                *this_meet = index;
                *other_meet = match;
                found = 1;
                break;
            }
        }
    }
    fmpz_clear(next_degree);

    if (tree != NULL)
        _fq_nmod_poly_tree_free(tree, (slong)count, ctx);
    for (i = 0; i < count * coefficient_count; i++)
        fq_nmod_clear(coefficient_values + i, ctx);
    for (i = 0; i < count; i++) {
        fq_nmod_poly_clear(evaluated + i, ctx);
        fq_nmod_clear(points + i, ctx);
    }
    nmod_poly_clear(representative);
    key_array_clear(root_keys, count * root_stride);
    free(root_counts);
    free(root_keys);
    free(evaluated);
    free(coefficient_values);
    free(points);
    return found;
}

static int process_ordered_batch_big(
    ordered_ball_t *ball, const size_t *frontier, size_t count,
    const fmpz_t radius, const modular_poly_t *phi,
    const fq_ctx_t ctx, unsigned n_threads,
    const search_graph_t *other, index_vector_t *next_frontier,
    size_t *this_meet, size_t *other_meet,
    const atomic_int *cancel_requested, int *cancelled_out)
{
    search_graph_t *graph = ball->graph;
    fq_struct *points, *coefficient_values;
    fq_poly_struct *evaluated;
    fq_poly_struct **tree = NULL;
    jkey_t *root_keys;
    size_t *root_counts;
    fmpz_poly_t representative;
    fmpz_t next_degree;
    size_t i, root_task_count;
    size_t coefficient_count = phi->n_coefficients;
    size_t root_stride = (size_t)phi->ell + 1;
    int use_fast, found = 0;
    big_evaluation_tasks_t evaluation_tasks;
    big_specialization_tasks_t specialization_tasks;
    big_root_tasks_t root_tasks;

    if (count == 0)
        return 0;
    if (count > (size_t)LONG_MAX)
        die("multipoint batch is too large for FLINT's slong interface");
    if (coefficient_count != 0 &&
        count > SIZE_MAX / coefficient_count)
        die("multipoint coefficient buffer size overflow");
    if (root_stride != 0 && count > SIZE_MAX / root_stride)
        die("multipoint root buffer size overflow");
    if (cancel_requested != NULL &&
        atomic_load_explicit(cancel_requested, memory_order_relaxed)) {
        if (cancelled_out != NULL)
            *cancelled_out = 1;
        return 0;
    }

    points = xrealloc(NULL, count, sizeof(*points));
    coefficient_values =
        xrealloc(NULL, count * coefficient_count,
                 sizeof(*coefficient_values));
    evaluated = xrealloc(NULL, count, sizeof(*evaluated));
    root_keys = xrealloc(NULL, count * root_stride, sizeof(*root_keys));
    root_counts = xrealloc(NULL, count, sizeof(*root_counts));
    fmpz_poly_init(representative);
    key_array_init(root_keys, count * root_stride);
    for (i = 0; i < count; i++) {
        fq_init(points + i, ctx);
        fq_poly_init(evaluated + i, ctx);
        key_to_fq(points + i, graph->nodes[frontier[i]].key,
                  representative, ctx);
    }
    for (i = 0; i < count * coefficient_count; i++)
        fq_init(coefficient_values + i, ctx);

    use_fast = ordered_force_multipoint ||
               (count >= MULTIPOINT_THRESHOLD &&
                phi->max_x_degree >= MULTIPOINT_MIN_DEGREE);
    if (use_fast) {
        tree = _fq_poly_tree_alloc((slong)count, ctx);
        _fq_poly_tree_build(tree, points, (slong)count, ctx);
    }

    evaluation_tasks.phi = phi;
    evaluation_tasks.points = points;
    evaluation_tasks.count = count;
    evaluation_tasks.coefficient_values = coefficient_values;
    evaluation_tasks.tree = tree;
    evaluation_tasks.use_fast = use_fast;
    evaluation_tasks.ctx = ctx;
    flint_parallel_do(evaluate_coefficient_task_big, &evaluation_tasks,
                      (slong)coefficient_count, (int)n_threads,
                      FLINT_PARALLEL_DYNAMIC);
    if (use_fast)
        ball->stats.fast_batches++;
    else
        ball->stats.horner_batches++;

    specialization_tasks.phi = phi;
    specialization_tasks.count = count;
    specialization_tasks.coefficient_values = coefficient_values;
    specialization_tasks.evaluated = evaluated;
    specialization_tasks.ctx = ctx;
    flint_parallel_do(build_specialization_task_big,
                      &specialization_tasks, (slong)count,
                      (int)n_threads, FLINT_PARALLEL_UNIFORM);

    ball->stats.specializations += (uint64_t)count;
    ball->stats.root_finds += (uint64_t)count;
    ball->expanded_states += count;

    root_task_count = count;
    if (root_task_count > (size_t)n_threads * 4)
        root_task_count = (size_t)n_threads * 4;
    root_tasks.evaluated = evaluated;
    root_tasks.count = count;
    root_tasks.task_count = root_task_count;
    root_tasks.root_stride = root_stride;
    root_tasks.root_keys = root_keys;
    root_tasks.root_counts = root_counts;
    root_tasks.ctx = ctx;
    flint_parallel_do(find_roots_task_big, &root_tasks,
                      (slong)root_task_count, (int)n_threads,
                      FLINT_PARALLEL_DYNAMIC);

    fmpz_init(next_degree);
    for (i = 0; i < count && !found; i++) {
        size_t root_index;

        fmpz_mul_ui(next_degree, graph->nodes[frontier[i]].degree,
                    phi->ell);
        if (fmpz_cmp(next_degree, radius) > 0)
            continue;
        for (root_index = 0; root_index < root_counts[i]; root_index++) {
            jkey_t root_key = root_keys[i * root_stride + root_index];
            const search_node_t *current =
                graph->nodes + frontier[i];
            size_t index, match;

            if (current->parent != NO_INDEX &&
                current->parent_ell == phi->ell &&
                key_equal(root_key,
                          graph->nodes[current->parent].key)) {
                ball->stats.skipped_backtracks++;
                continue;
            }
            if (!ordered_relax_node(graph, root_key, next_degree,
                                    frontier[i], phi->ell, &index))
                continue;
            index_vector_append(next_frontier, index);
            match = graph_find(other, graph->nodes[index].key);
            if (match != NO_INDEX) {
                *this_meet = index;
                *other_meet = match;
                found = 1;
                break;
            }
        }
    }
    fmpz_clear(next_degree);

    if (tree != NULL)
        _fq_poly_tree_free(tree, (slong)count, ctx);
    for (i = 0; i < count * coefficient_count; i++)
        fq_clear(coefficient_values + i, ctx);
    for (i = 0; i < count; i++) {
        fq_poly_clear(evaluated + i, ctx);
        fq_clear(points + i, ctx);
    }
    fmpz_poly_clear(representative);
    key_array_clear(root_keys, count * root_stride);
    free(root_counts);
    free(root_keys);
    free(evaluated);
    free(coefficient_values);
    free(points);
    return found;
}

static int process_ordered_batch(
    ordered_ball_t *ball, const size_t *frontier, size_t count,
    const fmpz_t radius, const modular_poly_t *phi,
    const field_context_t *field, unsigned n_threads,
    const search_graph_t *other, index_vector_t *next_frontier,
    size_t *this_meet, size_t *other_meet,
    const atomic_int *cancel_requested, int *cancelled_out)
{
    if (field->kind == BACKEND_FQ_NMOD)
        return process_ordered_batch_word(
            ball, frontier, count, radius, phi, field->context.word,
            n_threads, other, next_frontier, this_meet, other_meet,
            cancel_requested, cancelled_out);
    return process_ordered_batch_big(
        ball, frontier, count, radius, phi, field->context.big,
        n_threads, other, next_frontier, this_meet, other_meet,
        cancel_requested, cancelled_out);
}

static size_t ordered_batch_count(size_t remaining)
{
    if (ordered_multipoint_batch_limit != 0 &&
        remaining > ordered_multipoint_batch_limit)
        return ordered_multipoint_batch_limit;
    return remaining;
}

/*
 * Generate all canonical ell-power extensions of the endpoints produced by
 * smaller primes. Both balls advance one batch at a time. The k-th wave
 * depends on the (k-1)-st, but sources of different current degrees in one
 * wave share a single Phi_ell multipoint tree.
 */
static int ordered_expand_prime_stage(
    ordered_ball_t *left_ball, ordered_ball_t *right_ball,
    const fmpz_t radius, const modular_poly_t *phi,
    const field_context_t *field, unsigned n_threads,
    size_t *left_meet, size_t *right_meet,
    const atomic_int *cancel_requested, int *cancelled_out)
{
    index_vector_t left_current = {0}, right_current = {0};
    index_vector_t left_next = {0}, right_next = {0};
    int found = 0, prefer_left = 1;

    index_vector_fill_graph(&left_current, left_ball->graph);
    index_vector_fill_graph(&right_current, right_ball->graph);

    while (left_current.length != 0 || right_current.length != 0) {
        size_t left_offset = 0, right_offset = 0;

        index_vector_filter_eligible(
            &left_current, left_ball->graph, phi->ell, radius);
        index_vector_filter_eligible(
            &right_current, right_ball->graph, phi->ell, radius);
        if (left_current.length == 0 && right_current.length == 0)
            break;
        left_next.length = 0;
        right_next.length = 0;

        while (!found &&
               (left_offset < left_current.length ||
                right_offset < right_current.length)) {
            int expand_left;
            size_t count;
            double started;

            if (cancel_requested != NULL &&
                atomic_load_explicit(cancel_requested,
                                     memory_order_relaxed)) {
                if (cancelled_out != NULL)
                    *cancelled_out = 1;
                goto cleanup;
            }
            if (right_offset >= right_current.length)
                expand_left = 1;
            else if (left_offset >= left_current.length)
                expand_left = 0;
            else {
                expand_left = prefer_left;
                prefer_left = !prefer_left;
            }

            started = wall_seconds();
            if (expand_left) {
                count = ordered_batch_count(
                    left_current.length - left_offset);
                found = process_ordered_batch(
                    left_ball, left_current.items + left_offset, count,
                    radius, phi, field, n_threads, right_ball->graph,
                    &left_next, left_meet, right_meet,
                    cancel_requested, cancelled_out);
                left_offset += count;
                left_ball->active_seconds += wall_seconds() - started;
            } else {
                count = ordered_batch_count(
                    right_current.length - right_offset);
                found = process_ordered_batch(
                    right_ball, right_current.items + right_offset,
                    count, radius, phi, field, n_threads,
                    left_ball->graph, &right_next,
                    right_meet, left_meet,
                    cancel_requested, cancelled_out);
                right_offset += count;
                right_ball->active_seconds += wall_seconds() - started;
            }
            if (cancelled_out != NULL && *cancelled_out)
                goto cleanup;
        }
        if (found)
            break;

        index_vector_sort_unique(&left_next);
        index_vector_sort_unique(&right_next);
        {
            index_vector_t temporary = left_current;

            left_current = left_next;
            left_next = temporary;
            left_next.length = 0;
        }
        {
            index_vector_t temporary = right_current;

            right_current = right_next;
            right_next = temporary;
            right_next.length = 0;
        }
    }

cleanup:
    index_vector_clear(&right_next);
    index_vector_clear(&left_next);
    index_vector_clear(&right_current);
    index_vector_clear(&left_current);
    return found;
}

static void ordered_ball_report(const ordered_ball_t *ball,
                                int stopped_early)
{
    fprintf(stderr,
            "%s ball: %lu j-invariants, %lu source expansions in %.3f s%s "
            "(%llu specializations, %llu fast-eval batches, "
            "%llu Horner batches, %llu backtracks skipped)\n",
            ball->label,
            (unsigned long)ball->graph->length,
            (unsigned long)ball->expanded_states,
            ball->active_seconds,
            stopped_early ? ", stopped at first intersection" : "",
            (unsigned long long)ball->stats.specializations,
            (unsigned long long)ball->stats.fast_batches,
            (unsigned long long)ball->stats.horner_batches,
            (unsigned long long)ball->stats.skipped_backtracks);
}

static int ordered_search_between(
    search_graph_t *left, search_graph_t *right,
    jkey_t start, jkey_t target, const fmpz_t degree_bound,
    const modular_poly_t *polynomials, size_t n_polynomials,
    const field_context_t *field, unsigned n_threads,
    const char *left_label, const char *right_label,
    size_t *left_meet, size_t *right_meet,
    fmpz_t path_degree, int report,
    const atomic_int *cancel_requested, int *cancelled_out)
{
    ordered_ball_t left_ball, right_ball;
    fmpz_t radius, left_degree, right_degree;
    size_t polynomial_index;
    int found = 0;

    memset(&left_ball, 0, sizeof(left_ball));
    memset(&right_ball, 0, sizeof(right_ball));
    left_ball.graph = left;
    right_ball.graph = right;
    left_ball.label = left_label;
    right_ball.label = right_label;
    if (cancelled_out != NULL)
        *cancelled_out = 0;

    fmpz_init(radius);
    fmpz_init(left_degree);
    fmpz_init(right_degree);
    fmpz_sqrt(radius, degree_bound);
    graph_init(left);
    graph_init(right);
    fmpz_one(left_degree);
    graph_add(left, start, left_degree, NO_INDEX, 0);
    graph_add(right, target, left_degree, NO_INDEX, 0);

    if (report) {
        fprintf(stderr, "Search radius floor(sqrt(bound)) = ");
        fmpz_fprint(stderr, radius);
        fprintf(stderr, "\nEnumeration order = canonical prime blocks");
        if (ordered_multipoint_batch_limit == 0)
            fprintf(stderr,
                    "; multipoint batch cap = unlimited "
                    "(one complete eligible layer)\n");
        else
            fprintf(stderr, "; multipoint batch cap = %lu\n",
                    (unsigned long)ordered_multipoint_batch_limit);
        fprintf(stderr, "Ordered-search evaluation policy = %s\n",
                ordered_force_multipoint
                    ? "multipoint only (forced for every nonempty batch)"
                    : "adaptive Horner/multipoint crossover");
    }

    if (key_equal(start, target)) {
        *left_meet = 0;
        *right_meet = 0;
        fmpz_one(path_degree);
        found = 1;
    }

    for (polynomial_index = 0;
         polynomial_index < n_polynomials && !found;
         polynomial_index++) {
        if (report)
            fprintf(stderr, "Prime stage ell = %u\n",
                    polynomials[polynomial_index].ell);
        found = ordered_expand_prime_stage(
            &left_ball, &right_ball, radius,
            polynomials + polynomial_index, field, n_threads,
            left_meet, right_meet, cancel_requested, cancelled_out);
        if (cancelled_out != NULL && *cancelled_out)
            break;
    }

    if (found) {
        graph_chain_degree(left_degree, left, *left_meet);
        graph_chain_degree(right_degree, right, *right_meet);
        if (fmpz_cmp(left_degree, radius) > 0 ||
            fmpz_cmp(right_degree, radius) > 0)
            die("internal error: meeting endpoint exceeds search radius");
        fmpz_mul(path_degree, left_degree, right_degree);
        if (fmpz_cmp(path_degree, degree_bound) > 0)
            die("internal error: meeting path exceeds degree bound");
    }

    if (report) {
        ordered_ball_report(&left_ball, found);
        ordered_ball_report(&right_ball, found);
    }
    fmpz_clear(right_degree);
    fmpz_clear(left_degree);
    fmpz_clear(radius);
    return found;
}

typedef struct {
    search_graph_t left, right;
    jkey_t start, target;
    size_t endpoint;
    size_t attempt_number;
    size_t layer_position;
    size_t left_meet, right_meet;
    fmpz_t path_degree;
    const fmpz *degree_bound;
    const modular_poly_t *polynomials;
    size_t n_polynomials;
    const field_context_t *field;
    unsigned search_threads;
    const char *right_label;
    atomic_int *cancel_requested;
    atomic_int *winner_slot;
    int slot;
    int found;
    int cancelled;
    double elapsed_seconds;
} rerandomization_search_job_t;

static void *run_rerandomization_search(void *argument)
{
    rerandomization_search_job_t *job = argument;
    double started = wall_seconds();

    job->found = ordered_search_between(
        &job->left, &job->right, job->start, job->target,
        job->degree_bound, job->polynomials, job->n_polynomials,
        job->field, job->search_threads, "j'", job->right_label,
        &job->left_meet, &job->right_meet, job->path_degree, 0,
        job->cancel_requested, &job->cancelled);
    job->elapsed_seconds = wall_seconds() - started;

    if (job->found) {
        int expected = -1;

        if (atomic_compare_exchange_strong_explicit(
                job->winner_slot, &expected, job->slot,
                memory_order_relaxed, memory_order_relaxed))
            atomic_store_explicit(job->cancel_requested, 1,
                                  memory_order_relaxed);
    }
    return NULL;
}

static void print_path(const search_graph_t *left, size_t left_meet,
                       const search_graph_t *right, size_t right_meet,
                       const fmpz_t total_degree)
{
    size_t *reverse_path = NULL;
    size_t length = 0, right_length = 0, alloc = 0, i, current;

    current = left_meet;
    while (current != NO_INDEX) {
        if (length == alloc) {
            alloc = alloc ? 2 * alloc : 16;
            reverse_path = xrealloc(reverse_path, alloc,
                                    sizeof(*reverse_path));
        }
        reverse_path[length++] = current;
        current = left->nodes[current].parent;
    }

    current = right_meet;
    while (right->nodes[current].parent != NO_INDEX) {
        right_length++;
        current = right->nodes[current].parent;
    }

    printf("Path found. Total degree = ");
    fmpz_print(total_degree);
    printf("\nNumber of steps = %lu\n\n",
           (unsigned long)(length - 1 + right_length));

    print_key(left->nodes[reverse_path[length - 1]].key);
    printf("\n");
    for (i = length - 1; i > 0; i--) {
        const search_node_t *child = left->nodes + reverse_path[i - 1];

        printf("  --[%u]--> ", child->parent_ell);
        print_key(child->key);
        printf("\n");
    }

    current = right_meet;
    while (right->nodes[current].parent != NO_INDEX) {
        printf("  --[%u]--> ", right->nodes[current].parent_ell);
        current = right->nodes[current].parent;
        print_key(right->nodes[current].key);
        printf("\n");
    }
    free(reverse_path);
}

static void print_rerandomized_path(const search_graph_t *left,
                                    size_t left_meet,
                                    const search_graph_t *right,
                                    size_t right_meet,
                                    const jkey_t *walk,
                                    size_t walk_length,
                                    unsigned rerandomization_ell,
                                    const fmpz_t middle_degree,
                                    int append_frobenius_return)
{
    size_t *reverse_path = NULL;
    size_t length = 0, right_length = 0, alloc = 0, i, current;
    size_t rerandomization_depth = walk_length - 1;
    fmpz_t total_degree;

    current = left_meet;
    while (current != NO_INDEX) {
        if (length == alloc) {
            alloc = alloc ? 2 * alloc : 16;
            reverse_path = xrealloc(reverse_path, alloc,
                                    sizeof(*reverse_path));
        }
        reverse_path[length++] = current;
        current = left->nodes[current].parent;
    }

    current = right_meet;
    while (right->nodes[current].parent != NO_INDEX) {
        right_length++;
        current = right->nodes[current].parent;
    }

    fmpz_init(total_degree);
    fmpz_set(total_degree, middle_degree);
    for (i = 0;
         i < rerandomization_depth *
                 (append_frobenius_return ? 2 : 1);
         i++)
        fmpz_mul_ui(total_degree, total_degree, rerandomization_ell);

    printf("Path found after rerandomization. Total degree = ");
    fmpz_print(total_degree);
    printf("\nSmooth middle degree = ");
    fmpz_print(middle_degree);
    printf(" (bound N)\nNumber of steps = %lu\n",
           (unsigned long)(length - 1 + right_length +
                           rerandomization_depth *
                               (append_frobenius_return ? 2 : 1)));
    printf("Rerandomization prime ell' = %u\n", rerandomization_ell);
    printf("Rerandomization depth = %lu\n\n",
           (unsigned long)rerandomization_depth);

    print_key(walk[0]);
    printf("\n");
    for (i = 1; i < walk_length; i++) {
        printf("  --[%u]--> ", rerandomization_ell);
        print_key(walk[i]);
        printf("\n");
    }
    for (i = length - 1; i > 0; i--) {
        const search_node_t *child = left->nodes + reverse_path[i - 1];

        printf("  --[%u]--> ", child->parent_ell);
        print_key(child->key);
        printf("\n");
    }

    current = right_meet;
    while (right->nodes[current].parent != NO_INDEX) {
        printf("  --[%u]--> ", right->nodes[current].parent_ell);
        current = right->nodes[current].parent;
        print_key(right->nodes[current].key);
        printf("\n");
    }
    if (append_frobenius_return) {
        jkey_t conjugate;

        key_init(&conjugate);
        for (i = walk_length - 1; i > 0; i--) {
            printf("  --[%u]--> ", rerandomization_ell);
            frobenius_conjugate(&conjugate, walk + i - 1);
            print_key(conjugate);
            printf("\n");
        }
        key_clear(&conjugate);
    }

    fmpz_clear(total_degree);
    free(reverse_path);
}

/* ------------------------------- main --------------------------------- */

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

static void parse_field_coefficient(fmpz_t value, const char *s)
{
    if (fmpz_set_str(value, s, 10) != 0 || fmpz_sgn(value) < 0 ||
        fmpz_cmp(value, field_characteristic) >= 0)
        die("field coefficients must be decimal integers in [0,p)");
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

/* Return floor(N^(1/4)), subject to the program's supported B range. */
static unsigned default_bound_from_N(const fmpz_t N)
{
    fmpz_t fourth_root;
    unsigned B;

    fmpz_init(fourth_root);
    (void)fmpz_root(fourth_root, N, 4);
    if (fmpz_cmp_ui(fourth_root, 2) < 0 ||
        fmpz_cmp_ui(fourth_root, MAX_SMOOTHNESS_BOUND) > 0) {
        fmpz_clear(fourth_root);
        die("default floor(N^(1/4)) must satisfy 2 <= B <= 10000; "
            "supply --B explicitly");
    }
    B = (unsigned)fmpz_get_ui(fourth_root);
    fmpz_clear(fourth_root);
    return B;
}

static void choose_quadratic_nonresidue(fmpz_t q, const fmpz_t p)
{
    /* This choice keeps the historical a^2 + 1 basis whenever possible. */
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

static unsigned parse_bound(const char *s)
{
    char *end;
    unsigned long value;

    errno = 0;
    value = strtoul(s, &end, 10);
    if (errno || *end != '\0' || value < 2 ||
        value > MAX_SMOOTHNESS_BOUND)
        die("B must satisfy 2 <= B <= 10000");
    return (unsigned)value;
}

static unsigned parse_thread_count(const char *s)
{
    char *end;
    unsigned long value;

    errno = 0;
    value = strtoul(s, &end, 10);
    if (errno || *end != '\0' || value < 1 || value > 256)
        die("thread count must satisfy 1 <= n <= 256");
    return (unsigned)value;
}

static size_t parse_multipoint_batch_limit(const char *s)
{
    char *end;
    unsigned long long value;

    errno = 0;
    value = strtoull(s, &end, 10);
    if (s[0] == '-' || errno || *end != '\0' ||
        value > (unsigned long long)SIZE_MAX)
        die("multipoint batch limit must be a nonnegative decimal "
            "integer");
    return (size_t)value;
}

static unsigned parse_rerandomization_job_count(const char *s)
{
    char *end;
    unsigned long value;

    errno = 0;
    value = strtoul(s, &end, 10);
    if (errno || *end != '\0' || value < 1 ||
        value > MAX_RERANDOMIZATION_JOBS)
        die("rerandomization job count must satisfy 1 <= n <= 256");
    return (unsigned)value;
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

static int parse_on_off(const char *s)
{
    if (strcmp(s, "on") == 0)
        return 1;
    if (strcmp(s, "off") == 0)
        return 0;
    die("rerandomization must be either 'on' or 'off'");
    return 0;
}

static int is_prime(unsigned n)
{
    unsigned d;

    if (n < 2)
        return 0;
    for (d = 2; d * d <= n; d++)
        if (n % d == 0)
            return 0;
    return 1;
}

static unsigned smallest_prime_larger_than(unsigned n)
{
    unsigned candidate = n + 1;

    while (!is_prime(candidate))
        candidate++;
    return candidate;
}

static void frobenius_conjugate(jkey_t *result, const jkey_t *j)
{
    key_set(result, j);
    if (active_backend == BACKEND_FQ_NMOD) {
        if (result->value.word.im != 0)
            result->value.word.im =
                field_characteristic_word - result->value.word.im;
    } else if (!fmpz_is_zero(&result->value.big.im)) {
        fmpz_sub(&result->value.big.im, field_characteristic,
                 &result->value.big.im);
    }
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
                die("clock_gettime failed while seeding random generation");
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

static size_t ceil_log2_characteristic(void)
{
    return (size_t)fmpz_bits(field_characteristic);
}

static void generate_random_start_j(
    jkey_t *result, const modular_poly_t *phi_2,
    const field_context_t *field,
    random_state_t *random_state, size_t walk_length)
{
    jkey_t previous, current;
    int have_previous = 0;
    size_t step;

    if (phi_2->ell != 2)
        die("the first smooth modular polynomial is not Phi_2");
    key_init(&previous);
    key_init(&current);
    if (active_backend == BACKEND_FQ_NMOD) {
        key_set_ui(&current, UWORD(1728) % field_characteristic_word, 0);
    } else {
        fmpz_set_ui(&current.value.big.re, UWORD(1728));
        fmpz_mod(&current.value.big.re, &current.value.big.re,
                 field_characteristic);
        fmpz_zero(&current.value.big.im);
    }

    for (step = 0; step < walk_length; step++) {
        jkey_t *edges = NULL;
        size_t edge_count, available, edge_index, forbidden = NO_INDEX;
        size_t choice, selected = NO_INDEX;

        edge_count =
            modular_poly_edges_at(&edges, current, phi_2, field);
        if (edge_count != 3)
            die("Phi_2 specialization did not split into three "
                "F_{p^2}-rational edges during j-generation");
        available = edge_count;

        if (have_previous) {
            for (edge_index = 0; edge_index < edge_count; edge_index++) {
                if (key_equal(edges[edge_index], previous)) {
                    forbidden = edge_index;
                    break;
                }
            }
            if (forbidden != NO_INDEX)
                available--;
            if (available != 2)
                die("could not identify the dual 2-isogeny edge during "
                    "j-generation");
        }

        choice = random_below(random_state, available);
        for (edge_index = 0; edge_index < edge_count; edge_index++) {
            if (edge_index == forbidden)
                continue;
            if (choice-- == 0) {
                selected = edge_index;
                break;
            }
        }
        if (selected == NO_INDEX)
            die("internal error selecting a random 2-isogeny edge");
        key_set(&previous, &current);
        key_set(&current, edges + selected);
        have_previous = 1;
        key_array_clear(edges, (size_t)phi_2->ell + 1);
        free(edges);
    }
    key_set(result, &current);
    key_clear(&current);
    key_clear(&previous);
}

static size_t rerandomization_tree_add(rerandomization_tree_t *tree,
                                       jkey_t key, size_t parent)
{
    size_t index;

    if (tree->length == tree->alloc) {
        if (tree->alloc > SIZE_MAX / 2)
            die("rerandomization tree is too large");
        tree->alloc = tree->alloc ? 2 * tree->alloc : 256;
        tree->nodes = xrealloc(tree->nodes, tree->alloc,
                               sizeof(*tree->nodes));
    }
    index = tree->length++;
    key_init(&tree->nodes[index].key);
    key_set(&tree->nodes[index].key, &key);
    tree->nodes[index].parent = parent;
    tree->nodes[index].depth =
        parent == NO_INDEX ? 0 : tree->nodes[parent].depth + 1;
    return index;
}

static void append_index(size_t **indices, size_t *length, size_t *alloc,
                         size_t index)
{
    if (*length == *alloc) {
        if (*alloc > SIZE_MAX / 2)
            die("rerandomization layer is too large");
        *alloc = *alloc ? 2 * *alloc : 256;
        *indices = xrealloc(*indices, *alloc, sizeof(**indices));
    }
    (*indices)[(*length)++] = index;
}

static void expand_rerandomization_layer(
    rerandomization_tree_t *tree,
    const size_t *current_layer, size_t current_length,
    size_t **next_layer_out, size_t *next_length_out,
    const modular_poly_t *phi, const field_context_t *field)
{
    size_t *next_layer = NULL;
    size_t next_length = 0, next_alloc = 0;
    size_t i;

    for (i = 0; i < current_length; i++) {
        size_t parent_index = current_layer[i];
        jkey_t parent_key = tree->nodes[parent_index].key;
        size_t grandparent = tree->nodes[parent_index].parent;
        jkey_t *edges = NULL;
        size_t edge_count, edge_index, children = 0;
        int removed_dual = 0;

        edge_count =
            modular_poly_edges_at(&edges, parent_key, phi, field);
        if (edge_count != (size_t)phi->ell + 1)
            die("rerandomization specialization did not split into ell'+1 "
                "F_{p^2}-rational edges");

        for (edge_index = 0; edge_index < edge_count; edge_index++) {
            size_t child_index;

            if (grandparent != NO_INDEX && !removed_dual &&
                key_equal(edges[edge_index],
                          tree->nodes[grandparent].key)) {
                removed_dual = 1;
                continue;
            }
            child_index =
                rerandomization_tree_add(tree, edges[edge_index],
                                         parent_index);
            append_index(&next_layer, &next_length, &next_alloc,
                         child_index);
            children++;
        }
        key_array_clear(edges, (size_t)phi->ell + 1);
        free(edges);

        if (grandparent != NO_INDEX && !removed_dual)
            die("could not identify the dual rerandomization edge");
        if (children != (size_t)phi->ell +
                            (grandparent == NO_INDEX ? 1 : 0))
            die("incorrect non-backtracking rerandomization branch count");
    }

    *next_layer_out = next_layer;
    *next_length_out = next_length;
}

static jkey_t *reconstruct_rerandomization_walk(
    const rerandomization_tree_t *tree, size_t endpoint,
    size_t *walk_length_out)
{
    size_t walk_length = tree->nodes[endpoint].depth + 1;
    jkey_t *walk = xrealloc(NULL, walk_length, sizeof(*walk));
    size_t position = walk_length;

    key_array_init(walk, walk_length);
    while (endpoint != NO_INDEX) {
        position--;
        key_set(walk + position, &tree->nodes[endpoint].key);
        endpoint = tree->nodes[endpoint].parent;
    }
    if (position != 0)
        die("invalid rerandomization parent chain");
    *walk_length_out = walk_length;
    return walk;
}

static void rerandomization_tree_clear(rerandomization_tree_t *tree)
{
    size_t i;

    for (i = 0; i < tree->length; i++)
        key_clear(&tree->nodes[i].key);
    free(tree->nodes);
    tree->nodes = NULL;
    tree->length = tree->alloc = 0;
}

int main(int argc, char **argv)
{
    jkey_t j1, target;
    jkey_t *successful_walk = NULL;
    unsigned B, ell, rerandomization_ell = 0, n_threads = 1;
    unsigned rerandomization_jobs = 1;
    uint64_t supplied_seed = 0, generation_seed = 0;
    const char *positional[8];
    const char *N_string = NULL;
    const char *phi_directory = DEFAULT_PHI_DIRECTORY;
    const char *p_string = NULL, *pbits_string = NULL, *B_string = NULL;
    const char *N_option_string = NULL, *phi_directory_option = NULL;
    const char *j1_re_string = NULL, *j1_im_string = NULL;
    const char *j2_re_string = NULL, *j2_im_string = NULL;
    modular_poly_t *polynomials;
    modular_poly_t rerandomization_phi;
    size_t n_polynomials = 0, k = 0, n_positional = 0;
    size_t rerandomizations_attempted = 0;
    size_t generation_walk_length = 0;
    size_t *current_layer = NULL;
    size_t current_length = 0, current_alloc = 0;
    size_t successful_walk_length = 0, successful_endpoint = NO_INDEX;
    size_t successful_attempt_number = 0;
    size_t rerandomization_depth = 0;
    rerandomization_tree_t rerandomization_tree = {0};
    field_context_t field;
    search_graph_t left, right;
    fmpz_t N, path_degree, quadratic_nonresidue;
    size_t left_meet = 0, right_meet = 0;
    int found, graphs_live = 0, attempted_rerandomization = 0;
    int rerandomization_phi_loaded = 0;
    int automatic_j = 0, seed_supplied = 0, threads_supplied = 0;
    int multipoint_batch_supplied = 0, force_multipoint_supplied = 0;
    int rerandomization_jobs_supplied = 0;
    int p_supplied = 0, pbits_supplied = 0;
    int B_supplied = 0, B_defaulted = 0;
    int N_supplied = 0, phi_directory_supplied = 0;
    int j1_supplied = 0, j2_supplied = 0;
    int rerandomization_supplied = 0, rerandomization_value = 0;
    int named_configuration = 0, legacy_arguments = 0, j_mode = 0;
    int allow_rerandomization;
    int argument_index;
    double total_start = wall_seconds();
    double setup_end, load_start, load_end;
    double generation_start = 0.0, generation_end = 0.0;
    double initial_start, initial_end;
    double rerandomization_start = 0.0, rerandomization_setup_end = 0.0;
    double rerandomization_end = 0.0, output_end;
    double cleanup_start, total_end;
    double peak_memory_mib;

    fmpz_init(field_characteristic);
    fmpz_init(N);
    fmpz_init(path_degree);
    fmpz_init(quadratic_nonresidue);

    for (argument_index = 1; argument_index < argc; argument_index++) {
        if (strcmp(argv[argument_index], "--threads") == 0 ||
            strcmp(argv[argument_index], "-t") == 0) {
            if (threads_supplied || argument_index + 1 == argc)
                die("--threads must occur once and have a value");
            n_threads = parse_thread_count(argv[++argument_index]);
            threads_supplied = 1;
        } else if (strcmp(argv[argument_index],
                          "--multipoint-batch") == 0 ||
                   strcmp(argv[argument_index],
                          "--batch-size") == 0) {
            if (multipoint_batch_supplied ||
                argument_index + 1 == argc)
                die("--multipoint-batch must occur once and have a value");
            ordered_multipoint_batch_limit =
                parse_multipoint_batch_limit(
                    argv[++argument_index]);
            multipoint_batch_supplied = 1;
        } else if (strcmp(argv[argument_index],
                          "--force-multipoint") == 0) {
            if (force_multipoint_supplied)
                die("--force-multipoint must occur at most once");
            ordered_force_multipoint = 1;
            force_multipoint_supplied = 1;
        } else if (strcmp(argv[argument_index],
                          "--rerandomization-jobs") == 0 ||
                   strcmp(argv[argument_index],
                          "--parallel-rerandomizations") == 0) {
            if (rerandomization_jobs_supplied ||
                argument_index + 1 == argc)
                die("--rerandomization-jobs must occur once and have a "
                    "value");
            rerandomization_jobs =
                parse_rerandomization_job_count(
                    argv[++argument_index]);
            rerandomization_jobs_supplied = 1;
        } else if (strcmp(argv[argument_index], "--seed") == 0) {
            if (seed_supplied || argument_index + 1 == argc)
                die("--seed must occur once and have a value");
            supplied_seed = parse_seed(argv[++argument_index]);
            seed_supplied = 1;
        } else if (strcmp(argv[argument_index], "--p") == 0) {
            if (p_supplied || argument_index + 1 == argc)
                die("--p must occur once and have a value");
            p_string = argv[++argument_index];
            p_supplied = 1;
            named_configuration = 1;
        } else if (strcmp(argv[argument_index], "--pbits") == 0) {
            if (pbits_supplied || argument_index + 1 == argc)
                die("--pbits must occur once and have a value");
            pbits_string = argv[++argument_index];
            pbits_supplied = 1;
            named_configuration = 1;
        } else if (strcmp(argv[argument_index], "--B") == 0) {
            if (B_supplied || argument_index + 1 == argc)
                die("--B must occur once and have a value");
            B_string = argv[++argument_index];
            B_supplied = 1;
            named_configuration = 1;
        } else if (strcmp(argv[argument_index], "--N") == 0) {
            if (N_supplied || argument_index + 1 == argc)
                die("--N must occur once and have a value");
            N_option_string = argv[++argument_index];
            N_supplied = 1;
            named_configuration = 1;
        } else if (strcmp(argv[argument_index], "--phi-dir") == 0 ||
                   strcmp(argv[argument_index],
                          "--phi-directory") == 0) {
            if (phi_directory_supplied ||
                argument_index + 1 == argc)
                die("--phi-dir must occur once and have a value");
            phi_directory_option = argv[++argument_index];
            phi_directory_supplied = 1;
            named_configuration = 1;
        } else if (strcmp(argv[argument_index], "--j1") == 0) {
            if (j1_supplied || argument_index + 2 >= argc)
                die("--j1 must occur once and have two values");
            j1_re_string = argv[++argument_index];
            j1_im_string = argv[++argument_index];
            j1_supplied = 1;
            named_configuration = 1;
        } else if (strcmp(argv[argument_index], "--j2") == 0) {
            if (j2_supplied || argument_index + 2 >= argc)
                die("--j2 must occur once and have two values");
            j2_re_string = argv[++argument_index];
            j2_im_string = argv[++argument_index];
            j2_supplied = 1;
            named_configuration = 1;
        } else if (strcmp(argv[argument_index],
                          "--rerandomization") == 0 ||
                   strcmp(argv[argument_index],
                          "--rerandomize") == 0) {
            if (rerandomization_supplied ||
                argument_index + 1 == argc)
                die("--rerandomization must occur once and have a value");
            rerandomization_value =
                parse_on_off(argv[++argument_index]);
            rerandomization_supplied = 1;
        } else {
            if (argv[argument_index][0] == '-')
                die("unknown command-line option");
            if (n_positional == sizeof(positional) / sizeof(positional[0]))
                die("too many positional arguments");
            positional[n_positional++] = argv[argument_index];
        }
    }
    legacy_arguments =
        !named_configuration &&
        (n_positional == 4 || n_positional == 6 ||
         n_positional == 8);
    if (!legacy_arguments && n_positional > 2) {
        fprintf(stderr,
                "usage: %s [N [phi_directory]] "
                "[--p p | --pbits n] [--B B]\n"
                "          [--N N] [--phi-dir directory]\n"
                "          [--j1 re im [--j2 re im]] "
                "[--rerandomization on|off]\n"
                "          [--threads n] [--multipoint-batch n] "
                "[--force-multipoint]\n"
                "          "
                "[--rerandomization-jobs n] [--seed s]\n"
                "legacy: %s p B N phi_directory [options]\n"
                "        %s p j1_re j1_im B N phi_directory [options]\n"
                "        %s p j1_re j1_im j2_re j2_im B N "
                "phi_directory [options]\n",
                argv[0], argv[0], argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    if (!legacy_arguments) {
        if (j2_supplied && !j1_supplied)
            die("--j2 requires --j1");
        if (p_supplied && pbits_supplied)
            die("--p and --pbits are mutually exclusive");
        if (N_supplied && n_positional >= 1)
            die("specify N either positionally or with --N, not both");
        if (phi_directory_supplied && n_positional >= 2)
            die("specify the modular-polynomial directory either "
                "positionally or with --phi-dir, not both");
        if (!p_supplied && !pbits_supplied)
            p_string = DEFAULT_CHARACTERISTIC;
        if (B_supplied)
            B = parse_bound(B_string);
        else {
            B = 0;
            B_defaulted = 1;
        }
        if (N_supplied)
            N_string = N_option_string;
        else if (n_positional >= 1)
            N_string = positional[0];
        if (phi_directory_supplied)
            phi_directory = phi_directory_option;
        else if (n_positional >= 2)
            phi_directory = positional[1];
        j_mode = j2_supplied ? 2 : j1_supplied ? 1 : 0;
    } else {
        p_string = positional[0];
        if (n_positional == 4) {
            B = parse_bound(positional[1]);
            N_string = positional[2];
            phi_directory = positional[3];
            j_mode = 0;
        } else if (n_positional == 6) {
            j1_re_string = positional[1];
            j1_im_string = positional[2];
            B = parse_bound(positional[3]);
            N_string = positional[4];
            phi_directory = positional[5];
            j_mode = 1;
        } else {
            j1_re_string = positional[1];
            j1_im_string = positional[2];
            j2_re_string = positional[3];
            j2_im_string = positional[4];
            B = parse_bound(positional[5]);
            N_string = positional[6];
            phi_directory = positional[7];
            j_mode = 2;
        }
    }

    if (pbits_supplied) {
        ulong pbits = parse_characteristic_bits(pbits_string);

        largest_prime_3mod4_below_power_of_two(
            field_characteristic, pbits);
    } else {
        parse_characteristic(field_characteristic, p_string);
    }
    active_backend = fmpz_abs_fits_ui(field_characteristic)
                         ? BACKEND_FQ_NMOD
                         : BACKEND_FQ;
    key_init(&j1);
    key_init(&target);
    if (j_mode >= 1)
        key_set_strings(&j1, j1_re_string, j1_im_string);
    if (j_mode == 2)
        key_set_strings(&target, j2_re_string, j2_im_string);

    automatic_j = j_mode == 0;
    allow_rerandomization =
        rerandomization_supplied ? rerandomization_value : j_mode != 2;
    if (!allow_rerandomization && rerandomization_jobs_supplied)
        die("--rerandomization-jobs requires rerandomization to be on");
    if (!automatic_j && seed_supplied)
        die("--seed is only valid when no j-invariant is supplied");
    if (automatic_j) {
        if (fmpz_fdiv_ui(field_characteristic, UWORD(4)) != UWORD(3))
            die("automatic j-generation from 1728 requires p = 3 (mod 4)");
        generation_walk_length = 2 * ceil_log2_characteristic();
    }
    flint_set_num_threads((int)n_threads);

    if (N_string != NULL) {
        if (fmpz_set_str(N, N_string, 10) != 0 ||
            fmpz_sgn(N) <= 0)
            die("N must be a positive decimal integer");
    } else {
        fmpz_t half_p;

        fmpz_init_set(half_p, field_characteristic);
        fmpz_fdiv_q_2exp(half_p, half_p, 1);
        (void)fmpz_root(N, half_p, 3);
        fmpz_clear(half_p);
    }
    if (B_defaulted)
        B = default_bound_from_N(N);
    if (allow_rerandomization)
        rerandomization_ell = smallest_prime_larger_than(B);
    if (fmpz_cmp_ui(field_characteristic, (ulong)B) <= 0)
        die("p must be greater than B so every ell-isogeny degree is "
            "coprime to p");
    if (allow_rerandomization &&
        fmpz_equal_ui(field_characteristic,
                      (ulong)rerandomization_ell))
        die("the rerandomization degree ell' must be different from p");

    choose_quadratic_nonresidue(quadratic_nonresidue,
                                field_characteristic);
    field_context_init(&field, field_characteristic,
                       quadratic_nonresidue);
    setup_end = wall_seconds();

    for (ell = 2; ell <= B; ell++)
        if (is_prime(ell))
            n_polynomials++;
    polynomials = xrealloc(NULL, n_polynomials, sizeof(*polynomials));

    load_start = wall_seconds();
    for (ell = 2; ell <= B; ell++) {
        if (!is_prime(ell))
            continue;
        modular_poly_init(polynomials + k, ell, &field);
        modular_poly_load(polynomials + k, phi_directory, &field);
        k++;
    }
    load_end = wall_seconds();
    fprintf(stderr,
            "Loaded %lu smooth modular polynomials once, reduced modulo p, "
            "in %.3f s\n",
            (unsigned long)n_polynomials, load_end - load_start);

    if (automatic_j) {
        random_state_t random_state;

        generation_start = wall_seconds();
        generation_seed =
            random_state_init(&random_state, seed_supplied, supplied_seed);
        generate_random_start_j(&j1, polynomials, &field, &random_state,
                                generation_walk_length);
        generation_end = wall_seconds();
        fprintf(stderr,
                "Generated j by a %lu-step non-backtracking 2-isogeny "
                "walk from 1728\n",
                (unsigned long)generation_walk_length);
        fprintf(stderr, "Random-generation seed = %llu\n",
                (unsigned long long)generation_seed);
    }
    if (j_mode != 2)
        frobenius_conjugate(&target, &j1);

    fprintf(stderr, "Field characteristic p = ");
    fmpz_fprint(stderr, field_characteristic);
    if (pbits_supplied)
        fprintf(stderr,
                " (largest prime below 2^%s congruent to 3 modulo 4)",
                pbits_string);
    fprintf(stderr, "\nArithmetic backend = %s\n",
            active_backend == BACKEND_FQ_NMOD
                ? "fq_nmod (single-word characteristic)"
                : "fq (arbitrary-precision characteristic)");
    fprintf(stderr, "Degree bound N = ");
    fmpz_fprint(stderr, N);
    fprintf(stderr, "%s\n",
            N_string == NULL
                ? " (default floor(cuberoot(floor(p/2))))"
                : "");
    fprintf(stderr, "Modular-polynomial directory = %s%s\n",
            phi_directory,
            (phi_directory_supplied || n_positional >= 2 ||
             legacy_arguments)
                ? ""
                : " (default)");
    if (fmpz_cmp_ui(quadratic_nonresidue, 0) > 0) {
        fmpz_t p_minus_one;

        fmpz_init(p_minus_one);
        fmpz_sub_ui(p_minus_one, field_characteristic, 1);
        if (fmpz_equal(quadratic_nonresidue, p_minus_one))
            fprintf(stderr, "Field basis: F_p[a]/(a^2 + 1)\n");
        else {
            fprintf(stderr, "Field basis: F_p[a]/(a^2 - ");
            fmpz_fprint(stderr, quadratic_nonresidue);
            fprintf(stderr, ")\n");
        }
        fmpz_clear(p_minus_one);
    } else {
        die("internal error constructing the quadratic nonresidue");
    }
    fprintf(stderr, "%s j1 = ",
            automatic_j ? "Generated" : "Input");
    fprint_key(stderr, j1);
    fprintf(stderr, "\n");
    fprintf(stderr, "%s = ",
            j_mode == 2 ? "Input j2" : "Frobenius target j1^p");
    fprint_key(stderr, target);
    fprintf(stderr, "\n");
    if (allow_rerandomization)
        fprintf(stderr,
                "Smoothness bound B = %u%s; rerandomization = on "
                "(ell' = %u)%s\n",
                B, B_defaulted ? " [default floor(N^(1/4))]" : "",
                rerandomization_ell,
                rerandomization_supplied ? " [explicit]" : " [default]");
    else
        fprintf(stderr,
                "Smoothness bound B = %u%s; rerandomization = off%s\n",
                B, B_defaulted ? " [default floor(N^(1/4))]" : "",
                rerandomization_supplied ? " [explicit]" : " [default]");
    fprintf(stderr, "Worker threads = %u\n", n_threads);
    if (ordered_multipoint_batch_limit == 0)
        fprintf(stderr,
                "Multipoint batch limit = unlimited "
                "(complete eligible prime-power layer)\n");
    else
        fprintf(stderr, "Multipoint batch limit = %lu%s\n",
                (unsigned long)ordered_multipoint_batch_limit,
                multipoint_batch_supplied ? " [explicit]" : "");
    fprintf(stderr, "Force multipoint evaluation = %s%s\n",
            ordered_force_multipoint ? "on" : "off",
            force_multipoint_supplied ? " [explicit]" : " [default]");
    if (allow_rerandomization)
        fprintf(stderr,
                "Concurrent rerandomization jobs = %u%s\n",
                rerandomization_jobs,
                rerandomization_jobs_supplied ? " [explicit]"
                                              : " [default]");

    initial_start = wall_seconds();
    found = ordered_search_between(
        &left, &right, j1, target, N,
        polynomials, n_polynomials, &field, n_threads,
        "j1", j_mode == 2 ? "j2" : "j1^p",
        &left_meet, &right_meet, path_degree, 1, NULL, NULL);
    graphs_live = 1;
    initial_end = wall_seconds();

    if (found) {
        print_path(&left, left_meet, &right, right_meet, path_degree);
    } else if (!allow_rerandomization) {
        printf("No B-smooth isogeny from j1 to %s of degree at most ",
               j_mode == 2 ? "j2" : "j1^p");
        fmpz_fprint(stdout, N);
        printf(" was found.\n");
    } else {
        printf("Initial B-smooth search failed; rerandomizing with "
               "unbounded non-backtracking ell' = %u walks.\n",
               rerandomization_ell);
        fflush(stdout);

        rerandomization_start = wall_seconds();
        modular_poly_init(&rerandomization_phi, rerandomization_ell,
                          &field);
        modular_poly_load(&rerandomization_phi, phi_directory, &field);
        rerandomization_phi_loaded = 1;
        append_index(&current_layer, &current_length, &current_alloc,
                     rerandomization_tree_add(&rerandomization_tree, j1,
                                              NO_INDEX));
        rerandomization_setup_end = wall_seconds();

        fprintf(stderr, "Middle smooth-path degree bound N = ");
        fmpz_fprint(stderr, N);
        fprintf(stderr, "\n");

        graph_clear(&left);
        graph_clear(&right);
        graphs_live = 0;

        found = 0;
        for (;;) {
            size_t *next_layer = NULL;
            size_t next_length = 0;
            size_t rerandomization_index;

            expand_rerandomization_layer(
                &rerandomization_tree, current_layer, current_length,
                &next_layer, &next_length, &rerandomization_phi, &field);
            free(current_layer);
            current_layer = next_layer;
            current_length = next_length;
            current_alloc = next_length;
            rerandomization_depth++;

            fprintf(stderr,
                    "\nRerandomization depth %lu: %lu "
                    "non-backtracking walk%s\n",
                    (unsigned long)rerandomization_depth,
                    (unsigned long)current_length,
                    current_length == 1 ? "" : "s");

            for (rerandomization_index = 0;
                 rerandomization_index < current_length && !found;) {
                size_t batch_count =
                    current_length - rerandomization_index;
                size_t job_index;
                rerandomization_search_job_t *jobs;
                pthread_t *threads = NULL;
                atomic_int cancel_requested;
                atomic_int winner_slot;
                int winner;

                if (batch_count > (size_t)rerandomization_jobs)
                    batch_count = (size_t)rerandomization_jobs;
                jobs = xrealloc(NULL, batch_count, sizeof(*jobs));
                if (batch_count > 1)
                    threads = xrealloc(NULL, batch_count,
                                       sizeof(*threads));
                atomic_init(&cancel_requested, 0);
                atomic_init(&winner_slot, -1);

                fprintf(stderr,
                        "\nStarting rerandomization batch of %lu "
                        "candidate%s (walks %lu--%lu of %lu)\n",
                        (unsigned long)batch_count,
                        batch_count == 1 ? "" : "s",
                        (unsigned long)rerandomization_index + 1,
                        (unsigned long)(rerandomization_index +
                                        batch_count),
                        (unsigned long)current_length);

                for (job_index = 0; job_index < batch_count;
                     job_index++) {
                    rerandomization_search_job_t *job =
                        jobs + job_index;

                    memset(job, 0, sizeof(*job));
                    job->endpoint =
                        current_layer[rerandomization_index +
                                      job_index];
                    job->attempt_number =
                        ++rerandomizations_attempted;
                    job->layer_position =
                        rerandomization_index + job_index;
                    key_init(&job->start);
                    key_init(&job->target);
                    key_set(
                        &job->start,
                        &rerandomization_tree.nodes[job->endpoint].key);
                    if (j_mode == 2)
                        key_set(&job->target, &target);
                    else
                        frobenius_conjugate(&job->target,
                                            &job->start);
                    fmpz_init(job->path_degree);
                    job->degree_bound = N;
                    job->polynomials = polynomials;
                    job->n_polynomials = n_polynomials;
                    job->field = &field;
                    job->search_threads =
                        rerandomization_jobs == 1 ? n_threads : 1;
                    job->right_label =
                        j_mode == 2 ? "j2" : "(j')^p";
                    job->cancel_requested = &cancel_requested;
                    job->winner_slot = &winner_slot;
                    job->slot = (int)job_index;

                    fprintf(stderr,
                            "  attempt %lu, walk %lu: j' = ",
                            (unsigned long)job->attempt_number,
                            (unsigned long)job->layer_position + 1);
                    fprint_key(stderr, job->start);
                    fprintf(stderr, ", target = ");
                    fprint_key(stderr, job->target);
                    fprintf(stderr, "\n");
                }

                if (batch_count == 1) {
                    run_rerandomization_search(jobs);
                } else {
                    size_t created = 0;
                    int create_failed = 0;

                    for (job_index = 0; job_index < batch_count;
                         job_index++) {
                        if (pthread_create(threads + job_index, NULL,
                                           run_rerandomization_search,
                                           jobs + job_index) != 0) {
                            create_failed = 1;
                            atomic_store_explicit(
                                &cancel_requested, 1,
                                memory_order_relaxed);
                            break;
                        }
                        created++;
                    }
                    for (job_index = 0; job_index < created;
                         job_index++)
                        if (pthread_join(threads[job_index], NULL) != 0)
                            die("pthread_join failed during "
                                "rerandomization");
                    if (create_failed)
                        die("pthread_create failed during "
                            "rerandomization");
                }

                attempted_rerandomization = 1;
                winner = atomic_load_explicit(
                    &winner_slot, memory_order_relaxed);
                for (job_index = 0; job_index < batch_count;
                     job_index++) {
                    rerandomization_search_job_t *job =
                        jobs + job_index;
                    const char *status;

                    if ((int)job_index == winner)
                        status = "succeeded";
                    else if (job->found)
                        status = "also found a path";
                    else if (job->cancelled)
                        status = "cancelled after another success";
                    else
                        status = "failed";
                    fprintf(stderr,
                            "Rerandomization attempt %lu %s in "
                            "%.3f s (left ball %lu, right ball %lu)\n",
                            (unsigned long)job->attempt_number,
                            status, job->elapsed_seconds,
                            (unsigned long)job->left.length,
                            (unsigned long)job->right.length);
                }

                if (winner >= 0) {
                    rerandomization_search_job_t *winning_job =
                        jobs + winner;

                    found = 1;
                    successful_endpoint = winning_job->endpoint;
                    successful_attempt_number =
                        winning_job->attempt_number;
                    left = winning_job->left;
                    right = winning_job->right;
                    left_meet = winning_job->left_meet;
                    right_meet = winning_job->right_meet;
                    fmpz_set(path_degree,
                             winning_job->path_degree);
                    graphs_live = 1;
                }

                for (job_index = 0; job_index < batch_count;
                     job_index++) {
                    rerandomization_search_job_t *job =
                        jobs + job_index;

                    if ((int)job_index != winner) {
                        graph_clear(&job->left);
                        graph_clear(&job->right);
                    }
                    fmpz_clear(job->path_degree);
                    key_clear(&job->target);
                    key_clear(&job->start);
                }
                free(threads);
                free(jobs);
                rerandomization_index += batch_count;
            }
            if (found)
                break;
        }
        rerandomization_end = wall_seconds();
        successful_walk = reconstruct_rerandomization_walk(
            &rerandomization_tree, successful_endpoint,
            &successful_walk_length);
        printf("Rerandomization attempt %lu succeeded at depth %lu "
               "(%lu candidate%s started).\n",
               (unsigned long)successful_attempt_number,
               (unsigned long)rerandomization_depth,
               (unsigned long)rerandomizations_attempted,
               rerandomizations_attempted == 1 ? "" : "s");
        print_rerandomized_path(
            &left, left_meet, &right, right_meet, successful_walk,
            successful_walk_length, rerandomization_ell, path_degree,
            j_mode != 2);
    }
    fflush(stdout);
    output_end = wall_seconds();

    cleanup_start = wall_seconds();
    if (graphs_live) {
        graph_clear(&left);
        graph_clear(&right);
    }
    for (k = 0; k < n_polynomials; k++)
        modular_poly_clear(polynomials + k, &field);
    if (rerandomization_phi_loaded)
        modular_poly_clear(&rerandomization_phi, &field);
    rerandomization_tree_clear(&rerandomization_tree);
    free(current_layer);
    key_array_clear(successful_walk, successful_walk_length);
    free(successful_walk);
    free(polynomials);
    key_clear(&target);
    key_clear(&j1);
    field_context_clear(&field);
    fmpz_clear(quadratic_nonresidue);
    fmpz_clear(N);
    fmpz_clear(path_degree);
    fmpz_clear(field_characteristic);
    flint_cleanup();
    total_end = wall_seconds();
    peak_memory_mib = peak_resident_memory_mib();

    printf("\nTiming (wall clock, %u FLINT worker thread%s",
           n_threads, n_threads == 1 ? "" : "s");
    if (allow_rerandomization)
        printf(", %u concurrent rerandomization job%s",
               rerandomization_jobs,
               rerandomization_jobs == 1 ? "" : "s");
    printf(")\n");
    printf("  setup:                  %.3f s\n", setup_end - total_start);
    printf("  modular polynomial load: %.3f s\n", load_end - load_start);
    if (automatic_j)
        printf("  random j generation:    %.3f s\n",
               generation_end - generation_start);
    printf("  initial search:         %.3f s\n", initial_end - initial_start);
    if (rerandomization_start != 0.0) {
        printf("  rerandomization setup:  %.3f s\n",
               rerandomization_setup_end - rerandomization_start);
    }
    if (attempted_rerandomization)
        printf("  rerandomized searches:  %.3f s "
               "(%lu candidate%s started)\n",
               rerandomization_end - rerandomization_setup_end,
               (unsigned long)rerandomizations_attempted,
               rerandomizations_attempted == 1 ? "" : "s");
    printf("  path output:            %.3f s\n",
           output_end - (attempted_rerandomization
                             ? rerandomization_end
                             : rerandomization_start != 0.0
                                   ? rerandomization_setup_end
                                   : initial_end));
    printf("  cleanup:                %.3f s\n", total_end - cleanup_start);
    printf("  total runtime:          %.3f s\n", total_end - total_start);
    if (peak_memory_mib >= 0.0)
        printf("  peak resident memory:   %.2f MiB\n", peak_memory_mib);
    else
        printf("  peak resident memory:   unavailable\n");

    return found ? EXIT_SUCCESS : 2;
}
