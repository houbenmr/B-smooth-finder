#define _POSIX_C_SOURCE 200809L

/*
 * Prime-ordered smooth-isogeny search between supersingular j-invariants.
 * With zero or one input j-invariant, the target is the p-power Frobenius
 * conjugate. In these modes, only the ball from j is stored and completely
 * enumerated; its conjugate ball is represented implicitly, and collisions
 * with it are tested only after enumeration finishes. With two supplied
 * j-invariants, two balls are enumerated alternately and checked online.
 * Rerandomization defaults to enabled with zero or one supplied j-invariant
 * and disabled with two supplied j-invariants; the default can be overridden
 * explicitly.
 *
 * Field: F_{p^2} = F_p[a]/(a^2 - q), where p is an input prime smaller
 * than 2^127 and q is a deterministic quadratic nonresidue modulo p. For
 * p = 3 (mod 4), q = -1, preserving the representation
 * F_p[a]/(a^2 + 1). An element is entered and printed as (real, imag),
 * meaning real + imag*a.
 *
 * The program selects one of two FLINT backends at runtime:
 *   - fq_nmod for p values that fit in one machine word;
 *   - fq for larger p values represented by fmpz.
 *
 * Canonical usage:
 *   ./ordered_prime_search_flint34 \
 *       [N [phi_directory]] [options]
 *
 * Options:
 *   --p p         field characteristic (default: 4294967311)
 *   --pbits n     instead choose the largest prime p < 2^n with p = 3
 *                 (mod 4), 2 <= n <= 127; mutually exclusive with --p
 *   --B B         smoothness bound, 2 <= B <= 10000
 *                 (default: floor(N^(1/4)))
 *   --N N         middle smooth-path degree bound
 *                 (default: floor(cuberoot(floor(p/2))))
 *   --phi-dir d   modular-polynomial directory (default: mod_pols)
 *   --j1 re im    use the supplied starting j-invariant
 *   --j2 re im    use this target; requires --j1
 *   --rerandomization on|off
 *                 override the mode-dependent rerandomization default
 *   --threads n   total FLINT worker budget, 1 <= n <= 256
 *   --multipoint-batch n
 *                 cap one multipoint batch at n source points; n = 0
 *                 means no cap (default: no cap)
 *   --evaluation-shards n
 *                 use at most n independent point shards for evaluation;
 *                 n = 0 selects automatically from --threads (default)
 *   --force-multipoint
 *                 use a FLINT multipoint tree for every nonempty
 *                 ordered-search batch (default: off)
 *   --rerandomization-jobs n
 *                 test up to n rerandomization endpoints concurrently,
 *                 1 <= n <= 256 (default: 1)
 *   --root-profile report root-extraction CPU time by residual degree
 *   --parallel-factor-tree on|off
 *                 cooperatively split underfilled high-degree root batches
 *                 across factor-tree levels (default: on)
 *   --cubic-kernel flint|batched|compare
 *                 select the current FLINT Rabin path, an experimental
 *                 fixed-degree batched cubic path, or benchmark both on
 *                 every actual cubic specialization (default: batched)
 *   --seed s      seed automatic j-generation with the 64-bit integer s
 *
 * Legacy positional forms remain accepted:
 *   ./ordered_prime_search_flint34 \
 *       p B N phi_directory \
 *       [options]
 *   ./ordered_prime_search_flint34 \
 *       p j1_re j1_im B N \
 *       phi_directory [options]
 *   ./ordered_prime_search_flint34 \
 *       p j1_re j1_im \
 *       j2_re j2_im B N \
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
 *       ordered_prime_search_flint34.c \
 *       -pthread -lflint -lmpfr -lgmp \
 *       -o ordered_prime_search_flint34
 *
 * Performance design:
 *   - Every Phi_ell is read exactly once. Its coefficient polynomials in X
 *     are reduced modulo p and retained as fq_nmod_poly or fq_poly objects.
 *   - The balls are enumerated in descending canonical prime-factor order.
 *     Expensive high-degree Phi_ell root extractions therefore run while
 *     the frontier is still small; the larger later frontiers use cheaper
 *     low-degree modular polynomials. Coprime isogeny diamonds allow the
 *     prime factors of every smooth path to be placed in this fixed order.
 *   - At an ell-stage, endpoints of different current degrees are evaluated
 *     together with the same Phi_ell, subject to the individual condition
 *     degree*ell <= floor(sqrt(N)). Consecutive ell-power layers remain
 *     sequential because one layer supplies the points for the next.
 *   - In the zero- and one-j modes, one complete ball is enumerated from j.
 *     A post-enumeration scan looks for x and x^p in that same graph; this is
 *     exactly an intersection with the implicit conjugate ball.
 *   - In the two-j mode, the two explicit balls alternate after every batch.
 *     Newly discovered vertices are checked immediately against the other
 *     ball.
 *   - By default, a complete eligible ell-layer forms one batch. The
 *     optional --multipoint-batch n cap trades batching efficiency for
 *     lower temporary memory and more frequent intersection checks.
 *   - By default, smaller batches use Horner evaluation, since a multipoint
 *     tree costs more than it saves at that size. --force-multipoint
 *     disables this crossover for ordered-search batches.
 *   - Evaluation is two-dimensionally tiled over both coefficient
 *     polynomials and source points. Fast multipoint batches are split into
 *     independent product-tree shards, whose trees are built, evaluated and
 *     freed in parallel. This avoids the old ell+2-worker ceiling (only four
 *     evaluation workers for Phi_2) and parallelizes product-tree
 *     construction on many-core machines.
 *   - In one-ball Frobenius mode with high-degree modular polynomials, the
 *     root-j specializations for all primes are computed concurrently,
 *     largest primes first, and cached for their respective stages. This
 *     removes the otherwise serial first high-degree root at every stage.
 *   - FLINT's worker pool factors independent specializations in parallel.
 *     Root extraction uses
 *     an adaptive atomic work queue because FLINT 3.x declares, but does not
 *     implement, dynamic scheduling in flint_parallel_do. Each worker owns
 *     and reuses its exponent, random state, factor stack, field scratch, and
 *     polynomial scratch for a complete batch; there are no locks in the
 *     root hot path. This remains suitable for thread counts up to 256.
 *     Underfilled high-degree batches also expose successive Rabin
 *     factor-tree levels to the shared pool; this can be disabled with
 *     --parallel-factor-tree off.
 *   - When the current path arrived by an ell-isogeny, the known dual root
 *     is divided out before root extraction. Linear and quadratic residual
 *     polynomials are solved directly. Higher degrees use a persistent Rabin
 *     splitter assembled from public FLINT 3.4 polynomial operations, so a
 *     worker does not rebuild FLINT_BITS temporary polynomials for every
 *     specialization.
 *   - For 65--127-bit p with p = 3 (mod 4), graph keys use fixed two-limb
 *     coordinates and quadratic leaves use FLINT's optimized two-limb
 *     fmpz_mod arithmetic directly. Residual cubics use an eight-lane
 *     projective splitter with batched inversion. Other field bases retain
 *     the generic fq operations.
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
 *     --rerandomization-jobs. This option controls the number of resident
 *     balls, not a second thread pool. At each prime-power wave, source
 *     points from every resident ball are interleaved into shared evaluation
 *     and root batches. Thus --threads is one total worker budget and remains
 *     fully available even when several rerandomizations are active.
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
#include <stdint.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <time.h>

/*
 * FLINT 3.4 deprecates the historical unseparated RNG names. Retain source
 * compatibility with FLINT 3.0 while compiling warning-free against 3.4.
 */
#if __FLINT_RELEASE >= __FLINT_RELEASE_NUM(3, 4, 0)
#define ISOGENY_RAND_INIT(state) flint_rand_init(state)
#define ISOGENY_RAND_CLEAR(state) flint_rand_clear(state)
#define ISOGENY_RAND_SET_SEED(state, seed1, seed2) \
    flint_rand_set_seed((state), (seed1), (seed2))
#else
#define ISOGENY_RAND_INIT(state) flint_randinit(state)
#define ISOGENY_RAND_CLEAR(state) flint_randclear(state)
#define ISOGENY_RAND_SET_SEED(state, seed1, seed2) \
    flint_randseed((state), (seed1), (seed2))
#endif

/*
 * FLINT 3.4 no longer exports the two private fq[_nmod] Rabin splitters.
 * This version reconstructs the persistent splitter from public operations;
 * the macro only retains compatibility with older FLINT installations.
 */
#ifndef ISOGENY_USE_PRIVATE_RABIN
#if __FLINT_RELEASE < __FLINT_RELEASE_NUM(3, 4, 0)
#define ISOGENY_USE_PRIVATE_RABIN 1
#else
#define ISOGENY_USE_PRIVATE_RABIN 0
#endif
#endif

#define NO_INDEX ((size_t)-1)
/* Retained only by the unused reference Dijkstra implementation below. */
#define MULTIPOINT_BATCH_SIZE 256
#define MULTIPOINT_THRESHOLD 24
#define MULTIPOINT_MIN_DEGREE 12
#define DEFAULT_CHARACTERISTIC "4294967311"
#define DEFAULT_PHI_DIRECTORY "mod_pols"
#define MAX_SMOOTHNESS_BOUND 10000U
#define MAX_CHARACTERISTIC_BITS 127U
#define MAX_RERANDOMIZATION_JOBS 256
#define ROOT_FACTOR_STACK_SIZE (FLINT_BITS + 3)
#define ROOT_QUEUE_CHUNKS_PER_WORKER 64
#define ROOT_QUEUE_MAX_CLAIM 64
#define MIN_QUADRATIC_ROOTS_PER_TASK 64
#define EVALUATION_MIN_POINTS_PER_SHARD 512
#define HORNER_MIN_POINTS_PER_SHARD 128
#define EVALUATION_TASKS_PER_WORKER 4
#define PARALLEL_FACTOR_MIN_DEGREE 24
#define PARALLEL_LIFECYCLE_MIN_THREADS 16
#define PARALLEL_LIFECYCLE_MIN_VALUES 4096
#define ROOT_PROFILE_DEGREES (MAX_SMOOTHNESS_BOUND + 2U)
#ifndef CUBIC_BATCH_LANES
#define CUBIC_BATCH_LANES 8
#endif
#define CUBIC_MAX_SPLIT_ATTEMPTS 64

/*
 * FLINT releases before 3.4 export these two low-level Cantor--Zassenhaus
 * splitters but do not declare them in public headers. Calling them lets a
 * worker retain half the field order, its random state, and all scratch
 * polynomials across thousands of specializations. FLINT 3.4 removed the
 * symbols, so the declarations and calls are compiled out there.
 */
#if ISOGENY_USE_PRIVATE_RABIN
void _fq_nmod_poly_split_rabin(
    fq_nmod_poly_t a, fq_nmod_poly_t b, const fq_nmod_poly_t f,
    const fmpz_t halfq, fq_nmod_poly_t t, fq_nmod_poly_t t2,
    flint_rand_t randstate, const fq_nmod_ctx_t ctx);

void _fq_poly_split_rabin(
    fq_poly_t a, fq_poly_t b, const fq_poly_t f,
    const fmpz_t halfq, fq_poly_t t, fq_poly_t t2,
    flint_rand_t randstate, const fq_ctx_t ctx);
#endif

typedef enum {
    BACKEND_FQ_NMOD,
    BACKEND_FQ
} backend_kind_t;

typedef enum {
    CUBIC_KERNEL_FLINT,
    CUBIC_KERNEL_BATCHED,
    CUBIC_KERNEL_COMPARE
} cubic_kernel_kind_t;

typedef struct {
    fmpz_t scratch[14];
    fmpz_t zero;
    fmpz_t sqrt_exponent;
    fmpz_t inverse_two;
} fast_fp2_workspace_t;

static void fast_fp2_workspace_init(
    fast_fp2_workspace_t *workspace, const fq_ctx_t ctx)
{
    size_t index;

    for (index = 0; index < 14; index++)
        fmpz_init(workspace->scratch[index]);
    fmpz_init(workspace->zero);
    fmpz_init(workspace->sqrt_exponent);
    fmpz_add_ui(workspace->sqrt_exponent, fq_ctx_prime(ctx), 1);
    fmpz_fdiv_q_2exp(
        workspace->sqrt_exponent, workspace->sqrt_exponent, 2);
    fmpz_init(workspace->inverse_two);
    fmpz_add_ui(workspace->inverse_two, fq_ctx_prime(ctx), 1);
    fmpz_fdiv_q_2exp(
        workspace->inverse_two, workspace->inverse_two, 1);
}

static void fast_fp2_workspace_clear(
    fast_fp2_workspace_t *workspace)
{
    size_t index;

    fmpz_clear(workspace->inverse_two);
    fmpz_clear(workspace->sqrt_exponent);
    fmpz_clear(workspace->zero);
    for (index = 14; index-- > 0;)
        fmpz_clear(workspace->scratch[index]);
}

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
static int quadratic_basis_is_minus_one;
static cubic_kernel_kind_t cubic_kernel = CUBIC_KERNEL_BATCHED;
static int root_profile_enabled = 0;

/* Summed thread CPU time; wall timers below cover isolated cubic passes. */
static atomic_uint_fast64_t root_profile_counts[ROOT_PROFILE_DEGREES];
static atomic_uint_fast64_t root_profile_cpu_ns[ROOT_PROFILE_DEGREES];
static atomic_uint_fast64_t cubic_batched_counts;
static atomic_uint_fast64_t cubic_batched_cpu_ns;
static atomic_uint_fast64_t cubic_batched_fallbacks;
static atomic_uint_fast64_t cubic_control_wall_ns;
static atomic_uint_fast64_t cubic_batched_wall_ns;

typedef struct {
    ulong lo, hi;
} big_key_component_t;

typedef struct {
    union {
        struct {
            ulong re, im;
        } word;
        struct {
            /*
             * The optimized FLINT-3.4 backend is deliberately restricted
             * to p < 2^128.  Fixed two-limb coordinates avoid two separately
             * allocated GMP integers per graph record.
             */
            big_key_component_t re, im;
        } big;
    } value;
} jkey_t;

static void frobenius_conjugate(jkey_t *result, const jkey_t *j);

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

/*
 * One logical ball participating in a shared rerandomization search.
 * All lanes at the same prime stage contribute their source points to one
 * multipoint batch, while graph ownership and parent chains remain separate.
 */
typedef struct {
    ordered_ball_t ball;
    index_vector_t current;
    index_vector_t next;
    search_graph_t *other;
    size_t *this_meet;
    size_t *other_meet;
    int *job_found;
    size_t job_slot;
    size_t batch_marker;
} ordered_shared_lane_t;

typedef struct {
    ordered_shared_lane_t *lane;
    size_t node_index;
} ordered_shared_source_t;

/* Zero means that a complete eligible prime-power layer is one batch. */
static size_t ordered_multipoint_batch_limit = 0;
/* Nonzero forces every nonempty ordered-search batch through a product tree. */
static int ordered_force_multipoint = 0;
/* Cooperatively split one large residual polynomial across factor levels. */
static int parallel_factor_tree_enabled = 1;
/*
 * Zero selects an automatic shard count bounded by the FLINT worker budget.
 * A nonzero value is an upper bound, not a request to create empty shards.
 */
static size_t evaluation_shard_limit = 0;

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
    const size_t *indices;
    size_t count;
    atomic_size_t next_index;
    size_t claim_size;
    int skip_degree_three;
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
    const word_evaluation_tasks_t *base;
    size_t shard_count;
    const size_t *offsets;
    const size_t *counts;
    fq_nmod_poly_struct ***trees;
    int initialize_values;
} word_sharded_evaluation_tasks_t;

typedef struct {
    const big_evaluation_tasks_t *base;
    size_t shard_count;
    const size_t *offsets;
    const size_t *counts;
    fq_poly_struct ***trees;
    int initialize_values;
} big_sharded_evaluation_tasks_t;

typedef struct {
    const fq_poly_struct *evaluated;
    const size_t *indices;
    size_t count;
    atomic_size_t next_index;
    size_t claim_size;
    int skip_degree_two;
    int skip_degree_three;
    size_t root_stride;
    jkey_t *root_keys;
    size_t *root_counts;
    const fq_ctx_struct *ctx;
} big_root_tasks_t;

typedef struct {
    const fq_poly_struct *evaluated;
    const size_t *indices;
    size_t count;
    atomic_size_t next_index;
    size_t claim_size;
    size_t root_stride;
    jkey_t *root_keys;
    size_t *root_counts;
    const fq_ctx_struct *ctx;
} quadratic_big_root_tasks_t;

typedef struct {
    const fq_nmod_poly_struct *evaluated;
    const size_t *indices;
    size_t count;
    atomic_size_t next_index;
    size_t lane_width;
    jkey_t *root_keys;
    size_t *root_counts;
    size_t root_stride;
    const fq_nmod_ctx_struct *ctx;
} cubic_root_tasks_t;

typedef struct {
    const fq_poly_struct *evaluated;
    const size_t *indices;
    size_t count;
    atomic_size_t next_index;
    size_t lane_width;
    jkey_t *root_keys;
    size_t *root_counts;
    size_t root_stride;
    const fq_ctx_struct *ctx;
} cubic_big_root_tasks_t;

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

static void *xcalloc(size_t count, size_t size)
{
    void *result;

    if (size != 0 && count > SIZE_MAX / size)
        die("allocation size overflow");
    result = calloc(count, size);
    if (result == NULL && count != 0)
        die("out of memory");
    return result;
}

static uint64_t thread_cpu_nanoseconds(void)
{
    struct timespec now;

    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &now) != 0)
        die("clock_gettime(CLOCK_THREAD_CPUTIME_ID) failed");
    return (uint64_t)now.tv_sec * UINT64_C(1000000000) +
           (uint64_t)now.tv_nsec;
}

static uint64_t elapsed_nanoseconds(double start, double end)
{
    double elapsed = end - start;

    if (elapsed <= 0.0)
        return 0;
    if (elapsed >= (double)UINT64_MAX / 1.0e9)
        return UINT64_MAX;
    return (uint64_t)(elapsed * 1.0e9 + 0.5);
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

static int big_component_equal(big_key_component_t x,
                               big_key_component_t y)
{
    return x.lo == y.lo && x.hi == y.hi;
}

static int big_component_compare(big_key_component_t x,
                                 big_key_component_t y)
{
    if (x.hi != y.hi)
        return x.hi < y.hi ? -1 : 1;
    if (x.lo != y.lo)
        return x.lo < y.lo ? -1 : 1;
    return 0;
}

static void big_component_from_fmpz(big_key_component_t *result,
                                    const fmpz_t source)
{
    if (fmpz_sgn(source) < 0 || fmpz_bits(source) > 128)
        die("internal field representative is outside [0, 2^128)");
    fmpz_get_uiui(&result->hi, &result->lo, source);
}

static void big_component_to_fmpz(fmpz_t result,
                                  big_key_component_t source)
{
    fmpz_set_uiui(result, source.hi, source.lo);
}

static int key_equal(jkey_t x, jkey_t y)
{
    if (active_backend == BACKEND_FQ_NMOD)
        return x.value.word.re == y.value.word.re &&
               x.value.word.im == y.value.word.im;
    return big_component_equal(x.value.big.re, y.value.big.re) &&
           big_component_equal(x.value.big.im, y.value.big.im);
}

static void key_init(jkey_t *key)
{
    if (active_backend == BACKEND_FQ_NMOD) {
        key->value.word.re = 0;
        key->value.word.im = 0;
    } else {
        key->value.big.re.lo = key->value.big.re.hi = 0;
        key->value.big.im.lo = key->value.big.im.hi = 0;
    }
}

static void key_clear(jkey_t *key)
{
    (void)key;
}

static void key_set(jkey_t *result, const jkey_t *source)
{
    *result = *source;
}

static void key_set_ui(jkey_t *key, ulong real, ulong imaginary)
{
    if (active_backend == BACKEND_FQ_NMOD) {
        key->value.word.re = real;
        key->value.word.im = imaginary;
    } else {
        key->value.big.re.lo = real;
        key->value.big.re.hi = 0;
        key->value.big.im.lo = imaginary;
        key->value.big.im.hi = 0;
    }
}

static void key_array_init(jkey_t *keys, size_t count)
{
    size_t i;

    if (active_backend == BACKEND_FQ)
        return;
    for (i = 0; i < count; i++)
        key_init(keys + i);
}

static void key_array_clear(jkey_t *keys, size_t count)
{
    size_t i;

    if (active_backend == BACKEND_FQ)
        return;
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
        uint64_t re_hash =
            mix64((uint64_t)key.value.big.re.lo) ^
            mix64((uint64_t)key.value.big.re.hi +
                  UINT64_C(0x9e3779b97f4a7c15));
        uint64_t im_hash =
            mix64((uint64_t)key.value.big.im.lo +
                  UINT64_C(0xd1b54a32d192ed03)) ^
            mix64((uint64_t)key.value.big.im.hi +
                  UINT64_C(0x94d049bb133111eb));

        return re_hash ^ ((im_hash << 1) | (im_hash >> 63));
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
    slong length;

    (void)representative;
    (void)ctx;
    fmpz_poly_fit_length(x, 2);
    big_component_to_fmpz(x->coeffs + 0, key.value.big.re);
    big_component_to_fmpz(x->coeffs + 1, key.value.big.im);
    if (key.value.big.im.lo != 0 || key.value.big.im.hi != 0)
        length = 2;
    else if (key.value.big.re.lo != 0 || key.value.big.re.hi != 0)
        length = 1;
    else
        length = 0;
    _fmpz_poly_set_length(x, length);
}

static void fq_to_key(jkey_t *key, const fq_t x,
                      fmpz_poly_t representative,
                      const fq_ctx_t ctx)
{
    (void)representative;
    (void)ctx;
    if (x->length > 2)
        die("noncanonical F_{p^2} representative");
    if (x->length > 0)
        big_component_from_fmpz(&key->value.big.re, x->coeffs + 0);
    else
        key->value.big.re.lo = key->value.big.re.hi = 0;
    if (x->length > 1)
        big_component_from_fmpz(&key->value.big.im, x->coeffs + 1);
    else
        key->value.big.im.lo = key->value.big.im.hi = 0;
}

static void fprint_key(FILE *stream, jkey_t key)
{
    if (active_backend == BACKEND_FQ_NMOD) {
        fprintf(stream, "(%lu, %lu)", (unsigned long)key.value.word.re,
                (unsigned long)key.value.word.im);
    } else {
        fmpz_t component;

        fmpz_init(component);
        fprintf(stream, "(");
        big_component_to_fmpz(component, key.value.big.re);
        fmpz_fprint(stream, component);
        fprintf(stream, ", ");
        big_component_to_fmpz(component, key.value.big.im);
        fmpz_fprint(stream, component);
        fprintf(stream, ")");
        fmpz_clear(component);
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
        int comparison = big_component_compare(
            a->value.big.re, b->value.big.re);

        if (comparison != 0)
            return comparison;
        return big_component_compare(
            a->value.big.im, b->value.big.im);
    }
}

/*
 * Keep many more queue claims than workers to absorb the random variation in
 * Rabin splitting, but amortize the cache-line transfer of next_index when a
 * batch is large.  The fixed upper bound prevents a late worker from owning a
 * long tail.  In particular, with 256 workers a large batch still exposes at
 * least 16384 independently claimable chunks.
 */
static size_t adaptive_root_claim_size(size_t count, size_t workers)
{
    size_t target_chunks, claim;

    if (workers <= 1)
        return count == 0 ? 1 : count;
    if (workers > SIZE_MAX / ROOT_QUEUE_CHUNKS_PER_WORKER)
        target_chunks = SIZE_MAX;
    else
        target_chunks = workers * ROOT_QUEUE_CHUNKS_PER_WORKER;
    claim = count / target_chunks;
    if (count % target_chunks != 0)
        claim++;
    if (claim == 0)
        claim = 1;
    if (claim > ROOT_QUEUE_MAX_CLAIM)
        claim = ROOT_QUEUE_MAX_CLAIM;
    return claim;
}

static size_t evaluation_shard_count(size_t point_count,
                                     unsigned n_threads,
                                     int use_fast)
{
    size_t worker_cap = (size_t)n_threads;
    size_t minimum_points =
        use_fast ? EVALUATION_MIN_POINTS_PER_SHARD
                 : HORNER_MIN_POINTS_PER_SHARD;
    size_t size_cap, shards;

    if (point_count == 0 ||
        n_threads < PARALLEL_LIFECYCLE_MIN_THREADS)
        return 1;
    if (evaluation_shard_limit != 0 &&
        worker_cap > evaluation_shard_limit)
        worker_cap = evaluation_shard_limit;
    if (worker_cap == 0)
        worker_cap = 1;

    /*
     * Keep every shard substantial enough to amortize one FLINT callback and,
     * for fast evaluation, one independent product tree. The quotient is
     * deliberately rounded down: a final short shard would otherwise impose
     * disproportionate setup overhead.
     */
    size_cap = point_count / minimum_points;
    if (size_cap == 0)
        size_cap = 1;
    shards = worker_cap < size_cap ? worker_cap : size_cap;
    if (shards > point_count)
        shards = point_count;
    return shards == 0 ? 1 : shards;
}

static int parallel_coefficient_lifecycle(size_t point_count,
                                          size_t coefficient_count,
                                          unsigned n_threads)
{
    if (n_threads < PARALLEL_LIFECYCLE_MIN_THREADS ||
        coefficient_count == 0 ||
        point_count > SIZE_MAX / coefficient_count)
        return 0;
    return point_count * coefficient_count >=
           PARALLEL_LIFECYCLE_MIN_VALUES;
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

static void build_evaluation_tree_task_word(slong shard_index,
                                            void *argument)
{
    word_sharded_evaluation_tasks_t *tasks = argument;
    size_t shard = (size_t)shard_index;
    size_t count = tasks->counts[shard];

    tasks->trees[shard] =
        _fq_nmod_poly_tree_alloc((slong)count, tasks->base->ctx);
    _fq_nmod_poly_tree_build(
        tasks->trees[shard],
        tasks->base->points + tasks->offsets[shard],
        (slong)count, tasks->base->ctx);
}

static void free_evaluation_tree_task_word(slong shard_index,
                                           void *argument)
{
    word_sharded_evaluation_tasks_t *tasks = argument;
    size_t shard = (size_t)shard_index;

    _fq_nmod_poly_tree_free(
        tasks->trees[shard], (slong)tasks->counts[shard],
        tasks->base->ctx);
}

static void evaluate_coefficient_shard_task_word(slong task_index,
                                                 void *argument)
{
    word_sharded_evaluation_tasks_t *tasks = argument;
    size_t task = (size_t)task_index;
    size_t y = task / tasks->shard_count;
    size_t shard = task % tasks->shard_count;
    size_t offset = tasks->offsets[shard];
    size_t count = tasks->counts[shard];
    const fq_nmod_poly_struct *coefficient =
        tasks->base->phi->coefficient_x.word + y;
    fq_nmod_struct *values =
        tasks->base->coefficient_values + y * tasks->base->count + offset;
    size_t i;

    if (tasks->initialize_values)
        for (i = 0; i < count; i++)
            fq_nmod_init(values + i, tasks->base->ctx);

    if (fq_nmod_poly_is_zero(coefficient, tasks->base->ctx))
        return;
    if (tasks->base->use_fast) {
        _fq_nmod_poly_evaluate_fq_nmod_vec_fast_precomp(
            values, coefficient->coeffs, coefficient->length,
            tasks->trees[shard], (slong)count, tasks->base->ctx);
    } else {
        _fq_nmod_poly_evaluate_fq_nmod_vec_iter(
            values, coefficient->coeffs, coefficient->length,
            tasks->base->points + offset, (slong)count,
            tasks->base->ctx);
    }
}

static void clear_coefficient_shard_task_word(slong task_index,
                                              void *argument)
{
    word_sharded_evaluation_tasks_t *tasks = argument;
    size_t task = (size_t)task_index;
    size_t y = task / tasks->shard_count;
    size_t shard = task % tasks->shard_count;
    size_t offset = tasks->offsets[shard];
    size_t count = tasks->counts[shard];
    fq_nmod_struct *values =
        tasks->base->coefficient_values + y * tasks->base->count + offset;
    size_t i;

    for (i = 0; i < count; i++)
        fq_nmod_clear(values + i, tasks->base->ctx);
}

/*
 * Evaluate a whole bivariate specialization batch without limiting
 * concurrency to the number of Y coefficients. Fast evaluation gets one
 * product tree per point shard; Horner evaluation uses the same two-
 * dimensional coefficient-by-shard task grid without product trees.
 */
static void evaluate_batch_parallel_word(
    const word_evaluation_tasks_t *base, size_t coefficient_count,
    unsigned n_threads)
{
    word_sharded_evaluation_tasks_t tasks;
    size_t *storage;
    size_t shard_count, shard, base_size, remainder, offset = 0;
    size_t task_count;

    if (base->count == 0 || coefficient_count == 0)
        return;
    shard_count =
        evaluation_shard_count(base->count, n_threads, base->use_fast);
    if (coefficient_count > SIZE_MAX / shard_count)
        die("evaluation task count overflow");
    task_count = coefficient_count * shard_count;
    if (task_count > (size_t)LONG_MAX)
        die("evaluation task count exceeds FLINT's slong interface");

    storage = xrealloc(NULL, 2 * shard_count, sizeof(*storage));
    tasks.base = base;
    tasks.shard_count = shard_count;
    tasks.offsets = storage;
    tasks.counts = storage + shard_count;
    tasks.trees = NULL;
    tasks.initialize_values = parallel_coefficient_lifecycle(
        base->count, coefficient_count, n_threads);
    if (!tasks.initialize_values) {
        size_t i;

        for (i = 0; i < base->count * coefficient_count; i++)
            fq_nmod_init(
                base->coefficient_values + i, base->ctx);
    }
    base_size = base->count / shard_count;
    remainder = base->count % shard_count;
    for (shard = 0; shard < shard_count; shard++) {
        size_t count = base_size + (shard < remainder ? 1 : 0);

        ((size_t *)tasks.offsets)[shard] = offset;
        ((size_t *)tasks.counts)[shard] = count;
        offset += count;
    }
    if (offset != base->count)
        die("internal error partitioning word evaluation points");

    if (base->use_fast) {
        tasks.trees =
            xrealloc(NULL, shard_count, sizeof(*tasks.trees));
        if (shard_count == 1)
            build_evaluation_tree_task_word(0, &tasks);
        else
            flint_parallel_do(
                build_evaluation_tree_task_word, &tasks,
                (slong)shard_count, (int)n_threads,
                FLINT_PARALLEL_UNIFORM);
    }
    flint_parallel_do(
        evaluate_coefficient_shard_task_word, &tasks,
        (slong)task_count, (int)n_threads,
        FLINT_PARALLEL_UNIFORM);
    if (base->use_fast) {
        if (shard_count == 1)
            free_evaluation_tree_task_word(0, &tasks);
        else
            flint_parallel_do(
                free_evaluation_tree_task_word, &tasks,
                (slong)shard_count, (int)n_threads,
                FLINT_PARALLEL_UNIFORM);
        free(tasks.trees);
    }
    free(storage);
}

static void clear_batch_parallel_word(
    const word_evaluation_tasks_t *base, size_t coefficient_count,
    unsigned n_threads)
{
    word_sharded_evaluation_tasks_t tasks;
    size_t *storage;
    size_t shard_count, shard, base_size, remainder, offset = 0;
    size_t task_count;

    if (base->count == 0 || coefficient_count == 0)
        return;
    if (!parallel_coefficient_lifecycle(
            base->count, coefficient_count, n_threads)) {
        size_t i;

        for (i = 0; i < base->count * coefficient_count; i++)
            fq_nmod_clear(
                base->coefficient_values + i, base->ctx);
        return;
    }
    shard_count =
        evaluation_shard_count(base->count, n_threads, base->use_fast);
    if (coefficient_count > SIZE_MAX / shard_count)
        die("coefficient-clear task count overflow");
    task_count = coefficient_count * shard_count;
    storage = xrealloc(NULL, 2 * shard_count, sizeof(*storage));
    tasks.base = base;
    tasks.shard_count = shard_count;
    tasks.offsets = storage;
    tasks.counts = storage + shard_count;
    tasks.trees = NULL;
    tasks.initialize_values = 0;
    base_size = base->count / shard_count;
    remainder = base->count % shard_count;
    for (shard = 0; shard < shard_count; shard++) {
        size_t count = base_size + (shard < remainder ? 1 : 0);

        ((size_t *)tasks.offsets)[shard] = offset;
        ((size_t *)tasks.counts)[shard] = count;
        offset += count;
    }
    flint_parallel_do(
        clear_coefficient_shard_task_word, &tasks,
        (slong)task_count, (int)n_threads,
        FLINT_PARALLEL_UNIFORM);
    free(storage);
}

/*
 * Append all roots of a monic polynomial of degree at most two.  inverse_two
 * is prepared once per worker.  Returning zero means that an unexpected
 * nonsplit quadratic was encountered, in which case the caller preserves the
 * old FLINT fallback behavior.
 */
static int append_small_roots_word(
    jkey_t *roots, size_t *count, size_t capacity,
    const fq_nmod_poly_t polynomial, fq_nmod_t b, fq_nmod_t c,
    fq_nmod_t discriminant, fq_nmod_t square_root, fq_nmod_t root,
    const fq_nmod_t inverse_two, nmod_poly_t representative,
    const fq_nmod_ctx_t ctx)
{
    slong degree = fq_nmod_poly_degree(polynomial, ctx);

    if (degree < 0)
        return 1;
    if (degree == 0)
        return 1;
    if (degree == 1) {
        if (*count == capacity)
            return 0;
        fq_nmod_poly_get_coeff(c, polynomial, 0, ctx);
        fq_nmod_neg(root, c, ctx);
        fq_nmod_to_key(roots + (*count)++, root, representative, ctx);
        return 1;
    }
    if (degree != 2)
        return 0;

    fq_nmod_poly_get_coeff(c, polynomial, 0, ctx);
    fq_nmod_poly_get_coeff(b, polynomial, 1, ctx);
    fq_nmod_mul(discriminant, b, b, ctx);
    fq_nmod_add(root, c, c, ctx);
    fq_nmod_add(root, root, root, ctx);
    fq_nmod_sub(discriminant, discriminant, root, ctx);
    if (!fq_nmod_sqrt(square_root, discriminant, ctx))
        return 0;

    fq_nmod_neg(b, b, ctx);
    fq_nmod_add(root, b, square_root, ctx);
    fq_nmod_mul(root, root, inverse_two, ctx);
    if (*count == capacity)
        return 0;
    fq_nmod_to_key(roots + (*count)++, root, representative, ctx);
    if (!fq_nmod_is_zero(square_root, ctx)) {
        fq_nmod_sub(root, b, square_root, ctx);
        fq_nmod_mul(root, root, inverse_two, ctx);
        if (*count == capacity)
            return 0;
        fq_nmod_to_key(roots + (*count)++, root, representative, ctx);
    }
    return 1;
}

/*
 * FLINT 3.4 made its internal Rabin splitter private.  Calling
 * fq_nmod_poly_roots for every specialization would repeatedly allocate a
 * random state, (q-1)/2, and FLINT_BITS scratch polynomials.  This compact
 * public-API implementation retains all of that state in the surrounding
 * worker and performs one proper Las Vegas split.
 */
static void split_rabin_public_word(
    fq_nmod_poly_t factor_a, fq_nmod_poly_t factor_b,
    const fq_nmod_poly_t polynomial, const fmpz_t halfq,
    fq_nmod_poly_t random_linear, fq_nmod_poly_t inverse_series,
    flint_rand_t randstate, const fq_nmod_ctx_t ctx)
{
    slong degree, split_degree;

    degree = fq_nmod_poly_degree(polynomial, ctx);
    if (degree <= 1)
        die("internal Rabin split received a linear polynomial");

    fq_nmod_poly_reverse(
        random_linear, polynomial, polynomial->length, ctx);
    fq_nmod_poly_inv_series_newton(
        inverse_series, random_linear, random_linear->length, ctx);

    for (;;) {
        fq_nmod_poly_fit_length(factor_a, 2, ctx);
        fq_nmod_rand(factor_a->coeffs + 0, randstate, ctx);
        fq_nmod_rand_not_zero(
            factor_a->coeffs + 1, randstate, ctx);
        factor_a->length = 2;

        fq_nmod_poly_powmod_fmpz_sliding_preinv(
            random_linear, factor_a, halfq, 0, polynomial,
            inverse_series, ctx);
        fq_nmod_poly_add_si(
            random_linear, random_linear, -1, ctx);
        fq_nmod_poly_gcd(
            factor_a, random_linear, polynomial, ctx);
        split_degree = fq_nmod_poly_degree(factor_a, ctx);
        if (split_degree > 0 && split_degree < degree)
            break;
    }

    fq_nmod_poly_div(factor_b, polynomial, factor_a, ctx);
    if (fq_nmod_poly_degree(factor_a, ctx) <
        fq_nmod_poly_degree(factor_b, ctx))
        fq_nmod_poly_swap(factor_a, factor_b, ctx);
}

/*
 * Split a polynomial completely using one persistent worker workspace.
 * _fq_nmod_poly_split_rabin returns its larger factor first.  Pushing that
 * factor before the smaller one makes the smaller one the next LIFO item, so
 * the stack depth is logarithmic and ROOT_FACTOR_STACK_SIZE is sufficient.
 * Every quadratic produced at any recursion level is solved directly.
 *
 * SIZE_MAX is an internal failure marker used only for the compatibility
 * fallback in find_roots_task_word.
 */
static size_t fully_split_roots_word(
    jkey_t *roots, size_t capacity, const fq_nmod_poly_t polynomial,
    const fmpz_t halfq, flint_rand_t randstate,
    fq_nmod_poly_t factor_a, fq_nmod_poly_t factor_b,
    fq_nmod_poly_t temporary, fq_nmod_poly_t inverse_series,
    fq_nmod_poly_struct *stack, fq_nmod_t b, fq_nmod_t c,
    fq_nmod_t discriminant, fq_nmod_t square_root, fq_nmod_t root,
    const fq_nmod_t inverse_two, nmod_poly_t representative,
    const fq_nmod_ctx_t ctx)
{
    size_t count = 0, stack_length = 1;

    fq_nmod_poly_make_monic(stack, polynomial, ctx);
    while (stack_length != 0) {
        fq_nmod_poly_struct *factor = stack + --stack_length;
        slong degree = fq_nmod_poly_degree(factor, ctx);

        if (degree <= 2) {
            if (!append_small_roots_word(
                    roots, &count, capacity, factor, b, c, discriminant,
                    square_root, root, inverse_two, representative, ctx))
                return SIZE_MAX;
            continue;
        }
        if (stack_length + 2 > ROOT_FACTOR_STACK_SIZE)
            return SIZE_MAX;

#if ISOGENY_USE_PRIVATE_RABIN
        _fq_nmod_poly_split_rabin(
            factor_a, factor_b, factor, halfq, temporary, inverse_series,
            randstate, ctx);
#else
        split_rabin_public_word(
            factor_a, factor_b, factor, halfq, temporary, inverse_series,
            randstate, ctx);
#endif

        /*
         * FLINT orders factor_a by nonincreasing degree. Swaps retain all
         * allocated coefficient storage for reuse in later specializations.
         */
        fq_nmod_poly_swap(stack + stack_length, factor_a, ctx);
        stack_length++;
        fq_nmod_poly_swap(stack + stack_length, factor_b, ctx);
        stack_length++;
    }
    return count;
}

typedef struct {
    fq_nmod_struct c0[CUBIC_BATCH_LANES];
    fq_nmod_struct c1[CUBIC_BATCH_LANES];
    fq_nmod_struct c2[CUBIC_BATCH_LANES];
    fq_nmod_struct r40[CUBIC_BATCH_LANES];
    fq_nmod_struct r41[CUBIC_BATCH_LANES];
    fq_nmod_struct r42[CUBIC_BATCH_LANES];
    fq_nmod_struct accumulator[3][CUBIC_BATCH_LANES];
    fq_nmod_struct base[3][CUBIC_BATCH_LANES];
    fq_nmod_struct power_scratch[3][CUBIC_BATCH_LANES];
    fq_nmod_struct numerator[CUBIC_BATCH_LANES];
    fq_nmod_struct denominator[CUBIC_BATCH_LANES];
    fq_nmod_struct prefix[CUBIC_BATCH_LANES];
    fq_nmod_struct affine_root[CUBIC_BATCH_LANES];
    fq_nmod_t t0, t1, t2, t3, t4, s0, s1, one;
    fq_nmod_t product, inverse_product, inverse_two;
    fq_nmod_t square_root, discriminant, q0, q1;
    fq_nmod_poly_factor_t fallback_factors;
    fmpz_t halfq;
    flint_rand_t randstate;
    nmod_poly_t representative;
} cubic_word_workspace_t;

static void cubic_word_workspace_init(
    cubic_word_workspace_t *workspace, slong task_index,
    const fq_nmod_ctx_t ctx)
{
    fq_nmod_struct *vectors[] = {
        workspace->c0, workspace->c1, workspace->c2,
        workspace->r40, workspace->r41, workspace->r42,
        workspace->accumulator[0], workspace->accumulator[1],
        workspace->accumulator[2], workspace->base[0],
        workspace->base[1], workspace->base[2],
        workspace->power_scratch[0], workspace->power_scratch[1],
        workspace->power_scratch[2], workspace->numerator,
        workspace->denominator, workspace->prefix,
        workspace->affine_root
    };
    fq_nmod_struct *scalars[] = {
        workspace->t0, workspace->t1, workspace->t2,
        workspace->t3, workspace->t4, workspace->s0,
        workspace->s1, workspace->one, workspace->product,
        workspace->inverse_product, workspace->inverse_two,
        workspace->square_root, workspace->discriminant,
        workspace->q0, workspace->q1
    };
    size_t vector_index, scalar_index, lane;
    uint64_t seed;

    for (vector_index = 0;
         vector_index < sizeof(vectors) / sizeof(vectors[0]);
         vector_index++)
        for (lane = 0; lane < CUBIC_BATCH_LANES; lane++)
            fq_nmod_init(vectors[vector_index] + lane, ctx);
    for (scalar_index = 0;
         scalar_index < sizeof(scalars) / sizeof(scalars[0]);
         scalar_index++)
        fq_nmod_init(scalars[scalar_index], ctx);
    fq_nmod_one(workspace->one, ctx);
    fq_nmod_set_ui(workspace->inverse_two, 2, ctx);
    fq_nmod_inv(workspace->inverse_two, workspace->inverse_two, ctx);
    fq_nmod_poly_factor_init(workspace->fallback_factors, ctx);
    fmpz_init(workspace->halfq);
    fq_nmod_ctx_order(workspace->halfq, ctx);
    fmpz_sub_ui(workspace->halfq, workspace->halfq, 1);
    fmpz_fdiv_q_2exp(workspace->halfq, workspace->halfq, 1);
    ISOGENY_RAND_INIT(workspace->randstate);
    seed = mix64(UINT64_C(0xa0761d6478bd642f) +
                 (uint64_t)(task_index + 1));
    ISOGENY_RAND_SET_SEED(
        workspace->randstate, (ulong)seed,
        (ulong)mix64(seed ^ UINT64_C(0xe7037ed1a0b428db)));
    nmod_poly_init(workspace->representative,
                   field_characteristic_word);
}

static void cubic_word_workspace_clear(
    cubic_word_workspace_t *workspace, const fq_nmod_ctx_t ctx)
{
    fq_nmod_struct *vectors[] = {
        workspace->c0, workspace->c1, workspace->c2,
        workspace->r40, workspace->r41, workspace->r42,
        workspace->accumulator[0], workspace->accumulator[1],
        workspace->accumulator[2], workspace->base[0],
        workspace->base[1], workspace->base[2],
        workspace->power_scratch[0], workspace->power_scratch[1],
        workspace->power_scratch[2], workspace->numerator,
        workspace->denominator, workspace->prefix,
        workspace->affine_root
    };
    fq_nmod_struct *scalars[] = {
        workspace->t0, workspace->t1, workspace->t2,
        workspace->t3, workspace->t4, workspace->s0,
        workspace->s1, workspace->one, workspace->product,
        workspace->inverse_product, workspace->inverse_two,
        workspace->square_root, workspace->discriminant,
        workspace->q0, workspace->q1
    };
    size_t vector_index, scalar_index, lane;

    nmod_poly_clear(workspace->representative);
    ISOGENY_RAND_CLEAR(workspace->randstate);
    fmpz_clear(workspace->halfq);
    fq_nmod_poly_factor_clear(workspace->fallback_factors, ctx);
    for (scalar_index = sizeof(scalars) / sizeof(scalars[0]);
         scalar_index-- > 0;)
        fq_nmod_clear(scalars[scalar_index], ctx);
    for (vector_index = sizeof(vectors) / sizeof(vectors[0]);
         vector_index-- > 0;)
        for (lane = CUBIC_BATCH_LANES; lane-- > 0;)
            fq_nmod_clear(vectors[vector_index] + lane, ctx);
}

/*
 * Multiply two residues modulo
 *   Y^3 + c2 Y^2 + c1 Y + c0.
 * The degree-three and degree-four reduction coefficients are
 *   (-c0,-c1,-c2) and
 *   (c2*c0,c2*c1-c0,c2^2-c1).
 */
static void cubic_residue_mul_word(
    fq_nmod_t out0, fq_nmod_t out1, fq_nmod_t out2,
    const fq_nmod_t a0, const fq_nmod_t a1, const fq_nmod_t a2,
    const fq_nmod_t b0, const fq_nmod_t b1, const fq_nmod_t b2,
    const fq_nmod_t c0, const fq_nmod_t c1, const fq_nmod_t c2,
    const fq_nmod_t r40, const fq_nmod_t r41,
    const fq_nmod_t r42, cubic_word_workspace_t *workspace,
    const fq_nmod_ctx_t ctx)
{
    fq_nmod_mul(workspace->t0, a0, b0, ctx);
    fq_nmod_mul(workspace->t1, a0, b1, ctx);
    fq_nmod_mul(workspace->s0, a1, b0, ctx);
    fq_nmod_add(workspace->t1, workspace->t1, workspace->s0, ctx);
    fq_nmod_mul(workspace->t2, a0, b2, ctx);
    fq_nmod_mul(workspace->s0, a1, b1, ctx);
    fq_nmod_add(workspace->t2, workspace->t2, workspace->s0, ctx);
    fq_nmod_mul(workspace->s0, a2, b0, ctx);
    fq_nmod_add(workspace->t2, workspace->t2, workspace->s0, ctx);
    fq_nmod_mul(workspace->t3, a1, b2, ctx);
    fq_nmod_mul(workspace->s0, a2, b1, ctx);
    fq_nmod_add(workspace->t3, workspace->t3, workspace->s0, ctx);
    fq_nmod_mul(workspace->t4, a2, b2, ctx);

    fq_nmod_mul(workspace->s0, c0, workspace->t3, ctx);
    fq_nmod_sub(out0, workspace->t0, workspace->s0, ctx);
    fq_nmod_mul(workspace->s0, r40, workspace->t4, ctx);
    fq_nmod_add(out0, out0, workspace->s0, ctx);

    fq_nmod_mul(workspace->s0, c1, workspace->t3, ctx);
    fq_nmod_sub(out1, workspace->t1, workspace->s0, ctx);
    fq_nmod_mul(workspace->s0, r41, workspace->t4, ctx);
    fq_nmod_add(out1, out1, workspace->s0, ctx);

    fq_nmod_mul(workspace->s0, c2, workspace->t3, ctx);
    fq_nmod_sub(out2, workspace->t2, workspace->s0, ctx);
    fq_nmod_mul(workspace->s0, r42, workspace->t4, ctx);
    fq_nmod_add(out2, out2, workspace->s0, ctx);
}

static void cubic_residue_square_word(
    fq_nmod_t out0, fq_nmod_t out1, fq_nmod_t out2,
    const fq_nmod_t a0, const fq_nmod_t a1, const fq_nmod_t a2,
    const fq_nmod_t c0, const fq_nmod_t c1, const fq_nmod_t c2,
    const fq_nmod_t r40, const fq_nmod_t r41,
    const fq_nmod_t r42, cubic_word_workspace_t *workspace,
    const fq_nmod_ctx_t ctx)
{
    fq_nmod_mul(workspace->t0, a0, a0, ctx);
    fq_nmod_mul(workspace->t1, a0, a1, ctx);
    fq_nmod_add(workspace->t1, workspace->t1, workspace->t1, ctx);
    fq_nmod_mul(workspace->t2, a0, a2, ctx);
    fq_nmod_add(workspace->t2, workspace->t2, workspace->t2, ctx);
    fq_nmod_mul(workspace->s0, a1, a1, ctx);
    fq_nmod_add(workspace->t2, workspace->t2, workspace->s0, ctx);
    fq_nmod_mul(workspace->t3, a1, a2, ctx);
    fq_nmod_add(workspace->t3, workspace->t3, workspace->t3, ctx);
    fq_nmod_mul(workspace->t4, a2, a2, ctx);

    fq_nmod_mul(workspace->s0, c0, workspace->t3, ctx);
    fq_nmod_sub(out0, workspace->t0, workspace->s0, ctx);
    fq_nmod_mul(workspace->s0, r40, workspace->t4, ctx);
    fq_nmod_add(out0, out0, workspace->s0, ctx);

    fq_nmod_mul(workspace->s0, c1, workspace->t3, ctx);
    fq_nmod_sub(out1, workspace->t1, workspace->s0, ctx);
    fq_nmod_mul(workspace->s0, r41, workspace->t4, ctx);
    fq_nmod_add(out1, out1, workspace->s0, ctx);

    fq_nmod_mul(workspace->s0, c2, workspace->t3, ctx);
    fq_nmod_sub(out2, workspace->t2, workspace->s0, ctx);
    fq_nmod_mul(workspace->s0, r42, workspace->t4, ctx);
    fq_nmod_add(out2, out2, workspace->s0, ctx);
}

/*
 * Find one root of f projectively from h = g +/- 1.  For quadratic h,
 * the formulas below are the pseudo-remainder of f by h, scaled by lc(h)^2:
 *
 * r1 = A^2*c1 - A*C - A*c2*B + B^2,
 * r0 = A^2*c0 - A*c2*C + B*C.
 *
 * This avoids every inversion in the split step.
 */
static int cubic_projective_root_word(
    fq_nmod_t numerator, fq_nmod_t denominator,
    const fq_nmod_t c0, const fq_nmod_t c1, const fq_nmod_t c2,
    const fq_nmod_t h0, const fq_nmod_t h1, const fq_nmod_t h2,
    cubic_word_workspace_t *workspace, const fq_nmod_ctx_t ctx)
{
    if (fq_nmod_is_zero(h2, ctx)) {
        if (fq_nmod_is_zero(h1, ctx))
            return 0;
        fq_nmod_neg(numerator, h0, ctx);
        fq_nmod_set(denominator, h1, ctx);
        return 1;
    }

    /* h0 is allowed to alias workspace->t0 at the call site. */
    fq_nmod_set(workspace->s1, h0, ctx);
    fq_nmod_mul(workspace->t0, h2, h2, ctx);

    fq_nmod_mul(workspace->t1, workspace->t0, c1, ctx);
    fq_nmod_mul(workspace->s0, h2, workspace->s1, ctx);
    fq_nmod_sub(workspace->t1, workspace->t1, workspace->s0, ctx);
    fq_nmod_mul(workspace->s0, h2, c2, ctx);
    fq_nmod_mul(workspace->s0, workspace->s0, h1, ctx);
    fq_nmod_sub(workspace->t1, workspace->t1, workspace->s0, ctx);
    fq_nmod_mul(workspace->s0, h1, h1, ctx);
    fq_nmod_add(workspace->t1, workspace->t1, workspace->s0, ctx);

    fq_nmod_mul(workspace->t2, workspace->t0, c0, ctx);
    fq_nmod_mul(workspace->s0, h2, c2, ctx);
    fq_nmod_mul(workspace->s0, workspace->s0, workspace->s1, ctx);
    fq_nmod_sub(workspace->t2, workspace->t2, workspace->s0, ctx);
    fq_nmod_mul(workspace->s0, h1, workspace->s1, ctx);
    fq_nmod_add(workspace->t2, workspace->t2, workspace->s0, ctx);

    if (!fq_nmod_is_zero(workspace->t1, ctx)) {
        fq_nmod_neg(numerator, workspace->t2, ctx);
        fq_nmod_set(denominator, workspace->t1, ctx);
        return 1;
    }
    if (!fq_nmod_is_zero(workspace->t2, ctx))
        return 0;

    /* h divides f; the complementary linear factor has root B/A-c2. */
    fq_nmod_mul(workspace->s0, h2, c2, ctx);
    fq_nmod_sub(numerator, h1, workspace->s0, ctx);
    fq_nmod_set(denominator, h2, ctx);
    return 1;
}

static int cubic_projective_root_is_valid_word(
    const fq_nmod_t numerator, const fq_nmod_t denominator,
    const fq_nmod_t c0, const fq_nmod_t c1, const fq_nmod_t c2,
    cubic_word_workspace_t *workspace, const fq_nmod_ctx_t ctx)
{
    /* n^3 + c2*n^2*d + c1*n*d^2 + c0*d^3 */
    fq_nmod_mul(workspace->t0, numerator, numerator, ctx);
    fq_nmod_mul(workspace->t1, workspace->t0, numerator, ctx);
    fq_nmod_mul(workspace->s0, c2, workspace->t0, ctx);
    fq_nmod_mul(workspace->s0, workspace->s0, denominator, ctx);
    fq_nmod_add(workspace->t1, workspace->t1, workspace->s0, ctx);
    fq_nmod_mul(workspace->t2, denominator, denominator, ctx);
    fq_nmod_mul(workspace->s0, c1, numerator, ctx);
    fq_nmod_mul(workspace->s0, workspace->s0, workspace->t2, ctx);
    fq_nmod_add(workspace->t1, workspace->t1, workspace->s0, ctx);
    fq_nmod_mul(workspace->t2, workspace->t2, denominator, ctx);
    fq_nmod_mul(workspace->s0, c0, workspace->t2, ctx);
    fq_nmod_add(workspace->t1, workspace->t1, workspace->s0, ctx);
    return fq_nmod_is_zero(workspace->t1, ctx);
}

static size_t sort_unique_root_keys(jkey_t *roots, size_t count)
{
    size_t read_index, unique_count;

    qsort(roots, count, sizeof(*roots), compare_keys);
    if (count <= 1)
        return count;
    unique_count = 1;
    for (read_index = 1; read_index < count; read_index++) {
        if (!key_equal(roots[read_index], roots[unique_count - 1])) {
            if (unique_count != read_index)
                key_set(roots + unique_count, roots + read_index);
            unique_count++;
        }
    }
    return unique_count;
}

static int finish_cubic_roots_word(
    jkey_t *roots, size_t *count_out, size_t capacity,
    const fq_nmod_t root, const fq_nmod_t c0,
    const fq_nmod_t c1, const fq_nmod_t c2,
    cubic_word_workspace_t *workspace, const fq_nmod_ctx_t ctx)
{
    size_t count = 0;

    if (capacity < 3)
        return 0;
    /* Check f(root) before using synthetic division. */
    fq_nmod_mul(workspace->t0, root, root, ctx);
    fq_nmod_mul(workspace->t1, workspace->t0, root, ctx);
    fq_nmod_mul(workspace->s0, c2, workspace->t0, ctx);
    fq_nmod_add(workspace->t1, workspace->t1, workspace->s0, ctx);
    fq_nmod_mul(workspace->s0, c1, root, ctx);
    fq_nmod_add(workspace->t1, workspace->t1, workspace->s0, ctx);
    fq_nmod_add(workspace->t1, workspace->t1, c0, ctx);
    if (!fq_nmod_is_zero(workspace->t1, ctx))
        return 0;

    fq_nmod_to_key(roots + count++, root, workspace->representative, ctx);
    fq_nmod_add(workspace->q1, c2, root, ctx);
    fq_nmod_mul(workspace->q0, root, workspace->q1, ctx);
    fq_nmod_add(workspace->q0, workspace->q0, c1, ctx);
    fq_nmod_mul(workspace->discriminant,
                workspace->q1, workspace->q1, ctx);
    fq_nmod_add(workspace->t0, workspace->q0, workspace->q0, ctx);
    fq_nmod_add(workspace->t0, workspace->t0, workspace->t0, ctx);
    fq_nmod_sub(workspace->discriminant,
                workspace->discriminant, workspace->t0, ctx);
    if (!fq_nmod_sqrt(
            workspace->square_root, workspace->discriminant, ctx))
        return 0;

    fq_nmod_neg(workspace->q1, workspace->q1, ctx);
    fq_nmod_add(workspace->t0,
                workspace->q1, workspace->square_root, ctx);
    fq_nmod_mul(workspace->t0,
                workspace->t0, workspace->inverse_two, ctx);
    fq_nmod_to_key(
        roots + count++, workspace->t0, workspace->representative, ctx);
    if (!fq_nmod_is_zero(workspace->square_root, ctx)) {
        fq_nmod_sub(workspace->t0,
                    workspace->q1, workspace->square_root, ctx);
        fq_nmod_mul(workspace->t0,
                    workspace->t0, workspace->inverse_two, ctx);
        fq_nmod_to_key(
            roots + count++, workspace->t0,
            workspace->representative, ctx);
    }
    *count_out = sort_unique_root_keys(roots, count);
    return 1;
}

static size_t fallback_cubic_roots_word(
    jkey_t *roots, size_t capacity, const fq_nmod_poly_t polynomial,
    cubic_word_workspace_t *workspace, const fq_nmod_ctx_t ctx)
{
    size_t count = 0;
    slong factor_index;

    fq_nmod_poly_roots(
        workspace->fallback_factors, polynomial, 0, ctx);
    for (factor_index = 0;
         factor_index < workspace->fallback_factors->num;
         factor_index++) {
        if (fq_nmod_poly_degree(
                workspace->fallback_factors->poly + factor_index,
                ctx) != 1)
            continue;
        if (count == capacity)
            break;
        fq_nmod_poly_get_coeff(
            workspace->t0,
            workspace->fallback_factors->poly + factor_index, 0, ctx);
        fq_nmod_neg(workspace->t0, workspace->t0, ctx);
        fq_nmod_to_key(
            roots + count++, workspace->t0,
            workspace->representative, ctx);
    }
    return sort_unique_root_keys(roots, count);
}

static void find_cubic_roots_task_word(slong task_index, void *argument)
{
    cubic_root_tasks_t *tasks = argument;
    cubic_word_workspace_t workspace;
    size_t first, lane_count, lane, attempt;
    size_t actual_indices[CUBIC_BATCH_LANES];
    unsigned char resolved[CUBIC_BATCH_LANES];
    unsigned char usable[CUBIC_BATCH_LANES];
    uint64_t local_count = 0, local_cpu_ns = 0, local_fallbacks = 0;

    cubic_word_workspace_init(&workspace, task_index, tasks->ctx);
    for (;;) {
        uint64_t group_start;
        slong bit;

        first = atomic_fetch_add_explicit(
            &tasks->next_index, tasks->lane_width,
            memory_order_relaxed);
        if (first >= tasks->count)
            break;
        lane_count = tasks->count - first;
        if (lane_count > tasks->lane_width)
            lane_count = tasks->lane_width;
        group_start = thread_cpu_nanoseconds();

        for (lane = 0; lane < lane_count; lane++) {
            size_t actual_index = tasks->indices[first + lane];

            actual_indices[lane] = actual_index;
            resolved[lane] = 0;
            usable[lane] = 0;
            fq_nmod_poly_get_coeff(
                workspace.c0 + lane, tasks->evaluated + actual_index,
                0, tasks->ctx);
            fq_nmod_poly_get_coeff(
                workspace.c1 + lane, tasks->evaluated + actual_index,
                1, tasks->ctx);
            fq_nmod_poly_get_coeff(
                workspace.c2 + lane, tasks->evaluated + actual_index,
                2, tasks->ctx);
            fq_nmod_mul(workspace.r40 + lane,
                        workspace.c2 + lane, workspace.c0 + lane,
                        tasks->ctx);
            fq_nmod_mul(workspace.r41 + lane,
                        workspace.c2 + lane, workspace.c1 + lane,
                        tasks->ctx);
            fq_nmod_sub(workspace.r41 + lane,
                        workspace.r41 + lane, workspace.c0 + lane,
                        tasks->ctx);
            fq_nmod_mul(workspace.r42 + lane,
                        workspace.c2 + lane, workspace.c2 + lane,
                        tasks->ctx);
            fq_nmod_sub(workspace.r42 + lane,
                        workspace.r42 + lane, workspace.c1 + lane,
                        tasks->ctx);
        }

        for (attempt = 0; attempt < CUBIC_MAX_SPLIT_ATTEMPTS;
             attempt++) {
            size_t unresolved = 0;

            for (lane = 0; lane < lane_count; lane++) {
                if (resolved[lane])
                    continue;
                unresolved++;
                fq_nmod_rand(
                    workspace.base[0] + lane,
                    workspace.randstate, tasks->ctx);
                fq_nmod_rand_not_zero(
                    workspace.base[1] + lane,
                    workspace.randstate, tasks->ctx);
                fq_nmod_zero(workspace.base[2] + lane, tasks->ctx);
                fq_nmod_one(
                    workspace.accumulator[0] + lane, tasks->ctx);
                fq_nmod_zero(
                    workspace.accumulator[1] + lane, tasks->ctx);
                fq_nmod_zero(
                    workspace.accumulator[2] + lane, tasks->ctx);
            }
            if (unresolved == 0)
                break;

            for (bit = (slong)fmpz_bits(workspace.halfq) - 1;
                 bit >= 0; bit--) {
                for (lane = 0; lane < lane_count; lane++) {
                    if (resolved[lane])
                        continue;
                    cubic_residue_square_word(
                        workspace.power_scratch[0] + lane,
                        workspace.power_scratch[1] + lane,
                        workspace.power_scratch[2] + lane,
                        workspace.accumulator[0] + lane,
                        workspace.accumulator[1] + lane,
                        workspace.accumulator[2] + lane,
                        workspace.c0 + lane, workspace.c1 + lane,
                        workspace.c2 + lane, workspace.r40 + lane,
                        workspace.r41 + lane, workspace.r42 + lane,
                        &workspace, tasks->ctx);
                    fq_nmod_swap(
                        workspace.accumulator[0] + lane,
                        workspace.power_scratch[0] + lane, tasks->ctx);
                    fq_nmod_swap(
                        workspace.accumulator[1] + lane,
                        workspace.power_scratch[1] + lane, tasks->ctx);
                    fq_nmod_swap(
                        workspace.accumulator[2] + lane,
                        workspace.power_scratch[2] + lane, tasks->ctx);
                }
                if (!fmpz_tstbit(workspace.halfq, (ulong)bit))
                    continue;
                for (lane = 0; lane < lane_count; lane++) {
                    if (resolved[lane])
                        continue;
                    cubic_residue_mul_word(
                        workspace.power_scratch[0] + lane,
                        workspace.power_scratch[1] + lane,
                        workspace.power_scratch[2] + lane,
                        workspace.accumulator[0] + lane,
                        workspace.accumulator[1] + lane,
                        workspace.accumulator[2] + lane,
                        workspace.base[0] + lane,
                        workspace.base[1] + lane,
                        workspace.base[2] + lane,
                        workspace.c0 + lane, workspace.c1 + lane,
                        workspace.c2 + lane, workspace.r40 + lane,
                        workspace.r41 + lane, workspace.r42 + lane,
                        &workspace, tasks->ctx);
                    fq_nmod_swap(
                        workspace.accumulator[0] + lane,
                        workspace.power_scratch[0] + lane, tasks->ctx);
                    fq_nmod_swap(
                        workspace.accumulator[1] + lane,
                        workspace.power_scratch[1] + lane, tasks->ctx);
                    fq_nmod_swap(
                        workspace.accumulator[2] + lane,
                        workspace.power_scratch[2] + lane, tasks->ctx);
                }
            }

            for (lane = 0; lane < lane_count; lane++) {
                int found;

                if (resolved[lane])
                    continue;
                fq_nmod_sub_one(
                    workspace.t0, workspace.accumulator[0] + lane,
                    tasks->ctx);
                found = cubic_projective_root_word(
                    workspace.numerator + lane,
                    workspace.denominator + lane,
                    workspace.c0 + lane, workspace.c1 + lane,
                    workspace.c2 + lane, workspace.t0,
                    workspace.accumulator[1] + lane,
                    workspace.accumulator[2] + lane,
                    &workspace, tasks->ctx);
                if (!found) {
                    fq_nmod_add(
                        workspace.t0, workspace.accumulator[0] + lane,
                        workspace.one, tasks->ctx);
                    found = cubic_projective_root_word(
                        workspace.numerator + lane,
                        workspace.denominator + lane,
                        workspace.c0 + lane, workspace.c1 + lane,
                        workspace.c2 + lane, workspace.t0,
                        workspace.accumulator[1] + lane,
                        workspace.accumulator[2] + lane,
                        &workspace, tasks->ctx);
                }
                if (found &&
                    cubic_projective_root_is_valid_word(
                        workspace.numerator + lane,
                        workspace.denominator + lane,
                        workspace.c0 + lane, workspace.c1 + lane,
                        workspace.c2 + lane, &workspace, tasks->ctx))
                    resolved[lane] = usable[lane] = 1;
            }
        }

        /* Montgomery batch inversion of all successful denominators. */
        fq_nmod_one(workspace.product, tasks->ctx);
        for (lane = 0; lane < lane_count; lane++) {
            if (!usable[lane])
                continue;
            fq_nmod_set(
                workspace.prefix + lane, workspace.product, tasks->ctx);
            fq_nmod_mul(
                workspace.product, workspace.product,
                workspace.denominator + lane, tasks->ctx);
        }
        fq_nmod_inv(
            workspace.inverse_product, workspace.product, tasks->ctx);
        for (lane = lane_count; lane-- > 0;) {
            if (!usable[lane])
                continue;
            fq_nmod_mul(
                workspace.t0, workspace.inverse_product,
                workspace.prefix + lane, tasks->ctx);
            fq_nmod_mul(
                workspace.affine_root + lane,
                workspace.numerator + lane, workspace.t0, tasks->ctx);
            fq_nmod_mul(
                workspace.inverse_product, workspace.inverse_product,
                workspace.denominator + lane, tasks->ctx);
        }

        for (lane = 0; lane < lane_count; lane++) {
            size_t actual_index = actual_indices[lane];
            jkey_t *point_roots =
                tasks->root_keys + actual_index * tasks->root_stride;
            size_t count = 0;

            if (!usable[lane] ||
                !finish_cubic_roots_word(
                    point_roots, &count, tasks->root_stride,
                    workspace.affine_root + lane,
                    workspace.c0 + lane, workspace.c1 + lane,
                    workspace.c2 + lane, &workspace, tasks->ctx)) {
                count = fallback_cubic_roots_word(
                    point_roots, tasks->root_stride,
                    tasks->evaluated + actual_index,
                    &workspace, tasks->ctx);
                local_fallbacks++;
            }
            tasks->root_counts[actual_index] = count;
        }
        local_count += lane_count;
        local_cpu_ns += thread_cpu_nanoseconds() - group_start;
    }

    atomic_fetch_add_explicit(
        &cubic_batched_counts, local_count, memory_order_relaxed);
    atomic_fetch_add_explicit(
        &cubic_batched_cpu_ns, local_cpu_ns, memory_order_relaxed);
    atomic_fetch_add_explicit(
        &cubic_batched_fallbacks, local_fallbacks, memory_order_relaxed);
    if (root_profile_enabled &&
        cubic_kernel == CUBIC_KERNEL_BATCHED) {
        atomic_fetch_add_explicit(
            root_profile_counts + 3, local_count,
            memory_order_relaxed);
        atomic_fetch_add_explicit(
            root_profile_cpu_ns + 3, local_cpu_ns,
            memory_order_relaxed);
    }
    cubic_word_workspace_clear(&workspace, tasks->ctx);
}

typedef struct {
    fq_struct c0[CUBIC_BATCH_LANES];
    fq_struct c1[CUBIC_BATCH_LANES];
    fq_struct c2[CUBIC_BATCH_LANES];
    fq_struct r40[CUBIC_BATCH_LANES];
    fq_struct r41[CUBIC_BATCH_LANES];
    fq_struct r42[CUBIC_BATCH_LANES];
    fq_struct accumulator[3][CUBIC_BATCH_LANES];
    fq_struct base[3][CUBIC_BATCH_LANES];
    fq_struct power_scratch[3][CUBIC_BATCH_LANES];
    fq_struct numerator[CUBIC_BATCH_LANES];
    fq_struct denominator[CUBIC_BATCH_LANES];
    fq_struct prefix[CUBIC_BATCH_LANES];
    fq_struct affine_root[CUBIC_BATCH_LANES];
    fq_t t0, t1, t2, t3, t4, s0, s1, one;
    fq_t product, inverse_product, inverse_two;
    fq_t square_root, discriminant, q0, q1;
    fq_poly_factor_t fallback_factors;
    fmpz_t halfq;
    fast_fp2_workspace_t fast;
    flint_rand_t randstate;
    fmpz_poly_t representative;
} cubic_big_workspace_t;

FLINT_STATIC_NOINLINE void cubic_big_workspace_init(
    cubic_big_workspace_t *workspace, slong task_index,
    const fq_ctx_t ctx)
{
    fq_struct *vectors[] = {
        workspace->c0, workspace->c1, workspace->c2,
        workspace->r40, workspace->r41, workspace->r42,
        workspace->accumulator[0], workspace->accumulator[1],
        workspace->accumulator[2], workspace->base[0],
        workspace->base[1], workspace->base[2],
        workspace->power_scratch[0], workspace->power_scratch[1],
        workspace->power_scratch[2], workspace->numerator,
        workspace->denominator, workspace->prefix,
        workspace->affine_root
    };
    fq_struct *scalars[] = {
        workspace->t0, workspace->t1, workspace->t2,
        workspace->t3, workspace->t4, workspace->s0,
        workspace->s1, workspace->one, workspace->product,
        workspace->inverse_product, workspace->inverse_two,
        workspace->square_root, workspace->discriminant,
        workspace->q0, workspace->q1
    };
    size_t vector_index, scalar_index, lane;
    uint64_t seed;

    for (vector_index = 0;
         vector_index < sizeof(vectors) / sizeof(vectors[0]);
         vector_index++)
        for (lane = 0; lane < CUBIC_BATCH_LANES; lane++)
            fq_init(vectors[vector_index] + lane, ctx);
    for (scalar_index = 0;
         scalar_index < sizeof(scalars) / sizeof(scalars[0]);
         scalar_index++)
        fq_init(scalars[scalar_index], ctx);
    fast_fp2_workspace_init(&workspace->fast, ctx);
    fq_one(workspace->one, ctx);
    fq_set_ui(workspace->inverse_two, 2, ctx);
    fq_inv(workspace->inverse_two, workspace->inverse_two, ctx);
    fq_poly_factor_init(workspace->fallback_factors, ctx);
    fmpz_init(workspace->halfq);
    fq_ctx_order(workspace->halfq, ctx);
    fmpz_sub_ui(workspace->halfq, workspace->halfq, 1);
    fmpz_fdiv_q_2exp(workspace->halfq, workspace->halfq, 1);
    ISOGENY_RAND_INIT(workspace->randstate);
    seed = mix64(UINT64_C(0xa0761d6478bd642f) +
                 (uint64_t)(task_index + 1));
    ISOGENY_RAND_SET_SEED(
        workspace->randstate, (ulong)seed,
        (ulong)mix64(seed ^ UINT64_C(0xe7037ed1a0b428db)));
    fmpz_poly_init(workspace->representative);
}

static void cubic_big_workspace_clear(
    cubic_big_workspace_t *workspace, const fq_ctx_t ctx)
{
    fq_struct *vectors[] = {
        workspace->c0, workspace->c1, workspace->c2,
        workspace->r40, workspace->r41, workspace->r42,
        workspace->accumulator[0], workspace->accumulator[1],
        workspace->accumulator[2], workspace->base[0],
        workspace->base[1], workspace->base[2],
        workspace->power_scratch[0], workspace->power_scratch[1],
        workspace->power_scratch[2], workspace->numerator,
        workspace->denominator, workspace->prefix,
        workspace->affine_root
    };
    fq_struct *scalars[] = {
        workspace->t0, workspace->t1, workspace->t2,
        workspace->t3, workspace->t4, workspace->s0,
        workspace->s1, workspace->one, workspace->product,
        workspace->inverse_product, workspace->inverse_two,
        workspace->square_root, workspace->discriminant,
        workspace->q0, workspace->q1
    };
    size_t vector_index, scalar_index, lane;

    fmpz_poly_clear(workspace->representative);
    ISOGENY_RAND_CLEAR(workspace->randstate);
    fast_fp2_workspace_clear(&workspace->fast);
    fmpz_clear(workspace->halfq);
    fq_poly_factor_clear(workspace->fallback_factors, ctx);
    for (scalar_index = sizeof(scalars) / sizeof(scalars[0]);
         scalar_index-- > 0;)
        fq_clear(scalars[scalar_index], ctx);
    for (vector_index = sizeof(vectors) / sizeof(vectors[0]);
         vector_index-- > 0;)
        for (lane = CUBIC_BATCH_LANES; lane-- > 0;)
            fq_clear(vectors[vector_index] + lane, ctx);
}

static const fmpz *cubic_big_component(
    const fq_t value, slong component,
    const fast_fp2_workspace_t *workspace)
{
    if (component < value->length)
        return value->coeffs + component;
    return workspace->zero;
}

static void cubic_big_normalise(fq_t value)
{
    _fmpz_poly_set_length(value, 2);
    _fmpz_poly_normalise(value);
}

static void cubic_big_add(
    fq_t result, const fq_t left, const fq_t right,
    fast_fp2_workspace_t *workspace, const fq_ctx_t ctx)
{
    const fmpz *left_re, *left_im, *right_re, *right_im;

    if (!quadratic_basis_is_minus_one) {
        fq_add(result, left, right, ctx);
        return;
    }
    fmpz_poly_fit_length(result, 2);
    left_re = cubic_big_component(left, 0, workspace);
    left_im = cubic_big_component(left, 1, workspace);
    right_re = cubic_big_component(right, 0, workspace);
    right_im = cubic_big_component(right, 1, workspace);
    fmpz_mod_add(result->coeffs + 0, left_re, right_re, ctx->ctxp);
    fmpz_mod_add(result->coeffs + 1, left_im, right_im, ctx->ctxp);
    cubic_big_normalise(result);
}

static void cubic_big_sub(
    fq_t result, const fq_t left, const fq_t right,
    fast_fp2_workspace_t *workspace, const fq_ctx_t ctx)
{
    const fmpz *left_re, *left_im, *right_re, *right_im;

    if (!quadratic_basis_is_minus_one) {
        fq_sub(result, left, right, ctx);
        return;
    }
    fmpz_poly_fit_length(result, 2);
    left_re = cubic_big_component(left, 0, workspace);
    left_im = cubic_big_component(left, 1, workspace);
    right_re = cubic_big_component(right, 0, workspace);
    right_im = cubic_big_component(right, 1, workspace);
    fmpz_mod_sub(result->coeffs + 0, left_re, right_re, ctx->ctxp);
    fmpz_mod_sub(result->coeffs + 1, left_im, right_im, ctx->ctxp);
    cubic_big_normalise(result);
}

static void cubic_big_neg(
    fq_t result, const fq_t input,
    fast_fp2_workspace_t *workspace, const fq_ctx_t ctx)
{
    const fmpz *input_re, *input_im;

    if (!quadratic_basis_is_minus_one) {
        fq_neg(result, input, ctx);
        return;
    }
    fmpz_poly_fit_length(result, 2);
    input_re = cubic_big_component(input, 0, workspace);
    input_im = cubic_big_component(input, 1, workspace);
    fmpz_mod_neg(result->coeffs + 0, input_re, ctx->ctxp);
    fmpz_mod_neg(result->coeffs + 1, input_im, ctx->ctxp);
    cubic_big_normalise(result);
}

static void cubic_big_mul(
    fq_t result, const fq_t left, const fq_t right,
    fast_fp2_workspace_t *workspace, const fq_ctx_t ctx)
{
    const fmpz *left_re, *left_im, *right_re, *right_im;

    if (!quadratic_basis_is_minus_one) {
        fq_mul(result, left, right, ctx);
        return;
    }
    fmpz_poly_fit_length(result, 2);
    left_re = cubic_big_component(left, 0, workspace);
    left_im = cubic_big_component(left, 1, workspace);
    right_re = cubic_big_component(right, 0, workspace);
    right_im = cubic_big_component(right, 1, workspace);

    fmpz_mod_mul(workspace->scratch[0],
                 left_re, right_re, ctx->ctxp);
    fmpz_mod_mul(workspace->scratch[1],
                 left_im, right_im, ctx->ctxp);
    fmpz_mod_mul(workspace->scratch[2],
                 left_re, right_im, ctx->ctxp);
    fmpz_mod_mul(workspace->scratch[3],
                 left_im, right_re, ctx->ctxp);
    fmpz_mod_sub(result->coeffs + 0,
                 workspace->scratch[0],
                 workspace->scratch[1], ctx->ctxp);
    fmpz_mod_add(result->coeffs + 1,
                 workspace->scratch[2],
                 workspace->scratch[3], ctx->ctxp);
    cubic_big_normalise(result);
}

static void cubic_big_inv(
    fq_t result, const fq_t input,
    fast_fp2_workspace_t *workspace, const fq_ctx_t ctx)
{
    const fmpz *input_re, *input_im;

    if (!quadratic_basis_is_minus_one) {
        fq_inv(result, input, ctx);
        return;
    }
    fmpz_poly_fit_length(result, 2);
    input_re = cubic_big_component(input, 0, workspace);
    input_im = cubic_big_component(input, 1, workspace);
    fmpz_mod_mul(workspace->scratch[0],
                 input_re, input_re, ctx->ctxp);
    fmpz_mod_mul(workspace->scratch[1],
                 input_im, input_im, ctx->ctxp);
    fmpz_mod_add(workspace->scratch[2],
                 workspace->scratch[0],
                 workspace->scratch[1], ctx->ctxp);
    fmpz_mod_inv(workspace->scratch[3],
                 workspace->scratch[2], ctx->ctxp);
    fmpz_mod_mul(workspace->scratch[4],
                 input_re, workspace->scratch[3], ctx->ctxp);
    fmpz_mod_mul(workspace->scratch[5],
                 input_im, workspace->scratch[3], ctx->ctxp);
    fmpz_set(result->coeffs + 0, workspace->scratch[4]);
    fmpz_mod_neg(result->coeffs + 1,
                 workspace->scratch[5], ctx->ctxp);
    cubic_big_normalise(result);
}

static int cubic_big_base_sqrt(
    fmpz_t result, const fmpz_t input,
    fast_fp2_workspace_t *workspace, const fq_ctx_t ctx)
{
    (void)fmpz_mod_pow_fmpz(
        result, input, workspace->sqrt_exponent, ctx->ctxp);
    fmpz_mod_mul(
        workspace->scratch[13], result, result, ctx->ctxp);
    return fmpz_equal(workspace->scratch[13], input);
}

/*
 * Square root in F_p[i], p = 3 (mod 4).  Reducing the problem to the base
 * field avoids FLINT's generic extension-field square-root machinery.
 */
static int cubic_big_sqrt(
    fq_t result, const fq_t input,
    fast_fp2_workspace_t *workspace, const fq_ctx_t ctx)
{
    const fmpz *input_re, *input_im;

    if (!quadratic_basis_is_minus_one)
        return fq_sqrt(result, input, ctx);

    fmpz_poly_fit_length(result, 2);
    input_re = cubic_big_component(input, 0, workspace);
    input_im = cubic_big_component(input, 1, workspace);
    fmpz_set(workspace->scratch[10], input_re);
    fmpz_set(workspace->scratch[11], input_im);

    /* norm = a^2 + b^2 and s = sqrt(norm). */
    fmpz_mod_mul(workspace->scratch[0],
                 workspace->scratch[10],
                 workspace->scratch[10], ctx->ctxp);
    fmpz_mod_mul(workspace->scratch[1],
                 workspace->scratch[11],
                 workspace->scratch[11], ctx->ctxp);
    fmpz_mod_add(workspace->scratch[2],
                 workspace->scratch[0],
                 workspace->scratch[1], ctx->ctxp);
    if (!cubic_big_base_sqrt(
            workspace->scratch[3],
            workspace->scratch[2], workspace, ctx))
        return 0;

    /* Try x^2 = (a+s)/2, then the other sign if necessary. */
    fmpz_mod_add(workspace->scratch[5],
                 workspace->scratch[10],
                 workspace->scratch[3], ctx->ctxp);
    fmpz_mod_mul(workspace->scratch[5],
                 workspace->scratch[5],
                 workspace->inverse_two, ctx->ctxp);
    if (!cubic_big_base_sqrt(
            workspace->scratch[6],
            workspace->scratch[5], workspace, ctx)) {
        fmpz_mod_sub(workspace->scratch[5],
                     workspace->scratch[10],
                     workspace->scratch[3], ctx->ctxp);
        fmpz_mod_mul(workspace->scratch[5],
                     workspace->scratch[5],
                     workspace->inverse_two, ctx->ctxp);
        if (!cubic_big_base_sqrt(
                workspace->scratch[6],
                workspace->scratch[5], workspace, ctx))
            return 0;
    }

    if (fmpz_is_zero(workspace->scratch[6])) {
        if (!fmpz_is_zero(workspace->scratch[11]))
            return 0;
        fmpz_mod_neg(workspace->scratch[5],
                     workspace->scratch[10], ctx->ctxp);
        if (!cubic_big_base_sqrt(
                workspace->scratch[7],
                workspace->scratch[5], workspace, ctx))
            return 0;
    } else {
        fmpz_mod_add(workspace->scratch[8],
                     workspace->scratch[6],
                     workspace->scratch[6], ctx->ctxp);
        fmpz_mod_inv(workspace->scratch[9],
                     workspace->scratch[8], ctx->ctxp);
        fmpz_mod_mul(workspace->scratch[7],
                     workspace->scratch[11],
                     workspace->scratch[9], ctx->ctxp);
    }

    /* Verify (x+iy)^2 before publishing the result. */
    fmpz_mod_mul(workspace->scratch[0],
                 workspace->scratch[6],
                 workspace->scratch[6], ctx->ctxp);
    fmpz_mod_mul(workspace->scratch[1],
                 workspace->scratch[7],
                 workspace->scratch[7], ctx->ctxp);
    fmpz_mod_sub(workspace->scratch[2],
                 workspace->scratch[0],
                 workspace->scratch[1], ctx->ctxp);
    if (!fmpz_equal(workspace->scratch[2],
                    workspace->scratch[10]))
        return 0;
    fmpz_mod_mul(workspace->scratch[2],
                 workspace->scratch[6],
                 workspace->scratch[7], ctx->ctxp);
    fmpz_mod_add(workspace->scratch[2],
                 workspace->scratch[2],
                 workspace->scratch[2], ctx->ctxp);
    if (!fmpz_equal(workspace->scratch[2],
                    workspace->scratch[11]))
        return 0;

    fmpz_set(result->coeffs + 0, workspace->scratch[6]);
    fmpz_set(result->coeffs + 1, workspace->scratch[7]);
    cubic_big_normalise(result);
    return 1;
}

/*
 * Multiply two residues modulo
 *   Y^3 + c2 Y^2 + c1 Y + c0.
 * The degree-three and degree-four reduction coefficients are
 *   (-c0,-c1,-c2) and
 *   (c2*c0,c2*c1-c0,c2^2-c1).
 */
static void cubic_residue_mul_big(
    fq_t out0, fq_t out1, fq_t out2,
    const fq_t a0, const fq_t a1, const fq_t a2,
    const fq_t b0, const fq_t b1, const fq_t b2,
    const fq_t c0, const fq_t c1, const fq_t c2,
    const fq_t r40, const fq_t r41,
    const fq_t r42, cubic_big_workspace_t *workspace,
    const fq_ctx_t ctx)
{
    cubic_big_mul(workspace->t0, a0, b0, &workspace->fast, ctx);
    cubic_big_mul(workspace->t1, a0, b1, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, a1, b0, &workspace->fast, ctx);
    cubic_big_add(workspace->t1, workspace->t1, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->t2, a0, b2, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, a1, b1, &workspace->fast, ctx);
    cubic_big_add(workspace->t2, workspace->t2, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, a2, b0, &workspace->fast, ctx);
    cubic_big_add(workspace->t2, workspace->t2, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->t3, a1, b2, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, a2, b1, &workspace->fast, ctx);
    cubic_big_add(workspace->t3, workspace->t3, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->t4, a2, b2, &workspace->fast, ctx);

    cubic_big_mul(workspace->s0, c0, workspace->t3, &workspace->fast, ctx);
    cubic_big_sub(out0, workspace->t0, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, r40, workspace->t4, &workspace->fast, ctx);
    cubic_big_add(out0, out0, workspace->s0, &workspace->fast, ctx);

    cubic_big_mul(workspace->s0, c1, workspace->t3, &workspace->fast, ctx);
    cubic_big_sub(out1, workspace->t1, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, r41, workspace->t4, &workspace->fast, ctx);
    cubic_big_add(out1, out1, workspace->s0, &workspace->fast, ctx);

    cubic_big_mul(workspace->s0, c2, workspace->t3, &workspace->fast, ctx);
    cubic_big_sub(out2, workspace->t2, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, r42, workspace->t4, &workspace->fast, ctx);
    cubic_big_add(out2, out2, workspace->s0, &workspace->fast, ctx);
}

static void cubic_residue_square_big(
    fq_t out0, fq_t out1, fq_t out2,
    const fq_t a0, const fq_t a1, const fq_t a2,
    const fq_t c0, const fq_t c1, const fq_t c2,
    const fq_t r40, const fq_t r41,
    const fq_t r42, cubic_big_workspace_t *workspace,
    const fq_ctx_t ctx)
{
    cubic_big_mul(workspace->t0, a0, a0, &workspace->fast, ctx);
    cubic_big_mul(workspace->t1, a0, a1, &workspace->fast, ctx);
    cubic_big_add(workspace->t1, workspace->t1, workspace->t1, &workspace->fast, ctx);
    cubic_big_mul(workspace->t2, a0, a2, &workspace->fast, ctx);
    cubic_big_add(workspace->t2, workspace->t2, workspace->t2, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, a1, a1, &workspace->fast, ctx);
    cubic_big_add(workspace->t2, workspace->t2, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->t3, a1, a2, &workspace->fast, ctx);
    cubic_big_add(workspace->t3, workspace->t3, workspace->t3, &workspace->fast, ctx);
    cubic_big_mul(workspace->t4, a2, a2, &workspace->fast, ctx);

    cubic_big_mul(workspace->s0, c0, workspace->t3, &workspace->fast, ctx);
    cubic_big_sub(out0, workspace->t0, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, r40, workspace->t4, &workspace->fast, ctx);
    cubic_big_add(out0, out0, workspace->s0, &workspace->fast, ctx);

    cubic_big_mul(workspace->s0, c1, workspace->t3, &workspace->fast, ctx);
    cubic_big_sub(out1, workspace->t1, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, r41, workspace->t4, &workspace->fast, ctx);
    cubic_big_add(out1, out1, workspace->s0, &workspace->fast, ctx);

    cubic_big_mul(workspace->s0, c2, workspace->t3, &workspace->fast, ctx);
    cubic_big_sub(out2, workspace->t2, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, r42, workspace->t4, &workspace->fast, ctx);
    cubic_big_add(out2, out2, workspace->s0, &workspace->fast, ctx);
}

/*
 * Find one root of f projectively from h = g +/- 1.  For quadratic h,
 * the formulas below are the pseudo-remainder of f by h, scaled by lc(h)^2:
 *
 * r1 = A^2*c1 - A*C - A*c2*B + B^2,
 * r0 = A^2*c0 - A*c2*C + B*C.
 *
 * This avoids every inversion in the split step.
 */
static int cubic_projective_root_big(
    fq_t numerator, fq_t denominator,
    const fq_t c0, const fq_t c1, const fq_t c2,
    const fq_t h0, const fq_t h1, const fq_t h2,
    cubic_big_workspace_t *workspace, const fq_ctx_t ctx)
{
    if (fq_is_zero(h2, ctx)) {
        if (fq_is_zero(h1, ctx))
            return 0;
        cubic_big_neg(numerator, h0, &workspace->fast, ctx);
        fq_set(denominator, h1, ctx);
        return 1;
    }

    /* h0 is allowed to alias workspace->t0 at the call site. */
    fq_set(workspace->s1, h0, ctx);
    cubic_big_mul(workspace->t0, h2, h2, &workspace->fast, ctx);

    cubic_big_mul(workspace->t1, workspace->t0, c1, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, h2, workspace->s1, &workspace->fast, ctx);
    cubic_big_sub(workspace->t1, workspace->t1, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, h2, c2, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, workspace->s0, h1, &workspace->fast, ctx);
    cubic_big_sub(workspace->t1, workspace->t1, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, h1, h1, &workspace->fast, ctx);
    cubic_big_add(workspace->t1, workspace->t1, workspace->s0, &workspace->fast, ctx);

    cubic_big_mul(workspace->t2, workspace->t0, c0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, h2, c2, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, workspace->s0, workspace->s1, &workspace->fast, ctx);
    cubic_big_sub(workspace->t2, workspace->t2, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, h1, workspace->s1, &workspace->fast, ctx);
    cubic_big_add(workspace->t2, workspace->t2, workspace->s0, &workspace->fast, ctx);

    if (!fq_is_zero(workspace->t1, ctx)) {
        cubic_big_neg(numerator, workspace->t2, &workspace->fast, ctx);
        fq_set(denominator, workspace->t1, ctx);
        return 1;
    }
    if (!fq_is_zero(workspace->t2, ctx))
        return 0;

    /* h divides f; the complementary linear factor has root B/A-c2. */
    cubic_big_mul(workspace->s0, h2, c2, &workspace->fast, ctx);
    cubic_big_sub(numerator, h1, workspace->s0, &workspace->fast, ctx);
    fq_set(denominator, h2, ctx);
    return 1;
}

static int cubic_projective_root_is_valid_big(
    const fq_t numerator, const fq_t denominator,
    const fq_t c0, const fq_t c1, const fq_t c2,
    cubic_big_workspace_t *workspace, const fq_ctx_t ctx)
{
    /* n^3 + c2*n^2*d + c1*n*d^2 + c0*d^3 */
    cubic_big_mul(workspace->t0, numerator, numerator, &workspace->fast, ctx);
    cubic_big_mul(workspace->t1, workspace->t0, numerator, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, c2, workspace->t0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, workspace->s0, denominator, &workspace->fast, ctx);
    cubic_big_add(workspace->t1, workspace->t1, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->t2, denominator, denominator, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, c1, numerator, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, workspace->s0, workspace->t2, &workspace->fast, ctx);
    cubic_big_add(workspace->t1, workspace->t1, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->t2, workspace->t2, denominator, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, c0, workspace->t2, &workspace->fast, ctx);
    cubic_big_add(workspace->t1, workspace->t1, workspace->s0, &workspace->fast, ctx);
    return fq_is_zero(workspace->t1, ctx);
}

static int finish_cubic_roots_big(
    jkey_t *roots, size_t *count_out, size_t capacity,
    const fq_t root, const fq_t c0,
    const fq_t c1, const fq_t c2,
    cubic_big_workspace_t *workspace, const fq_ctx_t ctx)
{
    size_t count = 0;

    if (capacity < 3)
        return 0;
    /* Check f(root) before using synthetic division. */
    cubic_big_mul(workspace->t0, root, root, &workspace->fast, ctx);
    cubic_big_mul(workspace->t1, workspace->t0, root, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, c2, workspace->t0, &workspace->fast, ctx);
    cubic_big_add(workspace->t1, workspace->t1, workspace->s0, &workspace->fast, ctx);
    cubic_big_mul(workspace->s0, c1, root, &workspace->fast, ctx);
    cubic_big_add(workspace->t1, workspace->t1, workspace->s0, &workspace->fast, ctx);
    cubic_big_add(workspace->t1, workspace->t1, c0, &workspace->fast, ctx);
    if (!fq_is_zero(workspace->t1, ctx))
        return 0;

    fq_to_key(roots + count++, root, workspace->representative, ctx);
    cubic_big_add(workspace->q1, c2, root, &workspace->fast, ctx);
    cubic_big_mul(workspace->q0, root, workspace->q1, &workspace->fast, ctx);
    cubic_big_add(workspace->q0, workspace->q0, c1, &workspace->fast, ctx);
    cubic_big_mul(workspace->discriminant,
                workspace->q1, workspace->q1, &workspace->fast, ctx);
    cubic_big_add(workspace->t0, workspace->q0, workspace->q0, &workspace->fast, ctx);
    cubic_big_add(workspace->t0, workspace->t0, workspace->t0, &workspace->fast, ctx);
    cubic_big_sub(workspace->discriminant,
                workspace->discriminant, workspace->t0, &workspace->fast, ctx);
    if (!cubic_big_sqrt(
            workspace->square_root, workspace->discriminant, &workspace->fast, ctx))
        return 0;

    cubic_big_neg(workspace->q1, workspace->q1, &workspace->fast, ctx);
    cubic_big_add(workspace->t0,
                workspace->q1, workspace->square_root, &workspace->fast, ctx);
    cubic_big_mul(workspace->t0,
                workspace->t0, workspace->inverse_two, &workspace->fast, ctx);
    fq_to_key(
        roots + count++, workspace->t0, workspace->representative, ctx);
    if (!fq_is_zero(workspace->square_root, ctx)) {
        cubic_big_sub(workspace->t0,
                    workspace->q1, workspace->square_root, &workspace->fast, ctx);
        cubic_big_mul(workspace->t0,
                    workspace->t0, workspace->inverse_two, &workspace->fast, ctx);
        fq_to_key(
            roots + count++, workspace->t0,
            workspace->representative, ctx);
    }
    *count_out = sort_unique_root_keys(roots, count);
    return 1;
}

static size_t fallback_cubic_roots_big(
    jkey_t *roots, size_t capacity, const fq_poly_t polynomial,
    cubic_big_workspace_t *workspace, const fq_ctx_t ctx)
{
    size_t count = 0;
    slong factor_index;

    fq_poly_roots(
        workspace->fallback_factors, polynomial, 0, ctx);
    for (factor_index = 0;
         factor_index < workspace->fallback_factors->num;
         factor_index++) {
        if (fq_poly_degree(
                workspace->fallback_factors->poly + factor_index,
                ctx) != 1)
            continue;
        if (count == capacity)
            break;
        fq_poly_get_coeff(
            workspace->t0,
            workspace->fallback_factors->poly + factor_index, 0, ctx);
        cubic_big_neg(workspace->t0, workspace->t0, &workspace->fast, ctx);
        fq_to_key(
            roots + count++, workspace->t0,
            workspace->representative, ctx);
    }
    return sort_unique_root_keys(roots, count);
}

static void find_cubic_roots_task_big(slong task_index, void *argument)
{
    cubic_big_root_tasks_t *tasks = argument;
    cubic_big_workspace_t workspace;
    size_t first, lane_count, lane, attempt;
    size_t actual_indices[CUBIC_BATCH_LANES];
    unsigned char resolved[CUBIC_BATCH_LANES];
    unsigned char usable[CUBIC_BATCH_LANES];
    uint64_t local_count = 0, local_cpu_ns = 0, local_fallbacks = 0;

    cubic_big_workspace_init(&workspace, task_index, tasks->ctx);
    for (;;) {
        uint64_t group_start;
        slong bit;

        first = atomic_fetch_add_explicit(
            &tasks->next_index, tasks->lane_width,
            memory_order_relaxed);
        if (first >= tasks->count)
            break;
        lane_count = tasks->count - first;
        if (lane_count > tasks->lane_width)
            lane_count = tasks->lane_width;
        group_start = thread_cpu_nanoseconds();

        for (lane = 0; lane < lane_count; lane++) {
            size_t actual_index = tasks->indices[first + lane];

            actual_indices[lane] = actual_index;
            resolved[lane] = 0;
            usable[lane] = 0;
            fq_poly_get_coeff(
                workspace.c0 + lane, tasks->evaluated + actual_index,
                0, tasks->ctx);
            fq_poly_get_coeff(
                workspace.c1 + lane, tasks->evaluated + actual_index,
                1, tasks->ctx);
            fq_poly_get_coeff(
                workspace.c2 + lane, tasks->evaluated + actual_index,
                2, tasks->ctx);
            cubic_big_mul(workspace.r40 + lane,
                        workspace.c2 + lane, workspace.c0 + lane, &workspace.fast, tasks->ctx);
            cubic_big_mul(workspace.r41 + lane,
                        workspace.c2 + lane, workspace.c1 + lane, &workspace.fast, tasks->ctx);
            cubic_big_sub(workspace.r41 + lane,
                        workspace.r41 + lane, workspace.c0 + lane, &workspace.fast, tasks->ctx);
            cubic_big_mul(workspace.r42 + lane,
                        workspace.c2 + lane, workspace.c2 + lane, &workspace.fast, tasks->ctx);
            cubic_big_sub(workspace.r42 + lane,
                        workspace.r42 + lane, workspace.c1 + lane, &workspace.fast, tasks->ctx);
        }

        for (attempt = 0; attempt < CUBIC_MAX_SPLIT_ATTEMPTS;
             attempt++) {
            size_t unresolved = 0;

            for (lane = 0; lane < lane_count; lane++) {
                if (resolved[lane])
                    continue;
                unresolved++;
                fq_rand(
                    workspace.base[0] + lane,
                    workspace.randstate, tasks->ctx);
                fq_rand_not_zero(
                    workspace.base[1] + lane,
                    workspace.randstate, tasks->ctx);
                fq_zero(workspace.base[2] + lane, tasks->ctx);
                fq_one(
                    workspace.accumulator[0] + lane, tasks->ctx);
                fq_zero(
                    workspace.accumulator[1] + lane, tasks->ctx);
                fq_zero(
                    workspace.accumulator[2] + lane, tasks->ctx);
            }
            if (unresolved == 0)
                break;

            for (bit = (slong)fmpz_bits(workspace.halfq) - 1;
                 bit >= 0; bit--) {
                for (lane = 0; lane < lane_count; lane++) {
                    if (resolved[lane])
                        continue;
                    cubic_residue_square_big(
                        workspace.power_scratch[0] + lane,
                        workspace.power_scratch[1] + lane,
                        workspace.power_scratch[2] + lane,
                        workspace.accumulator[0] + lane,
                        workspace.accumulator[1] + lane,
                        workspace.accumulator[2] + lane,
                        workspace.c0 + lane, workspace.c1 + lane,
                        workspace.c2 + lane, workspace.r40 + lane,
                        workspace.r41 + lane, workspace.r42 + lane,
                        &workspace, tasks->ctx);
                    fq_swap(
                        workspace.accumulator[0] + lane,
                        workspace.power_scratch[0] + lane, tasks->ctx);
                    fq_swap(
                        workspace.accumulator[1] + lane,
                        workspace.power_scratch[1] + lane, tasks->ctx);
                    fq_swap(
                        workspace.accumulator[2] + lane,
                        workspace.power_scratch[2] + lane, tasks->ctx);
                }
                if (!fmpz_tstbit(workspace.halfq, (ulong)bit))
                    continue;
                for (lane = 0; lane < lane_count; lane++) {
                    if (resolved[lane])
                        continue;
                    cubic_residue_mul_big(
                        workspace.power_scratch[0] + lane,
                        workspace.power_scratch[1] + lane,
                        workspace.power_scratch[2] + lane,
                        workspace.accumulator[0] + lane,
                        workspace.accumulator[1] + lane,
                        workspace.accumulator[2] + lane,
                        workspace.base[0] + lane,
                        workspace.base[1] + lane,
                        workspace.base[2] + lane,
                        workspace.c0 + lane, workspace.c1 + lane,
                        workspace.c2 + lane, workspace.r40 + lane,
                        workspace.r41 + lane, workspace.r42 + lane,
                        &workspace, tasks->ctx);
                    fq_swap(
                        workspace.accumulator[0] + lane,
                        workspace.power_scratch[0] + lane, tasks->ctx);
                    fq_swap(
                        workspace.accumulator[1] + lane,
                        workspace.power_scratch[1] + lane, tasks->ctx);
                    fq_swap(
                        workspace.accumulator[2] + lane,
                        workspace.power_scratch[2] + lane, tasks->ctx);
                }
            }

            for (lane = 0; lane < lane_count; lane++) {
                int found;

                if (resolved[lane])
                    continue;
                fq_sub_one(
                    workspace.t0, workspace.accumulator[0] + lane,
                    tasks->ctx);
                found = cubic_projective_root_big(
                    workspace.numerator + lane,
                    workspace.denominator + lane,
                    workspace.c0 + lane, workspace.c1 + lane,
                    workspace.c2 + lane, workspace.t0,
                    workspace.accumulator[1] + lane,
                    workspace.accumulator[2] + lane,
                    &workspace, tasks->ctx);
                if (!found) {
                    cubic_big_add(
                        workspace.t0, workspace.accumulator[0] + lane,
                        workspace.one, &workspace.fast, tasks->ctx);
                    found = cubic_projective_root_big(
                        workspace.numerator + lane,
                        workspace.denominator + lane,
                        workspace.c0 + lane, workspace.c1 + lane,
                        workspace.c2 + lane, workspace.t0,
                        workspace.accumulator[1] + lane,
                        workspace.accumulator[2] + lane,
                        &workspace, tasks->ctx);
                }
                if (found &&
                    cubic_projective_root_is_valid_big(
                        workspace.numerator + lane,
                        workspace.denominator + lane,
                        workspace.c0 + lane, workspace.c1 + lane,
                        workspace.c2 + lane, &workspace, tasks->ctx))
                    resolved[lane] = usable[lane] = 1;
            }
        }

        /* Montgomery batch inversion of all successful denominators. */
        fq_one(workspace.product, tasks->ctx);
        for (lane = 0; lane < lane_count; lane++) {
            if (!usable[lane])
                continue;
            fq_set(
                workspace.prefix + lane, workspace.product, tasks->ctx);
            cubic_big_mul(
                workspace.product, workspace.product,
                workspace.denominator + lane, &workspace.fast, tasks->ctx);
        }
        cubic_big_inv(
            workspace.inverse_product, workspace.product, &workspace.fast, tasks->ctx);
        for (lane = lane_count; lane-- > 0;) {
            if (!usable[lane])
                continue;
            cubic_big_mul(
                workspace.t0, workspace.inverse_product,
                workspace.prefix + lane, &workspace.fast, tasks->ctx);
            cubic_big_mul(
                workspace.affine_root + lane,
                workspace.numerator + lane, workspace.t0, &workspace.fast, tasks->ctx);
            cubic_big_mul(
                workspace.inverse_product, workspace.inverse_product,
                workspace.denominator + lane, &workspace.fast, tasks->ctx);
        }

        for (lane = 0; lane < lane_count; lane++) {
            size_t actual_index = actual_indices[lane];
            jkey_t *point_roots =
                tasks->root_keys + actual_index * tasks->root_stride;
            size_t count = 0;

            if (!usable[lane] ||
                !finish_cubic_roots_big(
                    point_roots, &count, tasks->root_stride,
                    workspace.affine_root + lane,
                    workspace.c0 + lane, workspace.c1 + lane,
                    workspace.c2 + lane, &workspace, tasks->ctx)) {
                count = fallback_cubic_roots_big(
                    point_roots, tasks->root_stride,
                    tasks->evaluated + actual_index,
                    &workspace, tasks->ctx);
                local_fallbacks++;
            }
            tasks->root_counts[actual_index] = count;
        }
        local_count += lane_count;
        local_cpu_ns += thread_cpu_nanoseconds() - group_start;
    }

    atomic_fetch_add_explicit(
        &cubic_batched_counts, local_count, memory_order_relaxed);
    atomic_fetch_add_explicit(
        &cubic_batched_cpu_ns, local_cpu_ns, memory_order_relaxed);
    atomic_fetch_add_explicit(
        &cubic_batched_fallbacks, local_fallbacks, memory_order_relaxed);
    if (root_profile_enabled &&
        cubic_kernel == CUBIC_KERNEL_BATCHED) {
        atomic_fetch_add_explicit(
            root_profile_counts + 3, local_count,
            memory_order_relaxed);
        atomic_fetch_add_explicit(
            root_profile_cpu_ns + 3, local_cpu_ns,
            memory_order_relaxed);
    }
    cubic_big_workspace_clear(&workspace, tasks->ctx);
}

static void find_roots_task_word(slong task_index, void *argument)
{
    word_root_tasks_t *tasks = argument;
    fq_nmod_poly_factor_t roots;
    fq_nmod_poly_t factor_a, factor_b, temporary, inverse_series;
    fq_nmod_poly_struct factor_stack[ROOT_FACTOR_STACK_SIZE];
    fq_nmod_t b, c, discriminant, square_root, root, inverse_two;
    fmpz_t halfq;
    flint_rand_t randstate;
    nmod_poly_t representative;
    size_t first, last, i, stack_index;
    uint64_t *local_counts = NULL, *local_cpu_ns = NULL;
    uint64_t seed;

    if (root_profile_enabled) {
        local_counts = xcalloc(tasks->root_stride + 1,
                               sizeof(*local_counts));
        local_cpu_ns = xcalloc(tasks->root_stride + 1,
                               sizeof(*local_cpu_ns));
    }
    fq_nmod_poly_factor_init(roots, tasks->ctx);
    fq_nmod_poly_init(factor_a, tasks->ctx);
    fq_nmod_poly_init(factor_b, tasks->ctx);
    fq_nmod_poly_init(temporary, tasks->ctx);
    fq_nmod_poly_init(inverse_series, tasks->ctx);
    for (stack_index = 0; stack_index < ROOT_FACTOR_STACK_SIZE;
         stack_index++)
        fq_nmod_poly_init(factor_stack + stack_index, tasks->ctx);
    fq_nmod_init(b, tasks->ctx);
    fq_nmod_init(c, tasks->ctx);
    fq_nmod_init(discriminant, tasks->ctx);
    fq_nmod_init(square_root, tasks->ctx);
    fq_nmod_init(root, tasks->ctx);
    fq_nmod_init(inverse_two, tasks->ctx);
    fq_nmod_set_ui(inverse_two, 2, tasks->ctx);
    fq_nmod_inv(inverse_two, inverse_two, tasks->ctx);
    fmpz_init(halfq);
    fq_nmod_ctx_order(halfq, tasks->ctx);
    fmpz_sub_ui(halfq, halfq, 1);
    fmpz_fdiv_q_2exp(halfq, halfq, 1);
    ISOGENY_RAND_INIT(randstate);
    seed = mix64(UINT64_C(0x9e3779b97f4a7c15) +
                 (uint64_t)(task_index + 1));
    ISOGENY_RAND_SET_SEED(
        randstate, (ulong)seed,
        (ulong)mix64(seed ^ UINT64_C(0xd1b54a32d192ed03)));
    nmod_poly_init(representative, field_characteristic_word);

    for (;;) {
        first = atomic_fetch_add_explicit(
            &tasks->next_index, tasks->claim_size, memory_order_relaxed);
        if (first >= tasks->count)
            break;
        last = first + tasks->claim_size;
        if (last > tasks->count)
            last = tasks->count;
        for (i = first; i < last; i++) {
            size_t actual_index = tasks->indices == NULL
                                      ? i
                                      : tasks->indices[i];
            slong degree = fq_nmod_poly_degree(
                tasks->evaluated + actual_index, tasks->ctx);
            size_t count;
            uint64_t profile_start = 0;
            jkey_t *point_roots =
                tasks->root_keys + actual_index * tasks->root_stride;

            if (tasks->skip_degree_three && degree == 3)
                continue;
            if (root_profile_enabled)
                profile_start = thread_cpu_nanoseconds();
            count = fully_split_roots_word(
                point_roots, tasks->root_stride,
                tasks->evaluated + actual_index,
                halfq, randstate, factor_a, factor_b, temporary,
                inverse_series, factor_stack, b, c, discriminant,
                square_root, root, inverse_two, representative, tasks->ctx);
            if (count == SIZE_MAX) {
                slong factor_index;

                count = 0;
                fq_nmod_poly_roots(
                    roots, tasks->evaluated + actual_index, 0, tasks->ctx);
                for (factor_index = 0; factor_index < roots->num;
                     factor_index++) {
                    if (fq_nmod_poly_degree(
                            roots->poly + factor_index,
                            tasks->ctx) != 1)
                        continue;
                    if (count == tasks->root_stride)
                        break;
                    fq_nmod_poly_get_coeff(
                        c, roots->poly + factor_index, 0, tasks->ctx);
                    fq_nmod_neg(root, c, tasks->ctx);
                    fq_nmod_to_key(
                        point_roots + count++, root, representative,
                        tasks->ctx);
                }
            }
            qsort(point_roots, count, sizeof(*point_roots), compare_keys);
            if (count > 1) {
                size_t read_index, unique_count = 1;

                for (read_index = 1; read_index < count; read_index++) {
                    if (!key_equal(point_roots[read_index],
                                   point_roots[unique_count - 1])) {
                        if (unique_count != read_index)
                            key_set(point_roots + unique_count,
                                    point_roots + read_index);
                        unique_count++;
                    }
                }
                count = unique_count;
            }
            tasks->root_counts[actual_index] = count;
            if (root_profile_enabled && degree >= 0 &&
                (size_t)degree <= tasks->root_stride) {
                local_counts[degree]++;
                local_cpu_ns[degree] +=
                    thread_cpu_nanoseconds() - profile_start;
            }
        }
    }

    if (root_profile_enabled) {
        size_t degree;

        for (degree = 0;
             degree <= tasks->root_stride &&
             degree < ROOT_PROFILE_DEGREES;
             degree++) {
            if (local_counts[degree] == 0)
                continue;
            atomic_fetch_add_explicit(
                root_profile_counts + degree, local_counts[degree],
                memory_order_relaxed);
            atomic_fetch_add_explicit(
                root_profile_cpu_ns + degree, local_cpu_ns[degree],
                memory_order_relaxed);
        }
        free(local_cpu_ns);
        free(local_counts);
    }
    nmod_poly_clear(representative);
    ISOGENY_RAND_CLEAR(randstate);
    fmpz_clear(halfq);
    fq_nmod_clear(inverse_two, tasks->ctx);
    fq_nmod_clear(root, tasks->ctx);
    fq_nmod_clear(square_root, tasks->ctx);
    fq_nmod_clear(discriminant, tasks->ctx);
    fq_nmod_clear(c, tasks->ctx);
    fq_nmod_clear(b, tasks->ctx);
    for (stack_index = ROOT_FACTOR_STACK_SIZE; stack_index-- > 0;)
        fq_nmod_poly_clear(factor_stack + stack_index, tasks->ctx);
    fq_nmod_poly_clear(inverse_series, tasks->ctx);
    fq_nmod_poly_clear(temporary, tasks->ctx);
    fq_nmod_poly_clear(factor_b, tasks->ctx);
    fq_nmod_poly_clear(factor_a, tasks->ctx);
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

static void build_evaluation_tree_task_big(slong shard_index,
                                           void *argument)
{
    big_sharded_evaluation_tasks_t *tasks = argument;
    size_t shard = (size_t)shard_index;
    size_t count = tasks->counts[shard];

    tasks->trees[shard] =
        _fq_poly_tree_alloc((slong)count, tasks->base->ctx);
    _fq_poly_tree_build(
        tasks->trees[shard],
        tasks->base->points + tasks->offsets[shard],
        (slong)count, tasks->base->ctx);
}

static void free_evaluation_tree_task_big(slong shard_index,
                                          void *argument)
{
    big_sharded_evaluation_tasks_t *tasks = argument;
    size_t shard = (size_t)shard_index;

    _fq_poly_tree_free(
        tasks->trees[shard], (slong)tasks->counts[shard],
        tasks->base->ctx);
}

static void evaluate_coefficient_shard_task_big(slong task_index,
                                                void *argument)
{
    big_sharded_evaluation_tasks_t *tasks = argument;
    size_t task = (size_t)task_index;
    size_t y = task / tasks->shard_count;
    size_t shard = task % tasks->shard_count;
    size_t offset = tasks->offsets[shard];
    size_t count = tasks->counts[shard];
    const fq_poly_struct *coefficient =
        tasks->base->phi->coefficient_x.big + y;
    fq_struct *values =
        tasks->base->coefficient_values + y * tasks->base->count + offset;
    size_t i;

    if (tasks->initialize_values)
        for (i = 0; i < count; i++)
            fq_init(values + i, tasks->base->ctx);

    if (fq_poly_is_zero(coefficient, tasks->base->ctx))
        return;
    if (tasks->base->use_fast) {
        _fq_poly_evaluate_fq_vec_fast_precomp(
            values, coefficient->coeffs, coefficient->length,
            tasks->trees[shard], (slong)count, tasks->base->ctx);
    } else {
        _fq_poly_evaluate_fq_vec_iter(
            values, coefficient->coeffs, coefficient->length,
            tasks->base->points + offset, (slong)count,
            tasks->base->ctx);
    }
}

static void clear_coefficient_shard_task_big(slong task_index,
                                             void *argument)
{
    big_sharded_evaluation_tasks_t *tasks = argument;
    size_t task = (size_t)task_index;
    size_t y = task / tasks->shard_count;
    size_t shard = task % tasks->shard_count;
    size_t offset = tasks->offsets[shard];
    size_t count = tasks->counts[shard];
    fq_struct *values =
        tasks->base->coefficient_values + y * tasks->base->count + offset;
    size_t i;

    for (i = 0; i < count; i++)
        fq_clear(values + i, tasks->base->ctx);
}

static void evaluate_batch_parallel_big(
    const big_evaluation_tasks_t *base, size_t coefficient_count,
    unsigned n_threads)
{
    big_sharded_evaluation_tasks_t tasks;
    size_t *storage;
    size_t shard_count, shard, base_size, remainder, offset = 0;
    size_t task_count;

    if (base->count == 0 || coefficient_count == 0)
        return;
    shard_count =
        evaluation_shard_count(base->count, n_threads, base->use_fast);
    if (coefficient_count > SIZE_MAX / shard_count)
        die("evaluation task count overflow");
    task_count = coefficient_count * shard_count;
    if (task_count > (size_t)LONG_MAX)
        die("evaluation task count exceeds FLINT's slong interface");

    storage = xrealloc(NULL, 2 * shard_count, sizeof(*storage));
    tasks.base = base;
    tasks.shard_count = shard_count;
    tasks.offsets = storage;
    tasks.counts = storage + shard_count;
    tasks.trees = NULL;
    tasks.initialize_values = parallel_coefficient_lifecycle(
        base->count, coefficient_count, n_threads);
    if (!tasks.initialize_values) {
        size_t i;

        for (i = 0; i < base->count * coefficient_count; i++)
            fq_init(base->coefficient_values + i, base->ctx);
    }
    base_size = base->count / shard_count;
    remainder = base->count % shard_count;
    for (shard = 0; shard < shard_count; shard++) {
        size_t count = base_size + (shard < remainder ? 1 : 0);

        ((size_t *)tasks.offsets)[shard] = offset;
        ((size_t *)tasks.counts)[shard] = count;
        offset += count;
    }
    if (offset != base->count)
        die("internal error partitioning big evaluation points");

    if (base->use_fast) {
        tasks.trees =
            xrealloc(NULL, shard_count, sizeof(*tasks.trees));
        if (shard_count == 1)
            build_evaluation_tree_task_big(0, &tasks);
        else
            flint_parallel_do(
                build_evaluation_tree_task_big, &tasks,
                (slong)shard_count, (int)n_threads,
                FLINT_PARALLEL_UNIFORM);
    }
    flint_parallel_do(
        evaluate_coefficient_shard_task_big, &tasks,
        (slong)task_count, (int)n_threads,
        FLINT_PARALLEL_UNIFORM);
    if (base->use_fast) {
        if (shard_count == 1)
            free_evaluation_tree_task_big(0, &tasks);
        else
            flint_parallel_do(
                free_evaluation_tree_task_big, &tasks,
                (slong)shard_count, (int)n_threads,
                FLINT_PARALLEL_UNIFORM);
        free(tasks.trees);
    }
    free(storage);
}

static void clear_batch_parallel_big(
    const big_evaluation_tasks_t *base, size_t coefficient_count,
    unsigned n_threads)
{
    big_sharded_evaluation_tasks_t tasks;
    size_t *storage;
    size_t shard_count, shard, base_size, remainder, offset = 0;
    size_t task_count;

    if (base->count == 0 || coefficient_count == 0)
        return;
    if (!parallel_coefficient_lifecycle(
            base->count, coefficient_count, n_threads)) {
        size_t i;

        for (i = 0; i < base->count * coefficient_count; i++)
            fq_clear(base->coefficient_values + i, base->ctx);
        return;
    }
    shard_count =
        evaluation_shard_count(base->count, n_threads, base->use_fast);
    if (coefficient_count > SIZE_MAX / shard_count)
        die("coefficient-clear task count overflow");
    task_count = coefficient_count * shard_count;
    storage = xrealloc(NULL, 2 * shard_count, sizeof(*storage));
    tasks.base = base;
    tasks.shard_count = shard_count;
    tasks.offsets = storage;
    tasks.counts = storage + shard_count;
    tasks.trees = NULL;
    tasks.initialize_values = 0;
    base_size = base->count / shard_count;
    remainder = base->count % shard_count;
    for (shard = 0; shard < shard_count; shard++) {
        size_t count = base_size + (shard < remainder ? 1 : 0);

        ((size_t *)tasks.offsets)[shard] = offset;
        ((size_t *)tasks.counts)[shard] = count;
        offset += count;
    }
    flint_parallel_do(
        clear_coefficient_shard_task_big, &tasks,
        (slong)task_count, (int)n_threads,
        FLINT_PARALLEL_UNIFORM);
    free(storage);
}

static int append_small_roots_big(
    jkey_t *roots, size_t *count, size_t capacity,
    const fq_poly_t polynomial, fq_t b, fq_t c, fq_t discriminant,
    fq_t square_root, fq_t root, const fq_t inverse_two,
    fmpz_poly_t representative, fast_fp2_workspace_t *fast_workspace,
    const fq_ctx_t ctx)
{
    slong degree = fq_poly_degree(polynomial, ctx);

    if (degree < 0)
        return 1;
    if (degree == 0)
        return 1;
    if (degree == 1) {
        if (*count == capacity)
            return 0;
        fq_poly_get_coeff(c, polynomial, 0, ctx);
        cubic_big_neg(root, c, fast_workspace, ctx);
        fq_to_key(roots + (*count)++, root, representative, ctx);
        return 1;
    }
    if (degree != 2)
        return 0;

    fq_poly_get_coeff(c, polynomial, 0, ctx);
    fq_poly_get_coeff(b, polynomial, 1, ctx);
    cubic_big_mul(discriminant, b, b, fast_workspace, ctx);
    cubic_big_add(root, c, c, fast_workspace, ctx);
    cubic_big_add(root, root, root, fast_workspace, ctx);
    cubic_big_sub(discriminant, discriminant, root, fast_workspace, ctx);
    if (!cubic_big_sqrt(square_root, discriminant, fast_workspace, ctx))
        return 0;

    cubic_big_neg(b, b, fast_workspace, ctx);
    cubic_big_add(root, b, square_root, fast_workspace, ctx);
    cubic_big_mul(root, root, inverse_two, fast_workspace, ctx);
    if (*count == capacity)
        return 0;
    fq_to_key(roots + (*count)++, root, representative, ctx);
    if (!fq_is_zero(square_root, ctx)) {
        cubic_big_sub(root, b, square_root, fast_workspace, ctx);
        cubic_big_mul(root, root, inverse_two, fast_workspace, ctx);
        if (*count == capacity)
            return 0;
        fq_to_key(roots + (*count)++, root, representative, ctx);
    }
    return 1;
}

static void find_quadratic_roots_task_big(
    slong task_index, void *argument)
{
    quadratic_big_root_tasks_t *tasks = argument;
    fq_poly_factor_t fallback_factors;
    fq_t b, c, discriminant, square_root, root, inverse_two;
    fmpz_poly_t representative;
    fast_fp2_workspace_t fast_workspace;
    size_t first, last, i;
    uint64_t local_count = 0, local_cpu_ns = 0;

    (void)task_index;
    fq_poly_factor_init(fallback_factors, tasks->ctx);
    fq_init(b, tasks->ctx);
    fq_init(c, tasks->ctx);
    fq_init(discriminant, tasks->ctx);
    fq_init(square_root, tasks->ctx);
    fq_init(root, tasks->ctx);
    fq_init(inverse_two, tasks->ctx);
    fq_set_ui(inverse_two, 2, tasks->ctx);
    fq_inv(inverse_two, inverse_two, tasks->ctx);
    fmpz_poly_init(representative);
    fast_fp2_workspace_init(&fast_workspace, tasks->ctx);

    for (;;) {
        first = atomic_fetch_add_explicit(
            &tasks->next_index, tasks->claim_size,
            memory_order_relaxed);
        if (first >= tasks->count)
            break;
        last = first + tasks->claim_size;
        if (last > tasks->count)
            last = tasks->count;
        for (i = first; i < last; i++) {
            size_t actual_index = tasks->indices[i];
            jkey_t *point_roots =
                tasks->root_keys + actual_index * tasks->root_stride;
            size_t count = 0;
            uint64_t start = root_profile_enabled
                                 ? thread_cpu_nanoseconds()
                                 : 0;

            if (!append_small_roots_big(
                    point_roots, &count, tasks->root_stride,
                    tasks->evaluated + actual_index, b, c,
                    discriminant, square_root, root, inverse_two,
                    representative, &fast_workspace, tasks->ctx)) {
                slong factor_index;

                count = 0;
                fq_poly_roots(
                    fallback_factors,
                    tasks->evaluated + actual_index, 0, tasks->ctx);
                for (factor_index = 0;
                     factor_index < fallback_factors->num;
                     factor_index++) {
                    if (fq_poly_degree(
                            fallback_factors->poly + factor_index,
                            tasks->ctx) != 1)
                        continue;
                    if (count == tasks->root_stride)
                        break;
                    fq_poly_get_coeff(
                        c, fallback_factors->poly + factor_index,
                        0, tasks->ctx);
                    cubic_big_neg(
                        root, c, &fast_workspace, tasks->ctx);
                    fq_to_key(
                        point_roots + count++, root,
                        representative, tasks->ctx);
                }
            }
            tasks->root_counts[actual_index] =
                sort_unique_root_keys(point_roots, count);
            local_count++;
            if (root_profile_enabled)
                local_cpu_ns += thread_cpu_nanoseconds() - start;
        }
    }

    if (root_profile_enabled) {
        atomic_fetch_add_explicit(
            root_profile_counts + 2, local_count,
            memory_order_relaxed);
        atomic_fetch_add_explicit(
            root_profile_cpu_ns + 2, local_cpu_ns,
            memory_order_relaxed);
    }
    fast_fp2_workspace_clear(&fast_workspace);
    fmpz_poly_clear(representative);
    fq_clear(inverse_two, tasks->ctx);
    fq_clear(root, tasks->ctx);
    fq_clear(square_root, tasks->ctx);
    fq_clear(discriminant, tasks->ctx);
    fq_clear(c, tasks->ctx);
    fq_clear(b, tasks->ctx);
    fq_poly_factor_clear(fallback_factors, tasks->ctx);
}

static void split_rabin_public_big(
    fq_poly_t factor_a, fq_poly_t factor_b,
    const fq_poly_t polynomial, const fmpz_t halfq,
    fq_poly_t random_linear, fq_poly_t inverse_series,
    flint_rand_t randstate, const fq_ctx_t ctx)
{
    slong degree, split_degree;

    degree = fq_poly_degree(polynomial, ctx);
    if (degree <= 1)
        die("internal Rabin split received a linear polynomial");

    fq_poly_reverse(
        random_linear, polynomial, polynomial->length, ctx);
    fq_poly_inv_series_newton(
        inverse_series, random_linear, random_linear->length, ctx);

    for (;;) {
        fq_poly_fit_length(factor_a, 2, ctx);
        fq_rand(factor_a->coeffs + 0, randstate, ctx);
        fq_rand_not_zero(factor_a->coeffs + 1, randstate, ctx);
        factor_a->length = 2;

        fq_poly_powmod_fmpz_sliding_preinv(
            random_linear, factor_a, halfq, 0, polynomial,
            inverse_series, ctx);
        fq_poly_add_si(random_linear, random_linear, -1, ctx);
        fq_poly_gcd(factor_a, random_linear, polynomial, ctx);
        split_degree = fq_poly_degree(factor_a, ctx);
        if (split_degree > 0 && split_degree < degree)
            break;
    }

    fq_poly_div(factor_b, polynomial, factor_a, ctx);
    if (fq_poly_degree(factor_a, ctx) <
        fq_poly_degree(factor_b, ctx))
        fq_poly_swap(factor_a, factor_b, ctx);
}

static size_t fully_split_roots_big(
    jkey_t *roots, size_t capacity, const fq_poly_t polynomial,
    const fmpz_t halfq, flint_rand_t randstate,
    fq_poly_t factor_a, fq_poly_t factor_b, fq_poly_t temporary,
    fq_poly_t inverse_series, fq_poly_struct *stack, fq_t b, fq_t c,
    fq_t discriminant, fq_t square_root, fq_t root,
    const fq_t inverse_two, fmpz_poly_t representative,
    fast_fp2_workspace_t *fast_workspace, const fq_ctx_t ctx)
{
    size_t count = 0, stack_length = 1;

    fq_poly_make_monic(stack, polynomial, ctx);
    while (stack_length != 0) {
        fq_poly_struct *factor = stack + --stack_length;
        slong degree = fq_poly_degree(factor, ctx);

        if (degree <= 2) {
            if (!append_small_roots_big(
                    roots, &count, capacity, factor, b, c, discriminant,
                    square_root, root, inverse_two, representative,
                    fast_workspace, ctx))
                return SIZE_MAX;
            continue;
        }
        if (stack_length + 2 > ROOT_FACTOR_STACK_SIZE)
            return SIZE_MAX;

#if ISOGENY_USE_PRIVATE_RABIN
        _fq_poly_split_rabin(
            factor_a, factor_b, factor, halfq, temporary, inverse_series,
            randstate, ctx);
#else
        split_rabin_public_big(
            factor_a, factor_b, factor, halfq, temporary, inverse_series,
            randstate, ctx);
#endif
        fq_poly_swap(stack + stack_length, factor_a, ctx);
        stack_length++;
        fq_poly_swap(stack + stack_length, factor_b, ctx);
        stack_length++;
    }
    return count;
}

static void find_roots_task_big(slong task_index, void *argument)
{
    big_root_tasks_t *tasks = argument;
    fq_poly_factor_t roots;
    fq_poly_t factor_a, factor_b, temporary, inverse_series;
    fq_poly_struct factor_stack[ROOT_FACTOR_STACK_SIZE];
    fq_t b, c, discriminant, square_root, root, inverse_two;
    fmpz_t halfq;
    flint_rand_t randstate;
    fmpz_poly_t representative;
    fast_fp2_workspace_t fast_workspace;
    size_t first, last, i, stack_index;
    uint64_t *local_counts = NULL, *local_cpu_ns = NULL;
    uint64_t seed;

    if (root_profile_enabled) {
        local_counts = xcalloc(tasks->root_stride + 1,
                               sizeof(*local_counts));
        local_cpu_ns = xcalloc(tasks->root_stride + 1,
                               sizeof(*local_cpu_ns));
    }
    fq_poly_factor_init(roots, tasks->ctx);
    fq_poly_init(factor_a, tasks->ctx);
    fq_poly_init(factor_b, tasks->ctx);
    fq_poly_init(temporary, tasks->ctx);
    fq_poly_init(inverse_series, tasks->ctx);
    for (stack_index = 0; stack_index < ROOT_FACTOR_STACK_SIZE;
         stack_index++)
        fq_poly_init(factor_stack + stack_index, tasks->ctx);
    fq_init(b, tasks->ctx);
    fq_init(c, tasks->ctx);
    fq_init(discriminant, tasks->ctx);
    fq_init(square_root, tasks->ctx);
    fq_init(root, tasks->ctx);
    fq_init(inverse_two, tasks->ctx);
    fq_set_ui(inverse_two, 2, tasks->ctx);
    fq_inv(inverse_two, inverse_two, tasks->ctx);
    fmpz_init(halfq);
    fq_ctx_order(halfq, tasks->ctx);
    fmpz_sub_ui(halfq, halfq, 1);
    fmpz_fdiv_q_2exp(halfq, halfq, 1);
    ISOGENY_RAND_INIT(randstate);
    seed = mix64(UINT64_C(0x94d049bb133111eb) +
                 (uint64_t)(task_index + 1));
    ISOGENY_RAND_SET_SEED(
        randstate, (ulong)seed,
        (ulong)mix64(seed ^ UINT64_C(0xbf58476d1ce4e5b9)));
    fmpz_poly_init(representative);
    fast_fp2_workspace_init(&fast_workspace, tasks->ctx);

    for (;;) {
        first = atomic_fetch_add_explicit(
            &tasks->next_index, tasks->claim_size, memory_order_relaxed);
        if (first >= tasks->count)
            break;
        last = first + tasks->claim_size;
        if (last > tasks->count)
            last = tasks->count;
        for (i = first; i < last; i++) {
            size_t actual_index = tasks->indices == NULL
                                      ? i
                                      : tasks->indices[i];
            slong degree = fq_poly_degree(
                tasks->evaluated + actual_index, tasks->ctx);
            size_t count;
            uint64_t profile_start = 0;
            jkey_t *point_roots =
                tasks->root_keys + actual_index * tasks->root_stride;

            if ((tasks->skip_degree_two && degree == 2) ||
                (tasks->skip_degree_three && degree == 3))
                continue;
            if (root_profile_enabled)
                profile_start = thread_cpu_nanoseconds();
            count = fully_split_roots_big(
                point_roots, tasks->root_stride,
                tasks->evaluated + actual_index,
                halfq, randstate, factor_a, factor_b, temporary,
                inverse_series, factor_stack, b, c, discriminant,
                square_root, root, inverse_two, representative,
                &fast_workspace, tasks->ctx);
            if (count == SIZE_MAX) {
                slong factor_index;

                count = 0;
                fq_poly_roots(
                    roots, tasks->evaluated + actual_index, 0, tasks->ctx);
                for (factor_index = 0; factor_index < roots->num;
                     factor_index++) {
                    if (fq_poly_degree(
                            roots->poly + factor_index,
                            tasks->ctx) != 1)
                        continue;
                    if (count == tasks->root_stride)
                        break;
                    fq_poly_get_coeff(
                        c, roots->poly + factor_index, 0, tasks->ctx);
                    fq_neg(root, c, tasks->ctx);
                    fq_to_key(
                        point_roots + count++, root, representative,
                        tasks->ctx);
                }
            }
            qsort(point_roots, count, sizeof(*point_roots), compare_keys);
            if (count > 1) {
                size_t read_index, unique_count = 1;

                for (read_index = 1; read_index < count; read_index++) {
                    if (!key_equal(point_roots[read_index],
                                   point_roots[unique_count - 1])) {
                        if (unique_count != read_index)
                            key_set(point_roots + unique_count,
                                    point_roots + read_index);
                        unique_count++;
                    }
                }
                count = unique_count;
            }
            tasks->root_counts[actual_index] = count;
            if (root_profile_enabled && degree >= 0 &&
                (size_t)degree <= tasks->root_stride) {
                local_counts[degree]++;
                local_cpu_ns[degree] +=
                    thread_cpu_nanoseconds() - profile_start;
            }
        }
    }

    if (root_profile_enabled) {
        size_t degree;

        for (degree = 0;
             degree <= tasks->root_stride &&
             degree < ROOT_PROFILE_DEGREES;
             degree++) {
            if (local_counts[degree] == 0)
                continue;
            atomic_fetch_add_explicit(
                root_profile_counts + degree, local_counts[degree],
                memory_order_relaxed);
            atomic_fetch_add_explicit(
                root_profile_cpu_ns + degree, local_cpu_ns[degree],
                memory_order_relaxed);
        }
        free(local_cpu_ns);
        free(local_counts);
    }
    fast_fp2_workspace_clear(&fast_workspace);
    fmpz_poly_clear(representative);
    ISOGENY_RAND_CLEAR(randstate);
    fmpz_clear(halfq);
    fq_clear(inverse_two, tasks->ctx);
    fq_clear(root, tasks->ctx);
    fq_clear(square_root, tasks->ctx);
    fq_clear(discriminant, tasks->ctx);
    fq_clear(c, tasks->ctx);
    fq_clear(b, tasks->ctx);
    for (stack_index = ROOT_FACTOR_STACK_SIZE; stack_index-- > 0;)
        fq_poly_clear(factor_stack + stack_index, tasks->ctx);
    fq_poly_clear(inverse_series, tasks->ctx);
    fq_poly_clear(temporary, tasks->ctx);
    fq_poly_clear(factor_b, tasks->ctx);
    fq_poly_clear(factor_a, tasks->ctx);
    fq_poly_factor_clear(roots, tasks->ctx);
}

/*
 * A conventional batch assigns one complete specialization to one worker.
 * That is optimal once the batch contains many specializations, but leaves
 * most cores idle for the first high-prime stages, where a single polynomial
 * can have degree in the hundreds or thousands.  The cooperative path below
 * executes one Rabin-splitting level at a time.  Once the first factor has
 * split, its independent children fan out over the whole worker pool.
 */
typedef struct {
    fq_poly_struct polynomial;
    size_t source_slot;
} parallel_big_factor_t;

typedef struct {
    fq_poly_t factor_a, factor_b, temporary, inverse_series;
    fq_t b, c, discriminant, square_root, root, inverse_two;
    fmpz_t halfq;
    flint_rand_t randstate;
    fmpz_poly_t representative;
    fast_fp2_workspace_t fast;
} parallel_big_factor_workspace_t;

typedef struct {
    parallel_big_factor_workspace_t *workspaces;
    size_t count;
    const fq_ctx_struct *ctx;
} parallel_big_workspace_tasks_t;

typedef struct {
    parallel_big_factor_t *current;
    parallel_big_factor_t *children;
    unsigned char *child_counts;
    size_t count;
    atomic_size_t next_index;
    size_t claim_size;
    parallel_big_factor_workspace_t *workspaces;
    atomic_size_t *root_writes;
    atomic_uint_fast64_t *profile_cpu_ns;
    const size_t *actual_indices;
    jkey_t *root_keys;
    size_t root_stride;
    atomic_int failed;
    const fq_ctx_struct *ctx;
} parallel_big_factor_level_tasks_t;

static void parallel_big_workspace_init_task(slong task_index,
                                             void *argument)
{
    parallel_big_workspace_tasks_t *tasks = argument;
    parallel_big_factor_workspace_t *workspace =
        tasks->workspaces + (size_t)task_index;
    uint64_t seed;

    fq_poly_init(workspace->factor_a, tasks->ctx);
    fq_poly_init(workspace->factor_b, tasks->ctx);
    fq_poly_init(workspace->temporary, tasks->ctx);
    fq_poly_init(workspace->inverse_series, tasks->ctx);
    fq_init(workspace->b, tasks->ctx);
    fq_init(workspace->c, tasks->ctx);
    fq_init(workspace->discriminant, tasks->ctx);
    fq_init(workspace->square_root, tasks->ctx);
    fq_init(workspace->root, tasks->ctx);
    fq_init(workspace->inverse_two, tasks->ctx);
    fq_set_ui(workspace->inverse_two, 2, tasks->ctx);
    fq_inv(workspace->inverse_two, workspace->inverse_two, tasks->ctx);
    fmpz_init(workspace->halfq);
    fq_ctx_order(workspace->halfq, tasks->ctx);
    fmpz_sub_ui(workspace->halfq, workspace->halfq, 1);
    fmpz_fdiv_q_2exp(workspace->halfq, workspace->halfq, 1);
    ISOGENY_RAND_INIT(workspace->randstate);
    seed = mix64(UINT64_C(0x243f6a8885a308d3) +
                 (uint64_t)(task_index + 1));
    ISOGENY_RAND_SET_SEED(
        workspace->randstate, (ulong)seed,
        (ulong)mix64(seed ^ UINT64_C(0x13198a2e03707344)));
    fmpz_poly_init(workspace->representative);
    fast_fp2_workspace_init(&workspace->fast, tasks->ctx);
}

static void parallel_big_workspace_clear_task(slong task_index,
                                              void *argument)
{
    parallel_big_workspace_tasks_t *tasks = argument;
    parallel_big_factor_workspace_t *workspace =
        tasks->workspaces + (size_t)task_index;

    fast_fp2_workspace_clear(&workspace->fast);
    fmpz_poly_clear(workspace->representative);
    ISOGENY_RAND_CLEAR(workspace->randstate);
    fmpz_clear(workspace->halfq);
    fq_clear(workspace->inverse_two, tasks->ctx);
    fq_clear(workspace->root, tasks->ctx);
    fq_clear(workspace->square_root, tasks->ctx);
    fq_clear(workspace->discriminant, tasks->ctx);
    fq_clear(workspace->c, tasks->ctx);
    fq_clear(workspace->b, tasks->ctx);
    fq_poly_clear(workspace->inverse_series, tasks->ctx);
    fq_poly_clear(workspace->temporary, tasks->ctx);
    fq_poly_clear(workspace->factor_b, tasks->ctx);
    fq_poly_clear(workspace->factor_a, tasks->ctx);
}

static void parallel_big_factor_level_task(slong task_index,
                                           void *argument)
{
    parallel_big_factor_level_tasks_t *tasks = argument;
    parallel_big_factor_workspace_t *workspace =
        tasks->workspaces + (size_t)task_index;
    size_t first;

    for (;;) {
        size_t last, factor_index;

        first = atomic_fetch_add_explicit(
            &tasks->next_index, tasks->claim_size,
            memory_order_relaxed);
        if (first >= tasks->count)
            break;
        last = first + tasks->claim_size;
        if (last > tasks->count)
            last = tasks->count;
        for (factor_index = first; factor_index < last;
             factor_index++) {
            parallel_big_factor_t *factor =
                tasks->current + factor_index;
            size_t source_slot = factor->source_slot;
            slong degree = fq_poly_degree(
                &factor->polynomial, tasks->ctx);
            uint64_t started = root_profile_enabled
                                   ? thread_cpu_nanoseconds() : 0;

            tasks->child_counts[factor_index] = 0;
            if (degree <= 2) {
                jkey_t local_roots[2];
                size_t local_count = 0, root_index;

                key_init(local_roots + 0);
                key_init(local_roots + 1);
                if (!append_small_roots_big(
                        local_roots, &local_count, 2,
                        &factor->polynomial,
                        workspace->b, workspace->c,
                        workspace->discriminant,
                        workspace->square_root, workspace->root,
                        workspace->inverse_two,
                        workspace->representative,
                        &workspace->fast, tasks->ctx)) {
                    atomic_store_explicit(
                        &tasks->failed, 1, memory_order_relaxed);
                    continue;
                }
                for (root_index = 0; root_index < local_count;
                     root_index++) {
                    size_t position = atomic_fetch_add_explicit(
                        tasks->root_writes + source_slot, 1,
                        memory_order_relaxed);
                    size_t actual_index =
                        tasks->actual_indices[source_slot];

                    if (position >= tasks->root_stride) {
                        atomic_store_explicit(
                            &tasks->failed, 1,
                            memory_order_relaxed);
                        continue;
                    }
                    key_set(
                        tasks->root_keys +
                            actual_index * tasks->root_stride +
                            position,
                        local_roots + root_index);
                }
                key_clear(local_roots + 1);
                key_clear(local_roots + 0);
            } else {
#if ISOGENY_USE_PRIVATE_RABIN
                _fq_poly_split_rabin(
                    workspace->factor_a, workspace->factor_b,
                    &factor->polynomial, workspace->halfq,
                    workspace->temporary, workspace->inverse_series,
                    workspace->randstate, tasks->ctx);
#else
                split_rabin_public_big(
                    workspace->factor_a, workspace->factor_b,
                    &factor->polynomial, workspace->halfq,
                    workspace->temporary, workspace->inverse_series,
                    workspace->randstate, tasks->ctx);
#endif
                fq_poly_swap(
                    &tasks->children[2 * factor_index].polynomial,
                    workspace->factor_a, tasks->ctx);
                tasks->children[2 * factor_index].source_slot =
                    source_slot;
                fq_poly_swap(
                    &tasks->children[2 * factor_index + 1].polynomial,
                    workspace->factor_b, tasks->ctx);
                tasks->children[2 * factor_index + 1].source_slot =
                    source_slot;
                tasks->child_counts[factor_index] = 2;
            }
            if (root_profile_enabled) {
                atomic_fetch_add_explicit(
                    tasks->profile_cpu_ns + source_slot,
                    thread_cpu_nanoseconds() - started,
                    memory_order_relaxed);
            }
        }
    }
}

static int should_use_parallel_factor_tree_big(
    const fq_poly_struct *evaluated, const size_t *indices,
    size_t count, const fq_ctx_t ctx, unsigned n_threads)
{
    size_t i;
    slong maximum_degree = 0;

    if (!parallel_factor_tree_enabled || n_threads < 4 ||
        count == 0 || count >= (size_t)n_threads)
        return 0;
    for (i = 0; i < count; i++) {
        slong degree = fq_poly_degree(
            evaluated + indices[i], ctx);

        if (degree > maximum_degree)
            maximum_degree = degree;
    }
    return maximum_degree >= PARALLEL_FACTOR_MIN_DEGREE;
}

static void run_parallel_factor_tree_big(
    const fq_poly_struct *evaluated, const size_t *indices,
    size_t count, size_t root_stride, jkey_t *root_keys,
    size_t *root_counts, const fq_ctx_t ctx, unsigned n_threads)
{
    parallel_big_factor_t *current;
    parallel_big_factor_workspace_t *workspaces;
    parallel_big_workspace_tasks_t workspace_tasks;
    atomic_size_t *root_writes;
    atomic_uint_fast64_t *profile_cpu_ns;
    slong *initial_degrees;
    size_t workspace_count = (size_t)n_threads;
    size_t current_count = count;
    size_t i;

    current = xrealloc(NULL, count, sizeof(*current));
    root_writes = xrealloc(NULL, count, sizeof(*root_writes));
    profile_cpu_ns =
        xrealloc(NULL, count, sizeof(*profile_cpu_ns));
    initial_degrees =
        xrealloc(NULL, count, sizeof(*initial_degrees));
    for (i = 0; i < count; i++) {
        fq_poly_init(&current[i].polynomial, ctx);
        fq_poly_make_monic(
            &current[i].polynomial, evaluated + indices[i], ctx);
        current[i].source_slot = i;
        atomic_init(root_writes + i, 0);
        atomic_init(profile_cpu_ns + i, 0);
        initial_degrees[i] =
            fq_poly_degree(evaluated + indices[i], ctx);
    }

    workspaces =
        xrealloc(NULL, workspace_count, sizeof(*workspaces));
    workspace_tasks.workspaces = workspaces;
    workspace_tasks.count = workspace_count;
    workspace_tasks.ctx = ctx;
    flint_parallel_do(
        parallel_big_workspace_init_task, &workspace_tasks,
        (slong)workspace_count, (int)n_threads,
        FLINT_PARALLEL_UNIFORM);

    while (current_count != 0) {
        parallel_big_factor_t *children, *compacted;
        unsigned char *child_counts;
        parallel_big_factor_level_tasks_t level_tasks;
        size_t child_capacity, next_count = 0, write = 0;
        size_t task_count;

        if (current_count > SIZE_MAX / 2)
            die("parallel factor tree size overflow");
        child_capacity = 2 * current_count;
        children =
            xrealloc(NULL, child_capacity, sizeof(*children));
        child_counts = xcalloc(current_count, sizeof(*child_counts));
        for (i = 0; i < child_capacity; i++) {
            fq_poly_init(&children[i].polynomial, ctx);
            children[i].source_slot = 0;
        }

        task_count = current_count < workspace_count
                         ? current_count : workspace_count;
        level_tasks.current = current;
        level_tasks.children = children;
        level_tasks.child_counts = child_counts;
        level_tasks.count = current_count;
        atomic_init(&level_tasks.next_index, 0);
        level_tasks.claim_size =
            adaptive_root_claim_size(current_count, task_count);
        level_tasks.workspaces = workspaces;
        level_tasks.root_writes = root_writes;
        level_tasks.profile_cpu_ns = profile_cpu_ns;
        level_tasks.actual_indices = indices;
        level_tasks.root_keys = root_keys;
        level_tasks.root_stride = root_stride;
        atomic_init(&level_tasks.failed, 0);
        level_tasks.ctx = ctx;
        flint_parallel_do(
            parallel_big_factor_level_task, &level_tasks,
            (slong)task_count, (int)n_threads,
            FLINT_PARALLEL_UNIFORM);
        if (atomic_load_explicit(
                &level_tasks.failed, memory_order_relaxed))
            die("parallel factor tree failed to extract split roots");

        for (i = 0; i < current_count; i++)
            next_count += child_counts[i];
        compacted =
            xrealloc(NULL, next_count, sizeof(*compacted));
        for (i = 0; i < next_count; i++) {
            fq_poly_init(&compacted[i].polynomial, ctx);
            compacted[i].source_slot = 0;
        }
        for (i = 0; i < current_count; i++) {
            size_t child;

            for (child = 0; child < child_counts[i]; child++) {
                size_t source = 2 * i + child;

                fq_poly_swap(
                    &compacted[write].polynomial,
                    &children[source].polynomial, ctx);
                compacted[write].source_slot =
                    children[source].source_slot;
                write++;
            }
        }
        if (write != next_count)
            die("internal error compacting parallel factor level");

        for (i = 0; i < child_capacity; i++)
            fq_poly_clear(&children[i].polynomial, ctx);
        for (i = 0; i < current_count; i++)
            fq_poly_clear(&current[i].polynomial, ctx);
        free(child_counts);
        free(children);
        free(current);
        current = compacted;
        current_count = next_count;
    }
    free(current);

    flint_parallel_do(
        parallel_big_workspace_clear_task, &workspace_tasks,
        (slong)workspace_count, (int)n_threads,
        FLINT_PARALLEL_UNIFORM);
    free(workspaces);

    for (i = 0; i < count; i++) {
        size_t actual_index = indices[i];
        size_t root_count = atomic_load_explicit(
            root_writes + i, memory_order_relaxed);

        root_counts[actual_index] = sort_unique_root_keys(
            root_keys + actual_index * root_stride, root_count);
        if (root_profile_enabled &&
            initial_degrees[i] >= 0 &&
            (size_t)initial_degrees[i] < ROOT_PROFILE_DEGREES) {
            atomic_fetch_add_explicit(
                root_profile_counts + initial_degrees[i], 1,
                memory_order_relaxed);
            atomic_fetch_add_explicit(
                root_profile_cpu_ns + initial_degrees[i],
                atomic_load_explicit(
                    profile_cpu_ns + i, memory_order_relaxed),
                memory_order_relaxed);
        }
    }
    free(initial_degrees);
    free(profile_cpu_ns);
    free(root_writes);
}

static void run_word_root_pass(
    const fq_nmod_poly_struct *evaluated, const size_t *indices,
    size_t logical_count, int skip_degree_three,
    size_t root_stride, jkey_t *root_keys, size_t *root_counts,
    const fq_nmod_ctx_t ctx, unsigned n_threads)
{
    word_root_tasks_t tasks;
    size_t task_count;

    if (logical_count == 0)
        return;
    task_count = logical_count;
    if (task_count > (size_t)n_threads)
        task_count = (size_t)n_threads;
    tasks.evaluated = evaluated;
    tasks.indices = indices;
    tasks.count = logical_count;
    atomic_init(&tasks.next_index, 0);
    tasks.claim_size =
        adaptive_root_claim_size(logical_count, task_count);
    tasks.skip_degree_three = skip_degree_three;
    tasks.root_stride = root_stride;
    tasks.root_keys = root_keys;
    tasks.root_counts = root_counts;
    tasks.ctx = ctx;
    flint_parallel_do(find_roots_task_word, &tasks,
                      (slong)task_count, (int)n_threads,
                      FLINT_PARALLEL_UNIFORM);
}

static void run_big_root_pass(
    const fq_poly_struct *evaluated, const size_t *indices,
    size_t logical_count, int skip_degree_two,
    int skip_degree_three,
    size_t root_stride, jkey_t *root_keys, size_t *root_counts,
    const fq_ctx_t ctx, unsigned n_threads)
{
    big_root_tasks_t tasks;
    size_t task_count;

    if (logical_count == 0)
        return;
    task_count = logical_count;
    if (task_count > (size_t)n_threads)
        task_count = (size_t)n_threads;
    tasks.evaluated = evaluated;
    tasks.indices = indices;
    tasks.count = logical_count;
    atomic_init(&tasks.next_index, 0);
    tasks.claim_size =
        adaptive_root_claim_size(logical_count, task_count);
    tasks.skip_degree_two = skip_degree_two;
    tasks.skip_degree_three = skip_degree_three;
    tasks.root_stride = root_stride;
    tasks.root_keys = root_keys;
    tasks.root_counts = root_counts;
    tasks.ctx = ctx;
    flint_parallel_do(find_roots_task_big, &tasks,
                      (slong)task_count, (int)n_threads,
                      FLINT_PARALLEL_UNIFORM);
}

static void run_quadratic_big_root_pass(
    const fq_poly_struct *evaluated, const size_t *indices,
    size_t count, size_t root_stride, jkey_t *root_keys,
    size_t *root_counts, const fq_ctx_t ctx, unsigned n_threads)
{
    quadratic_big_root_tasks_t tasks;
    size_t task_count, useful_tasks;

    if (count == 0)
        return;
    useful_tasks =
        (count + MIN_QUADRATIC_ROOTS_PER_TASK - 1) /
        MIN_QUADRATIC_ROOTS_PER_TASK;
    if (useful_tasks == 0)
        useful_tasks = 1;
    task_count = useful_tasks;
    if (task_count > (size_t)n_threads)
        task_count = (size_t)n_threads;
    tasks.evaluated = evaluated;
    tasks.indices = indices;
    tasks.count = count;
    atomic_init(&tasks.next_index, 0);
    tasks.claim_size =
        adaptive_root_claim_size(count, task_count);
    tasks.root_stride = root_stride;
    tasks.root_keys = root_keys;
    tasks.root_counts = root_counts;
    tasks.ctx = ctx;
    flint_parallel_do(
        find_quadratic_roots_task_big, &tasks,
        (slong)task_count, (int)n_threads, FLINT_PARALLEL_UNIFORM);
}

static size_t adaptive_cubic_lane_width(size_t count,
                                        unsigned n_threads)
{
    size_t target_groups, width;

    if (count == 0 || n_threads <= 1)
        return CUBIC_BATCH_LANES;
    target_groups =
        (size_t)n_threads * EVALUATION_TASKS_PER_WORKER;
    width = count / target_groups;
    if (count % target_groups != 0)
        width++;
    if (width == 0)
        width = 1;
    if (width > CUBIC_BATCH_LANES)
        width = CUBIC_BATCH_LANES;
    return width;
}

static void run_batched_cubic_pass(
    const fq_nmod_poly_struct *evaluated, const size_t *indices,
    size_t cubic_count, size_t root_stride, jkey_t *root_keys,
    size_t *root_counts, const fq_nmod_ctx_t ctx, unsigned n_threads)
{
    cubic_root_tasks_t tasks;
    size_t task_count, lane_width;

    if (cubic_count == 0)
        return;
    lane_width = adaptive_cubic_lane_width(cubic_count, n_threads);
    task_count =
        (cubic_count + lane_width - 1) / lane_width;
    if (task_count > (size_t)n_threads)
        task_count = (size_t)n_threads;
    tasks.evaluated = evaluated;
    tasks.indices = indices;
    tasks.count = cubic_count;
    atomic_init(&tasks.next_index, 0);
    tasks.lane_width = lane_width;
    tasks.root_keys = root_keys;
    tasks.root_counts = root_counts;
    tasks.root_stride = root_stride;
    tasks.ctx = ctx;
    flint_parallel_do(find_cubic_roots_task_word, &tasks,
                      (slong)task_count, (int)n_threads,
                      FLINT_PARALLEL_UNIFORM);
}

static void run_batched_cubic_pass_big(
    const fq_poly_struct *evaluated, const size_t *indices,
    size_t cubic_count, size_t root_stride, jkey_t *root_keys,
    size_t *root_counts, const fq_ctx_t ctx, unsigned n_threads)
{
    cubic_big_root_tasks_t tasks;
    size_t task_count, lane_width;

    if (cubic_count == 0)
        return;
    lane_width = adaptive_cubic_lane_width(cubic_count, n_threads);
    task_count =
        (cubic_count + lane_width - 1) / lane_width;
    if (task_count > (size_t)n_threads)
        task_count = (size_t)n_threads;
    tasks.evaluated = evaluated;
    tasks.indices = indices;
    tasks.count = cubic_count;
    atomic_init(&tasks.next_index, 0);
    tasks.lane_width = lane_width;
    tasks.root_keys = root_keys;
    tasks.root_counts = root_counts;
    tasks.root_stride = root_stride;
    tasks.ctx = ctx;
    flint_parallel_do(find_cubic_roots_task_big, &tasks,
                      (slong)task_count, (int)n_threads,
                      FLINT_PARALLEL_UNIFORM);
}

/*
 * Run the selected word-backend root engine. Compare mode benchmarks the two
 * cubic paths in isolation on exactly the same specialization objects, checks
 * every sorted root set, and uses the FLINT result for graph construction.
 */
static void extract_roots_parallel_word(
    const fq_nmod_poly_struct *evaluated, size_t count,
    size_t root_stride, jkey_t *root_keys, size_t *root_counts,
    const fq_nmod_ctx_t ctx, unsigned n_threads)
{
    size_t *cubic_indices = NULL;
    size_t cubic_count = 0, i;

    if (cubic_kernel == CUBIC_KERNEL_FLINT) {
        run_word_root_pass(
            evaluated, NULL, count, 0, root_stride,
            root_keys, root_counts, ctx, n_threads);
        return;
    }

    cubic_indices = xrealloc(NULL, count, sizeof(*cubic_indices));
    for (i = 0; i < count; i++)
        if (fq_nmod_poly_degree(evaluated + i, ctx) == 3)
            cubic_indices[cubic_count++] = i;

    if (cubic_kernel == CUBIC_KERNEL_COMPARE && cubic_count != 0) {
        jkey_t *batched_keys =
            xrealloc(NULL, count * root_stride, sizeof(*batched_keys));
        size_t *batched_counts =
            xcalloc(count, sizeof(*batched_counts));
        double start, end;

        key_array_init(batched_keys, count * root_stride);
        start = wall_seconds();
        run_word_root_pass(
            evaluated, cubic_indices, cubic_count, 0,
            root_stride, root_keys, root_counts, ctx, n_threads);
        end = wall_seconds();
        atomic_fetch_add_explicit(
            &cubic_control_wall_ns, elapsed_nanoseconds(start, end),
            memory_order_relaxed);

        start = wall_seconds();
        run_batched_cubic_pass(
            evaluated, cubic_indices, cubic_count, root_stride,
            batched_keys, batched_counts, ctx, n_threads);
        end = wall_seconds();
        atomic_fetch_add_explicit(
            &cubic_batched_wall_ns, elapsed_nanoseconds(start, end),
            memory_order_relaxed);

        for (i = 0; i < cubic_count; i++) {
            size_t actual_index = cubic_indices[i];
            size_t root_index;

            if (root_counts[actual_index] !=
                batched_counts[actual_index])
                die("batched cubic kernel returned a different root count");
            for (root_index = 0;
                 root_index < root_counts[actual_index]; root_index++) {
                if (!key_equal(
                        root_keys[actual_index * root_stride + root_index],
                        batched_keys[
                            actual_index * root_stride + root_index]))
                    die("batched cubic kernel returned a different root");
            }
        }
        key_array_clear(batched_keys, count * root_stride);
        free(batched_counts);
        free(batched_keys);
    } else if (cubic_count != 0) {
        double start = wall_seconds();

        run_batched_cubic_pass(
            evaluated, cubic_indices, cubic_count, root_stride,
            root_keys, root_counts, ctx, n_threads);
        atomic_fetch_add_explicit(
            &cubic_batched_wall_ns,
            elapsed_nanoseconds(start, wall_seconds()),
            memory_order_relaxed);
    }

    /* Degree-three slots were filled above; process every other degree. */
    run_word_root_pass(
        evaluated, NULL, count, 1, root_stride,
        root_keys, root_counts, ctx, n_threads);
    free(cubic_indices);
}

static void extract_roots_parallel_big(
    const fq_poly_struct *evaluated, size_t count,
    size_t root_stride, jkey_t *root_keys, size_t *root_counts,
    const fq_ctx_t ctx, unsigned n_threads)
{
    size_t *indices = NULL;
    size_t *quadratic_indices, *cubic_indices, *generic_indices;
    size_t quadratic_count = 0, cubic_count = 0, generic_count = 0;
    size_t quadratic_write, cubic_write, generic_write, i;

    for (i = 0; i < count; i++) {
        slong degree = fq_poly_degree(evaluated + i, ctx);

        if (quadratic_basis_is_minus_one && degree == 2)
            quadratic_count++;
        else if (cubic_kernel != CUBIC_KERNEL_FLINT && degree == 3)
            cubic_count++;
        else
            generic_count++;
    }
    indices = xrealloc(NULL, count, sizeof(*indices));
    quadratic_indices = indices;
    cubic_indices = indices + quadratic_count;
    generic_indices = cubic_indices + cubic_count;
    quadratic_write = cubic_write = generic_write = 0;
    for (i = 0; i < count; i++) {
        slong degree = fq_poly_degree(evaluated + i, ctx);

        if (quadratic_basis_is_minus_one && degree == 2)
            quadratic_indices[quadratic_write++] = i;
        else if (cubic_kernel != CUBIC_KERNEL_FLINT && degree == 3)
            cubic_indices[cubic_write++] = i;
        else
            generic_indices[generic_write++] = i;
    }

    run_quadratic_big_root_pass(
        evaluated, quadratic_indices, quadratic_count, root_stride,
        root_keys, root_counts, ctx, n_threads);

    if (cubic_kernel == CUBIC_KERNEL_COMPARE && cubic_count != 0) {
        jkey_t *batched_keys =
            xrealloc(NULL, count * root_stride, sizeof(*batched_keys));
        size_t *batched_counts =
            xcalloc(count, sizeof(*batched_counts));
        double start, end;

        key_array_init(batched_keys, count * root_stride);
        start = wall_seconds();
        run_big_root_pass(
            evaluated, cubic_indices, cubic_count, 0, 0,
            root_stride, root_keys, root_counts, ctx, n_threads);
        end = wall_seconds();
        atomic_fetch_add_explicit(
            &cubic_control_wall_ns, elapsed_nanoseconds(start, end),
            memory_order_relaxed);

        start = wall_seconds();
        run_batched_cubic_pass_big(
            evaluated, cubic_indices, cubic_count, root_stride,
            batched_keys, batched_counts, ctx, n_threads);
        end = wall_seconds();
        atomic_fetch_add_explicit(
            &cubic_batched_wall_ns, elapsed_nanoseconds(start, end),
            memory_order_relaxed);

        for (i = 0; i < cubic_count; i++) {
            size_t actual_index = cubic_indices[i];
            size_t root_index;

            if (root_counts[actual_index] !=
                batched_counts[actual_index])
                die("batched cubic kernel returned a different root count");
            for (root_index = 0;
                 root_index < root_counts[actual_index]; root_index++) {
                if (!key_equal(
                        root_keys[actual_index * root_stride + root_index],
                        batched_keys[
                            actual_index * root_stride + root_index]))
                    die("batched cubic kernel returned a different root");
            }
        }
        key_array_clear(batched_keys, count * root_stride);
        free(batched_counts);
        free(batched_keys);
    } else if (cubic_count != 0) {
        double start = wall_seconds();

        run_batched_cubic_pass_big(
            evaluated, cubic_indices, cubic_count, root_stride,
            root_keys, root_counts, ctx, n_threads);
        atomic_fetch_add_explicit(
            &cubic_batched_wall_ns,
            elapsed_nanoseconds(start, wall_seconds()),
            memory_order_relaxed);
    }

    if (should_use_parallel_factor_tree_big(
            evaluated, generic_indices, generic_count,
            ctx, n_threads)) {
        run_parallel_factor_tree_big(
            evaluated, generic_indices, generic_count,
            root_stride, root_keys, root_counts, ctx, n_threads);
    } else {
        run_big_root_pass(
            evaluated, generic_indices, generic_count, 0, 0,
            root_stride, root_keys, root_counts, ctx, n_threads);
    }
    free(indices);
}

/*
 * Synthetic division by Y - known_root. The quotient and remainder scratch
 * objects are reused across the whole batch, avoiding a FLINT allocation for
 * every known dual edge.
 */
static int divide_out_known_root_word(
    fq_nmod_poly_t polynomial, fq_nmod_poly_t quotient,
    const fq_nmod_t known_root, fq_nmod_t remainder,
    const fq_nmod_ctx_t ctx)
{
    slong degree = fq_nmod_poly_degree(polynomial, ctx);
    slong k;

    if (degree < 1)
        return 0;
    fq_nmod_poly_fit_length(quotient, degree, ctx);
    fq_nmod_set(quotient->coeffs + degree - 1,
                polynomial->coeffs + degree, ctx);
    for (k = degree - 1; k > 0; k--) {
        fq_nmod_mul(quotient->coeffs + k - 1, quotient->coeffs + k,
                    known_root, ctx);
        fq_nmod_add(quotient->coeffs + k - 1,
                    quotient->coeffs + k - 1,
                    polynomial->coeffs + k, ctx);
    }
    _fq_nmod_poly_set_length(quotient, degree, ctx);
    fq_nmod_mul(remainder, quotient->coeffs, known_root, ctx);
    fq_nmod_add(remainder, remainder, polynomial->coeffs, ctx);
    if (!fq_nmod_is_zero(remainder, ctx))
        return 0;
    fq_nmod_poly_swap(polynomial, quotient, ctx);
    return 1;
}

static int divide_out_known_root_big(
    fq_poly_t polynomial, fq_poly_t quotient, const fq_t known_root,
    fq_t remainder, const fq_ctx_t ctx)
{
    slong degree = fq_poly_degree(polynomial, ctx);
    slong k;

    if (degree < 1)
        return 0;
    fq_poly_fit_length(quotient, degree, ctx);
    fq_set(quotient->coeffs + degree - 1,
           polynomial->coeffs + degree, ctx);
    for (k = degree - 1; k > 0; k--) {
        fq_mul(quotient->coeffs + k - 1, quotient->coeffs + k,
               known_root, ctx);
        fq_add(quotient->coeffs + k - 1,
               quotient->coeffs + k - 1,
               polynomial->coeffs + k, ctx);
    }
    _fq_poly_set_length(quotient, degree, ctx);
    fq_mul(remainder, quotient->coeffs, known_root, ctx);
    fq_add(remainder, remainder, polynomial->coeffs, ctx);
    if (!fq_is_zero(remainder, ctx))
        return 0;
    fq_poly_swap(polynomial, quotient, ctx);
    return 1;
}

typedef struct {
    fq_nmod_poly_struct *evaluated;
    const size_t *indices;
    const jkey_t *known_roots;
    size_t count;
    atomic_size_t next_index;
    size_t claim_size;
    atomic_int failed;
    const fq_nmod_ctx_struct *ctx;
} word_dual_division_tasks_t;

typedef struct {
    fq_poly_struct *evaluated;
    const size_t *indices;
    const jkey_t *known_roots;
    size_t count;
    atomic_size_t next_index;
    size_t claim_size;
    atomic_int failed;
    const fq_ctx_struct *ctx;
} big_dual_division_tasks_t;

static void divide_dual_roots_task_word(slong task_index, void *argument)
{
    word_dual_division_tasks_t *tasks = argument;
    fq_nmod_poly_t quotient;
    fq_nmod_t known_root, remainder;
    nmod_poly_t representative;
    size_t first;

    (void)task_index;
    fq_nmod_poly_init(quotient, tasks->ctx);
    fq_nmod_init(known_root, tasks->ctx);
    fq_nmod_init(remainder, tasks->ctx);
    nmod_poly_init(representative, field_characteristic_word);
    for (;;) {
        size_t last, i;

        first = atomic_fetch_add_explicit(
            &tasks->next_index, tasks->claim_size,
            memory_order_relaxed);
        if (first >= tasks->count)
            break;
        last = first + tasks->claim_size;
        if (last > tasks->count)
            last = tasks->count;
        for (i = first; i < last; i++) {
            key_to_fq_nmod(
                known_root, tasks->known_roots[i],
                representative, tasks->ctx);
            if (!divide_out_known_root_word(
                    tasks->evaluated + tasks->indices[i], quotient,
                    known_root, remainder, tasks->ctx))
                atomic_store_explicit(
                    &tasks->failed, 1, memory_order_relaxed);
        }
    }
    nmod_poly_clear(representative);
    fq_nmod_clear(remainder, tasks->ctx);
    fq_nmod_clear(known_root, tasks->ctx);
    fq_nmod_poly_clear(quotient, tasks->ctx);
}

static void divide_dual_roots_task_big(slong task_index, void *argument)
{
    big_dual_division_tasks_t *tasks = argument;
    fq_poly_t quotient;
    fq_t known_root, remainder;
    fmpz_poly_t representative;
    size_t first;

    (void)task_index;
    fq_poly_init(quotient, tasks->ctx);
    fq_init(known_root, tasks->ctx);
    fq_init(remainder, tasks->ctx);
    fmpz_poly_init(representative);
    for (;;) {
        size_t last, i;

        first = atomic_fetch_add_explicit(
            &tasks->next_index, tasks->claim_size,
            memory_order_relaxed);
        if (first >= tasks->count)
            break;
        last = first + tasks->claim_size;
        if (last > tasks->count)
            last = tasks->count;
        for (i = first; i < last; i++) {
            key_to_fq(
                known_root, tasks->known_roots[i],
                representative, tasks->ctx);
            if (!divide_out_known_root_big(
                    tasks->evaluated + tasks->indices[i], quotient,
                    known_root, remainder, tasks->ctx))
                atomic_store_explicit(
                    &tasks->failed, 1, memory_order_relaxed);
        }
    }
    fmpz_poly_clear(representative);
    fq_clear(remainder, tasks->ctx);
    fq_clear(known_root, tasks->ctx);
    fq_poly_clear(quotient, tasks->ctx);
}

static void divide_dual_roots_parallel_word(
    fq_nmod_poly_struct *evaluated, const size_t *indices,
    const jkey_t *known_roots, size_t count,
    const fq_nmod_ctx_t ctx, unsigned n_threads)
{
    word_dual_division_tasks_t tasks;
    size_t task_count;

    if (count == 0)
        return;
    task_count = count < (size_t)n_threads
                     ? count : (size_t)n_threads;
    tasks.evaluated = evaluated;
    tasks.indices = indices;
    tasks.known_roots = known_roots;
    tasks.count = count;
    atomic_init(&tasks.next_index, 0);
    tasks.claim_size =
        adaptive_root_claim_size(count, task_count);
    atomic_init(&tasks.failed, 0);
    tasks.ctx = ctx;
    flint_parallel_do(
        divide_dual_roots_task_word, &tasks, (slong)task_count,
        (int)n_threads, FLINT_PARALLEL_UNIFORM);
    if (atomic_load_explicit(&tasks.failed, memory_order_relaxed))
        die("known dual isogeny is not a root of the specialization");
}

static void divide_dual_roots_parallel_big(
    fq_poly_struct *evaluated, const size_t *indices,
    const jkey_t *known_roots, size_t count,
    const fq_ctx_t ctx, unsigned n_threads)
{
    big_dual_division_tasks_t tasks;
    size_t task_count;

    if (count == 0)
        return;
    task_count = count < (size_t)n_threads
                     ? count : (size_t)n_threads;
    tasks.evaluated = evaluated;
    tasks.indices = indices;
    tasks.known_roots = known_roots;
    tasks.count = count;
    atomic_init(&tasks.next_index, 0);
    tasks.claim_size =
        adaptive_root_claim_size(count, task_count);
    atomic_init(&tasks.failed, 0);
    tasks.ctx = ctx;
    flint_parallel_do(
        divide_dual_roots_task_big, &tasks, (slong)task_count,
        (int)n_threads, FLINT_PARALLEL_UNIFORM);
    if (atomic_load_explicit(&tasks.failed, memory_order_relaxed))
        die("known dual isogeny is not a root of the specialization");
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

typedef struct {
    jkey_t *edges;
    size_t count;
} root_edge_cache_entry_t;

typedef struct {
    root_edge_cache_entry_t *entries;
    jkey_t point;
    const modular_poly_t *polynomials;
    size_t count;
    atomic_size_t next_index;
    const field_context_t *field;
} root_edge_precompute_tasks_t;

static void precompute_root_edges_task(slong task_index, void *argument)
{
    root_edge_precompute_tasks_t *tasks = argument;
    size_t claimed;

    (void)task_index;
    for (;;) {
        size_t polynomial_index;
        root_edge_cache_entry_t *entry;

        claimed = atomic_fetch_add_explicit(
            &tasks->next_index, 1, memory_order_relaxed);
        if (claimed >= tasks->count)
            break;
        /* Start the expensive, high-degree specializations first. */
        polynomial_index = tasks->count - 1 - claimed;
        entry = tasks->entries + polynomial_index;
        entry->count = modular_poly_edges_at(
            &entry->edges, tasks->point,
            tasks->polynomials + polynomial_index, tasks->field);
    }
}

static root_edge_cache_entry_t *precompute_root_edges(
    jkey_t point, const modular_poly_t *polynomials,
    size_t count, const field_context_t *field,
    unsigned n_threads)
{
    root_edge_cache_entry_t *entries;
    root_edge_precompute_tasks_t tasks;
    size_t task_count;

    entries = xcalloc(count, sizeof(*entries));
    if (count == 0)
        return entries;
    task_count = count < (size_t)n_threads
                     ? count : (size_t)n_threads;
    tasks.entries = entries;
    tasks.point = point;
    tasks.polynomials = polynomials;
    tasks.count = count;
    atomic_init(&tasks.next_index, 0);
    tasks.field = field;
    flint_parallel_do(
        precompute_root_edges_task, &tasks, (slong)task_count,
        (int)n_threads, FLINT_PARALLEL_UNIFORM);
    return entries;
}

static void root_edge_cache_clear(root_edge_cache_entry_t *entries,
                                  size_t count)
{
    size_t i;

    for (i = 0; i < count; i++) {
        key_array_clear(entries[i].edges, entries[i].count);
        free(entries[i].edges);
    }
    free(entries);
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
        size_t first_new;
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

        extract_roots_parallel_word(
            evaluated, count, (size_t)phi->ell + 1,
            root_keys, root_counts, ctx, n_threads);

        first_new = graph->length;
        for (i = 0; i < count; i++) {
            size_t root_index;

            for (root_index = 0; root_index < root_counts[i]; root_index++) {
                jkey_t root_key =
                    root_keys[i * ((size_t)phi->ell + 1) + root_index];
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
        size_t first_new;
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

        extract_roots_parallel_big(
            evaluated, count, (size_t)phi->ell + 1,
            root_keys, root_counts, ctx, n_threads);

        first_new = graph->length;
        for (i = 0; i < count; i++) {
            size_t root_index;

            for (root_index = 0; root_index < root_counts[i]; root_index++) {
                jkey_t root_key =
                    root_keys[i * ((size_t)phi->ell + 1) + root_index];
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

typedef struct {
    size_t *items;
    size_t count;
    size_t task_count;
} index_fill_tasks_t;

typedef struct {
    const size_t *input;
    size_t *temporary;
    size_t *output;
    size_t count;
    size_t task_count;
    size_t *kept_counts;
    size_t *output_offsets;
    const search_graph_t *graph;
    const fmpz *radius;
    unsigned ell;
} index_filter_tasks_t;

static void balanced_partition(size_t count, size_t task_count,
                               size_t task, size_t *start,
                               size_t *length)
{
    size_t base = count / task_count;
    size_t remainder = count % task_count;

    *length = base + (task < remainder ? 1 : 0);
    *start = task * base + (task < remainder ? task : remainder);
}

static void index_fill_task(slong task_index, void *argument)
{
    index_fill_tasks_t *tasks = argument;
    size_t start, length, i;

    balanced_partition(
        tasks->count, tasks->task_count, (size_t)task_index,
        &start, &length);
    for (i = start; i < start + length; i++)
        tasks->items[i] = i;
}

static void index_filter_task(slong task_index, void *argument)
{
    index_filter_tasks_t *tasks = argument;
    size_t task = (size_t)task_index;
    size_t start, length, source, kept = 0;
    fmpz_t next_degree;

    balanced_partition(
        tasks->count, tasks->task_count, task, &start, &length);
    fmpz_init(next_degree);
    for (source = start; source < start + length; source++) {
        size_t index = tasks->input[source];

        fmpz_mul_ui(
            next_degree, tasks->graph->nodes[index].degree,
            tasks->ell);
        if (fmpz_cmp(next_degree, tasks->radius) <= 0)
            tasks->temporary[start + kept++] = index;
    }
    fmpz_clear(next_degree);
    tasks->kept_counts[task] = kept;
}

static void index_filter_copy_task(slong task_index, void *argument)
{
    index_filter_tasks_t *tasks = argument;
    size_t task = (size_t)task_index;
    size_t start, length;

    balanced_partition(
        tasks->count, tasks->task_count, task, &start, &length);
    (void)length;
    memcpy(
        tasks->output + tasks->output_offsets[task],
        tasks->temporary + start,
        tasks->kept_counts[task] * sizeof(*tasks->output));
}

static void index_vector_fill_graph(index_vector_t *vector,
                                    const search_graph_t *graph,
                                    unsigned n_threads)
{
    size_t task_count;
    index_fill_tasks_t tasks;

    vector->length = 0;
    if (vector->alloc < graph->length) {
        vector->alloc = graph->length;
        vector->items = xrealloc(vector->items, vector->alloc,
                                 sizeof(*vector->items));
    }
    vector->length = graph->length;
    if (graph->length < 4096 || n_threads <= 1) {
        size_t i;

        for (i = 0; i < graph->length; i++)
            vector->items[i] = i;
        return;
    }
    task_count = (size_t)n_threads;
    if (task_count > graph->length)
        task_count = graph->length;
    tasks.items = vector->items;
    tasks.count = graph->length;
    tasks.task_count = task_count;
    flint_parallel_do(
        index_fill_task, &tasks, (slong)task_count,
        (int)n_threads, FLINT_PARALLEL_UNIFORM);
}

static void index_vector_filter_eligible(index_vector_t *vector,
                                         const search_graph_t *graph,
                                         unsigned ell,
                                         const fmpz_t radius,
                                         unsigned n_threads)
{
    fmpz_t next_degree;
    size_t source, target = 0;

    if (vector->length >= 4096 && n_threads > 1) {
        index_filter_tasks_t tasks;
        size_t task, task_count = (size_t)n_threads;
        size_t *temporary =
            xrealloc(NULL, vector->length, sizeof(*temporary));
        size_t *metadata =
            xrealloc(NULL, 2 * task_count, sizeof(*metadata));

        if (task_count > vector->length)
            task_count = vector->length;
        tasks.input = vector->items;
        tasks.temporary = temporary;
        tasks.output = vector->items;
        tasks.count = vector->length;
        tasks.task_count = task_count;
        tasks.kept_counts = metadata;
        tasks.output_offsets = metadata + task_count;
        tasks.graph = graph;
        tasks.radius = radius;
        tasks.ell = ell;
        flint_parallel_do(
            index_filter_task, &tasks, (slong)task_count,
            (int)n_threads, FLINT_PARALLEL_UNIFORM);
        target = 0;
        for (task = 0; task < task_count; task++) {
            tasks.output_offsets[task] = target;
            target += tasks.kept_counts[task];
        }
        flint_parallel_do(
            index_filter_copy_task, &tasks, (slong)task_count,
            (int)n_threads, FLINT_PARALLEL_UNIFORM);
        vector->length = target;
        free(metadata);
        free(temporary);
        return;
    }

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
    jkey_t *root_keys;
    size_t *root_counts;
    nmod_poly_t representative;
    fmpz_t next_degree;
    size_t i;
    size_t coefficient_count = phi->n_coefficients;
    size_t root_stride = (size_t)phi->ell + 1;
    int use_fast, found = 0;
    word_evaluation_tasks_t evaluation_tasks;
    word_specialization_tasks_t specialization_tasks;

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
    use_fast = ordered_force_multipoint ||
               (count >= MULTIPOINT_THRESHOLD &&
                phi->max_x_degree >= MULTIPOINT_MIN_DEGREE);
    evaluation_tasks.phi = phi;
    evaluation_tasks.points = points;
    evaluation_tasks.count = count;
    evaluation_tasks.coefficient_values = coefficient_values;
    evaluation_tasks.tree = NULL;
    evaluation_tasks.use_fast = use_fast;
    evaluation_tasks.ctx = ctx;
    evaluate_batch_parallel_word(
        &evaluation_tasks, coefficient_count, n_threads);
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
    clear_batch_parallel_word(
        &evaluation_tasks, coefficient_count, n_threads);

    {
        size_t *dual_indices =
            xrealloc(NULL, count, sizeof(*dual_indices));
        jkey_t *dual_roots =
            xrealloc(NULL, count, sizeof(*dual_roots));
        size_t dual_count = 0;

        for (i = 0; i < count; i++) {
            const search_node_t *current = graph->nodes + frontier[i];

            if (current->parent == NO_INDEX ||
                current->parent_ell != phi->ell)
                continue;
            dual_indices[dual_count] = i;
            dual_roots[dual_count] =
                graph->nodes[current->parent].key;
            dual_count++;
        }
        divide_dual_roots_parallel_word(
            evaluated, dual_indices, dual_roots, dual_count,
            ctx, n_threads);
        ball->stats.skipped_backtracks += dual_count;
        free(dual_roots);
        free(dual_indices);
    }

    ball->stats.specializations += (uint64_t)count;
    ball->stats.root_finds += (uint64_t)count;
    ball->expanded_states += count;

    extract_roots_parallel_word(
        evaluated, count, root_stride, root_keys, root_counts, ctx,
        n_threads);

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
            if (other != NULL) {
                match = graph_find(other, graph->nodes[index].key);
                if (match != NO_INDEX) {
                    *this_meet = index;
                    *other_meet = match;
                    found = 1;
                    break;
                }
            }
        }
    }
    fmpz_clear(next_degree);

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
    jkey_t *root_keys;
    size_t *root_counts;
    fmpz_poly_t representative;
    fmpz_t next_degree;
    size_t i;
    size_t coefficient_count = phi->n_coefficients;
    size_t root_stride = (size_t)phi->ell + 1;
    int use_fast, found = 0;
    big_evaluation_tasks_t evaluation_tasks;
    big_specialization_tasks_t specialization_tasks;

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
    use_fast = ordered_force_multipoint ||
               (count >= MULTIPOINT_THRESHOLD &&
                phi->max_x_degree >= MULTIPOINT_MIN_DEGREE);
    evaluation_tasks.phi = phi;
    evaluation_tasks.points = points;
    evaluation_tasks.count = count;
    evaluation_tasks.coefficient_values = coefficient_values;
    evaluation_tasks.tree = NULL;
    evaluation_tasks.use_fast = use_fast;
    evaluation_tasks.ctx = ctx;
    evaluate_batch_parallel_big(
        &evaluation_tasks, coefficient_count, n_threads);
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
    clear_batch_parallel_big(
        &evaluation_tasks, coefficient_count, n_threads);

    {
        size_t *dual_indices =
            xrealloc(NULL, count, sizeof(*dual_indices));
        jkey_t *dual_roots =
            xrealloc(NULL, count, sizeof(*dual_roots));
        size_t dual_count = 0;

        for (i = 0; i < count; i++) {
            const search_node_t *current = graph->nodes + frontier[i];

            if (current->parent == NO_INDEX ||
                current->parent_ell != phi->ell)
                continue;
            dual_indices[dual_count] = i;
            dual_roots[dual_count] =
                graph->nodes[current->parent].key;
            dual_count++;
        }
        divide_dual_roots_parallel_big(
            evaluated, dual_indices, dual_roots, dual_count,
            ctx, n_threads);
        ball->stats.skipped_backtracks += dual_count;
        free(dual_roots);
        free(dual_indices);
    }

    ball->stats.specializations += (uint64_t)count;
    ball->stats.root_finds += (uint64_t)count;
    ball->expanded_states += count;

    extract_roots_parallel_big(
        evaluated, count, root_stride, root_keys, root_counts, ctx,
        n_threads);

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
            if (other != NULL) {
                match = graph_find(other, graph->nodes[index].key);
                if (match != NO_INDEX) {
                    *this_meet = index;
                    *other_meet = match;
                    found = 1;
                    break;
                }
            }
        }
    }
    fmpz_clear(next_degree);

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

/*
 * Evaluate one Phi_ell batch whose sources may belong to several independent
 * rerandomization balls.  The finite-field work is shared, but every result is
 * relaxed into the graph that owns its source.  This is the key to using one
 * total FLINT thread budget without nesting FLINT worker pools.
 */
static int process_ordered_shared_batch_word(
    const ordered_shared_source_t *sources, size_t count,
    const fmpz_t radius, const modular_poly_t *phi,
    const fq_nmod_ctx_t ctx, unsigned n_threads,
    size_t batch_marker, size_t *winner_slot)
{
    fq_nmod_struct *points, *coefficient_values;
    fq_nmod_poly_struct *evaluated;
    jkey_t *root_keys;
    size_t *root_counts;
    nmod_poly_t representative;
    fmpz_t next_degree;
    size_t i;
    size_t coefficient_count = phi->n_coefficients;
    size_t root_stride = (size_t)phi->ell + 1;
    int use_fast, found = 0;
    word_evaluation_tasks_t evaluation_tasks;
    word_specialization_tasks_t specialization_tasks;

    if (count == 0)
        return 0;
    if (count > (size_t)LONG_MAX)
        die("shared multipoint batch is too large for FLINT's slong "
            "interface");
    if (coefficient_count != 0 &&
        count > SIZE_MAX / coefficient_count)
        die("shared multipoint coefficient buffer size overflow");
    if (root_stride != 0 && count > SIZE_MAX / root_stride)
        die("shared multipoint root buffer size overflow");

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
        search_graph_t *graph = sources[i].lane->ball.graph;

        fq_nmod_init(points + i, ctx);
        fq_nmod_poly_init(evaluated + i, ctx);
        key_to_fq_nmod(
            points + i, graph->nodes[sources[i].node_index].key,
            representative, ctx);
    }
    use_fast = ordered_force_multipoint ||
               (count >= MULTIPOINT_THRESHOLD &&
                phi->max_x_degree >= MULTIPOINT_MIN_DEGREE);
    evaluation_tasks.phi = phi;
    evaluation_tasks.points = points;
    evaluation_tasks.count = count;
    evaluation_tasks.coefficient_values = coefficient_values;
    evaluation_tasks.tree = NULL;
    evaluation_tasks.use_fast = use_fast;
    evaluation_tasks.ctx = ctx;
    evaluate_batch_parallel_word(
        &evaluation_tasks, coefficient_count, n_threads);

    /*
     * A shared batch is one fast/Horner batch for each participating ball,
     * even though only one product tree was built globally.
     */
    for (i = 0; i < count; i++) {
        ordered_shared_lane_t *lane = sources[i].lane;

        if (lane->batch_marker == batch_marker)
            continue;
        lane->batch_marker = batch_marker;
        if (use_fast)
            lane->ball.stats.fast_batches++;
        else
            lane->ball.stats.horner_batches++;
    }

    specialization_tasks.phi = phi;
    specialization_tasks.count = count;
    specialization_tasks.coefficient_values = coefficient_values;
    specialization_tasks.evaluated = evaluated;
    specialization_tasks.ctx = ctx;
    flint_parallel_do(build_specialization_task_word,
                      &specialization_tasks, (slong)count,
                      (int)n_threads, FLINT_PARALLEL_UNIFORM);
    clear_batch_parallel_word(
        &evaluation_tasks, coefficient_count, n_threads);

    {
        size_t *dual_indices =
            xrealloc(NULL, count, sizeof(*dual_indices));
        jkey_t *dual_roots =
            xrealloc(NULL, count, sizeof(*dual_roots));
        size_t dual_count = 0;

        for (i = 0; i < count; i++) {
            ordered_shared_lane_t *lane = sources[i].lane;
            search_graph_t *graph = lane->ball.graph;
            const search_node_t *current =
                graph->nodes + sources[i].node_index;

            if (current->parent == NO_INDEX ||
                current->parent_ell != phi->ell)
                continue;
            dual_indices[dual_count] = i;
            dual_roots[dual_count] =
                graph->nodes[current->parent].key;
            dual_count++;
            lane->ball.stats.skipped_backtracks++;
        }
        divide_dual_roots_parallel_word(
            evaluated, dual_indices, dual_roots, dual_count,
            ctx, n_threads);
        free(dual_roots);
        free(dual_indices);
    }

    for (i = 0; i < count; i++) {
        ordered_shared_lane_t *lane = sources[i].lane;

        lane->ball.stats.specializations++;
        lane->ball.stats.root_finds++;
        lane->ball.expanded_states++;
    }

    extract_roots_parallel_word(
        evaluated, count, root_stride, root_keys, root_counts, ctx,
        n_threads);

    fmpz_init(next_degree);
    for (i = 0; i < count && !found; i++) {
        ordered_shared_lane_t *lane = sources[i].lane;
        search_graph_t *graph = lane->ball.graph;
        size_t source_index = sources[i].node_index;
        size_t root_index;

        fmpz_mul_ui(next_degree, graph->nodes[source_index].degree,
                    phi->ell);
        if (fmpz_cmp(next_degree, radius) > 0)
            continue;
        for (root_index = 0; root_index < root_counts[i]; root_index++) {
            jkey_t root_key = root_keys[i * root_stride + root_index];
            const search_node_t *current =
                graph->nodes + source_index;
            size_t index, match;

            if (current->parent != NO_INDEX &&
                current->parent_ell == phi->ell &&
                key_equal(root_key,
                          graph->nodes[current->parent].key)) {
                lane->ball.stats.skipped_backtracks++;
                continue;
            }
            if (!ordered_relax_node(graph, root_key, next_degree,
                                    source_index, phi->ell, &index))
                continue;
            index_vector_append(&lane->next, index);
            if (lane->other == NULL)
                continue;
            match = graph_find(lane->other, graph->nodes[index].key);
            if (match == NO_INDEX)
                continue;
            *lane->this_meet = index;
            *lane->other_meet = match;
            *lane->job_found = 1;
            *winner_slot = lane->job_slot;
            found = 1;
            break;
        }
    }
    fmpz_clear(next_degree);

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

static int process_ordered_shared_batch_big(
    const ordered_shared_source_t *sources, size_t count,
    const fmpz_t radius, const modular_poly_t *phi,
    const fq_ctx_t ctx, unsigned n_threads,
    size_t batch_marker, size_t *winner_slot)
{
    fq_struct *points, *coefficient_values;
    fq_poly_struct *evaluated;
    jkey_t *root_keys;
    size_t *root_counts;
    fmpz_poly_t representative;
    fmpz_t next_degree;
    size_t i;
    size_t coefficient_count = phi->n_coefficients;
    size_t root_stride = (size_t)phi->ell + 1;
    int use_fast, found = 0;
    big_evaluation_tasks_t evaluation_tasks;
    big_specialization_tasks_t specialization_tasks;

    if (count == 0)
        return 0;
    if (count > (size_t)LONG_MAX)
        die("shared multipoint batch is too large for FLINT's slong "
            "interface");
    if (coefficient_count != 0 &&
        count > SIZE_MAX / coefficient_count)
        die("shared multipoint coefficient buffer size overflow");
    if (root_stride != 0 && count > SIZE_MAX / root_stride)
        die("shared multipoint root buffer size overflow");

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
        search_graph_t *graph = sources[i].lane->ball.graph;

        fq_init(points + i, ctx);
        fq_poly_init(evaluated + i, ctx);
        key_to_fq(points + i,
                  graph->nodes[sources[i].node_index].key,
                  representative, ctx);
    }
    use_fast = ordered_force_multipoint ||
               (count >= MULTIPOINT_THRESHOLD &&
                phi->max_x_degree >= MULTIPOINT_MIN_DEGREE);
    evaluation_tasks.phi = phi;
    evaluation_tasks.points = points;
    evaluation_tasks.count = count;
    evaluation_tasks.coefficient_values = coefficient_values;
    evaluation_tasks.tree = NULL;
    evaluation_tasks.use_fast = use_fast;
    evaluation_tasks.ctx = ctx;
    evaluate_batch_parallel_big(
        &evaluation_tasks, coefficient_count, n_threads);

    for (i = 0; i < count; i++) {
        ordered_shared_lane_t *lane = sources[i].lane;

        if (lane->batch_marker == batch_marker)
            continue;
        lane->batch_marker = batch_marker;
        if (use_fast)
            lane->ball.stats.fast_batches++;
        else
            lane->ball.stats.horner_batches++;
    }

    specialization_tasks.phi = phi;
    specialization_tasks.count = count;
    specialization_tasks.coefficient_values = coefficient_values;
    specialization_tasks.evaluated = evaluated;
    specialization_tasks.ctx = ctx;
    flint_parallel_do(build_specialization_task_big,
                      &specialization_tasks, (slong)count,
                      (int)n_threads, FLINT_PARALLEL_UNIFORM);
    clear_batch_parallel_big(
        &evaluation_tasks, coefficient_count, n_threads);

    {
        size_t *dual_indices =
            xrealloc(NULL, count, sizeof(*dual_indices));
        jkey_t *dual_roots =
            xrealloc(NULL, count, sizeof(*dual_roots));
        size_t dual_count = 0;

        for (i = 0; i < count; i++) {
            ordered_shared_lane_t *lane = sources[i].lane;
            search_graph_t *graph = lane->ball.graph;
            const search_node_t *current =
                graph->nodes + sources[i].node_index;

            if (current->parent == NO_INDEX ||
                current->parent_ell != phi->ell)
                continue;
            dual_indices[dual_count] = i;
            dual_roots[dual_count] =
                graph->nodes[current->parent].key;
            dual_count++;
            lane->ball.stats.skipped_backtracks++;
        }
        divide_dual_roots_parallel_big(
            evaluated, dual_indices, dual_roots, dual_count,
            ctx, n_threads);
        free(dual_roots);
        free(dual_indices);
    }

    for (i = 0; i < count; i++) {
        ordered_shared_lane_t *lane = sources[i].lane;

        lane->ball.stats.specializations++;
        lane->ball.stats.root_finds++;
        lane->ball.expanded_states++;
    }

    extract_roots_parallel_big(
        evaluated, count, root_stride, root_keys, root_counts, ctx,
        n_threads);

    fmpz_init(next_degree);
    for (i = 0; i < count && !found; i++) {
        ordered_shared_lane_t *lane = sources[i].lane;
        search_graph_t *graph = lane->ball.graph;
        size_t source_index = sources[i].node_index;
        size_t root_index;

        fmpz_mul_ui(next_degree, graph->nodes[source_index].degree,
                    phi->ell);
        if (fmpz_cmp(next_degree, radius) > 0)
            continue;
        for (root_index = 0; root_index < root_counts[i]; root_index++) {
            jkey_t root_key = root_keys[i * root_stride + root_index];
            const search_node_t *current =
                graph->nodes + source_index;
            size_t index, match;

            if (current->parent != NO_INDEX &&
                current->parent_ell == phi->ell &&
                key_equal(root_key,
                          graph->nodes[current->parent].key)) {
                lane->ball.stats.skipped_backtracks++;
                continue;
            }
            if (!ordered_relax_node(graph, root_key, next_degree,
                                    source_index, phi->ell, &index))
                continue;
            index_vector_append(&lane->next, index);
            if (lane->other == NULL)
                continue;
            match = graph_find(lane->other, graph->nodes[index].key);
            if (match == NO_INDEX)
                continue;
            *lane->this_meet = index;
            *lane->other_meet = match;
            *lane->job_found = 1;
            *winner_slot = lane->job_slot;
            found = 1;
            break;
        }
    }
    fmpz_clear(next_degree);

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

static int process_ordered_shared_batch(
    const ordered_shared_source_t *sources, size_t count,
    const fmpz_t radius, const modular_poly_t *phi,
    const field_context_t *field, unsigned n_threads,
    size_t batch_marker, size_t *winner_slot)
{
    if (field->kind == BACKEND_FQ_NMOD)
        return process_ordered_shared_batch_word(
            sources, count, radius, phi, field->context.word,
            n_threads, batch_marker, winner_slot);
    return process_ordered_shared_batch_big(
        sources, count, radius, phi, field->context.big,
        n_threads, batch_marker, winner_slot);
}

static size_t ordered_batch_count(size_t remaining)
{
    if (ordered_multipoint_batch_limit != 0 &&
        remaining > ordered_multipoint_batch_limit)
        return ordered_multipoint_batch_limit;
    return remaining;
}

/*
 * Advance several independent balls through one prime stage in lockstep.
 * Sources are interleaved by lane before applying an optional manual batch
 * cap, preventing a small cap from repeatedly favoring the first job.
 */
static int ordered_expand_prime_stage_shared(
    ordered_shared_lane_t *lanes, size_t n_lanes,
    const fmpz_t radius, const modular_poly_t *phi,
    const field_context_t *field, unsigned n_threads,
    size_t *winner_slot, size_t *batch_marker)
{
    ordered_shared_source_t *sources = NULL;
    size_t lane_index, position;
    int found = 0;

    for (lane_index = 0; lane_index < n_lanes; lane_index++) {
        lanes[lane_index].current.length = 0;
        lanes[lane_index].next.length = 0;
        index_vector_fill_graph(&lanes[lane_index].current,
                                lanes[lane_index].ball.graph,
                                n_threads);
    }

    while (!found) {
        size_t total = 0, maximum_length = 0;
        size_t source_count = 0, offset = 0;

        for (lane_index = 0; lane_index < n_lanes; lane_index++) {
            ordered_shared_lane_t *lane = lanes + lane_index;

            index_vector_filter_eligible(
                &lane->current, lane->ball.graph, phi->ell, radius,
                n_threads);
            if (lane->current.length > SIZE_MAX - total)
                die("shared prime-stage source count overflow");
            total += lane->current.length;
            if (lane->current.length > maximum_length)
                maximum_length = lane->current.length;
            lane->next.length = 0;
        }
        if (total == 0)
            break;

        sources = xrealloc(sources, total, sizeof(*sources));
        for (position = 0; position < maximum_length; position++) {
            for (lane_index = 0; lane_index < n_lanes; lane_index++) {
                ordered_shared_lane_t *lane = lanes + lane_index;

                if (position >= lane->current.length)
                    continue;
                sources[source_count].lane = lane;
                sources[source_count].node_index =
                    lane->current.items[position];
                source_count++;
            }
        }
        if (source_count != total)
            die("internal error flattening shared prime-stage sources");

        while (offset < total && !found) {
            size_t count = ordered_batch_count(total - offset);

            (*batch_marker)++;
            if (*batch_marker == 0) {
                for (lane_index = 0; lane_index < n_lanes;
                     lane_index++)
                    lanes[lane_index].batch_marker = 0;
                *batch_marker = 1;
            }
            found = process_ordered_shared_batch(
                sources + offset, count, radius, phi, field,
                n_threads, *batch_marker, winner_slot);
            offset += count;
        }
        if (found)
            break;

        for (lane_index = 0; lane_index < n_lanes; lane_index++) {
            ordered_shared_lane_t *lane = lanes + lane_index;
            index_vector_t temporary;

            index_vector_sort_unique(&lane->next);
            temporary = lane->current;
            lane->current = lane->next;
            lane->next = temporary;
            lane->next.length = 0;
        }
    }

    free(sources);
    return found;
}

/*
 * Generate all canonical ell-power extensions of the endpoints produced by
 * already processed larger primes. Both balls advance one batch at a time.
 * The k-th wave
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

    index_vector_fill_graph(
        &left_current, left_ball->graph, n_threads);
    index_vector_fill_graph(
        &right_current, right_ball->graph, n_threads);

    while (left_current.length != 0 || right_current.length != 0) {
        size_t left_offset = 0, right_offset = 0;

        index_vector_filter_eligible(
            &left_current, left_ball->graph, phi->ell, radius,
            n_threads);
        index_vector_filter_eligible(
            &right_current, right_ball->graph, phi->ell, radius,
            n_threads);
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

/*
 * Enumerate one complete prime stage without testing intersections. This is
 * used by Frobenius mode, where the second ball is represented implicitly by
 * conjugating the completed first ball.
 */
static void ordered_expand_prime_stage_single(
    ordered_ball_t *ball, const fmpz_t radius,
    const modular_poly_t *phi, const field_context_t *field,
    unsigned n_threads, const atomic_int *cancel_requested,
    int *cancelled_out,
    const root_edge_cache_entry_t *cached_root)
{
    index_vector_t current = {0}, next = {0};
    int first_wave = 1;

    index_vector_fill_graph(&current, ball->graph, n_threads);
    while (current.length != 0) {
        size_t offset = 0;

        index_vector_filter_eligible(
            &current, ball->graph, phi->ell, radius, n_threads);
        if (current.length == 0)
            break;
        next.length = 0;

        if (first_wave && cached_root != NULL &&
            current.length != 0 && current.items[0] == 0) {
            fmpz_t cached_degree;
            size_t edge_index;

            fmpz_init_set_ui(cached_degree, phi->ell);
            for (edge_index = 0; edge_index < cached_root->count;
                 edge_index++) {
                size_t index;

                if (ordered_relax_node(
                        ball->graph, cached_root->edges[edge_index],
                        cached_degree, 0, phi->ell, &index))
                    index_vector_append(&next, index);
            }
            fmpz_clear(cached_degree);
            ball->stats.specializations++;
            ball->stats.root_finds++;
            ball->stats.horner_batches++;
            ball->expanded_states++;
            memmove(
                current.items, current.items + 1,
                (current.length - 1) * sizeof(*current.items));
            current.length--;
        }
        first_wave = 0;

        while (offset < current.length) {
            size_t count;
            double started;
            int unexpected_intersection;

            if (cancel_requested != NULL &&
                atomic_load_explicit(cancel_requested,
                                     memory_order_relaxed)) {
                if (cancelled_out != NULL)
                    *cancelled_out = 1;
                goto cleanup;
            }
            count = ordered_batch_count(current.length - offset);
            started = wall_seconds();
            unexpected_intersection = process_ordered_batch(
                ball, current.items + offset, count, radius, phi,
                field, n_threads, NULL, &next, NULL, NULL,
                cancel_requested, cancelled_out);
            ball->active_seconds += wall_seconds() - started;
            if (unexpected_intersection)
                die("internal error: intersection in one-ball enumeration");
            offset += count;
            if (cancelled_out != NULL && *cancelled_out)
                goto cleanup;
        }

        index_vector_sort_unique(&next);
        {
            index_vector_t temporary = current;

            current = next;
            next = temporary;
            next.length = 0;
        }
    }

cleanup:
    index_vector_clear(&next);
    index_vector_clear(&current);
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
        fprintf(stderr,
                "\nEnumeration order = descending canonical prime blocks");
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

    polynomial_index = n_polynomials;
    while (polynomial_index > 0 && !found) {
        polynomial_index--;
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
    const search_graph_t *graph;
    const fmpz *degree_bound;
    size_t *matches;
    atomic_size_t next_index;
    size_t claim_size;
} frobenius_collision_tasks_t;

static void frobenius_collision_task(slong task_index, void *argument)
{
    frobenius_collision_tasks_t *tasks = argument;
    fmpz_t candidate_degree;
    jkey_t conjugate;
    size_t first;

    (void)task_index;
    fmpz_init(candidate_degree);
    key_init(&conjugate);
    for (;;) {
        size_t last, index;

        first = atomic_fetch_add_explicit(
            &tasks->next_index, tasks->claim_size,
            memory_order_relaxed);
        if (first >= tasks->graph->length)
            break;
        last = first + tasks->claim_size;
        if (last > tasks->graph->length)
            last = tasks->graph->length;
        for (index = first; index < last; index++) {
            size_t match;

            frobenius_conjugate(
                &conjugate, &tasks->graph->nodes[index].key);
            match = graph_find(tasks->graph, conjugate);
            if (match != NO_INDEX) {
                fmpz_mul(
                    candidate_degree,
                    tasks->graph->nodes[index].degree,
                    tasks->graph->nodes[match].degree);
                if (fmpz_cmp(candidate_degree,
                             tasks->degree_bound) > 0)
                    match = NO_INDEX;
            }
            tasks->matches[index] = match;
        }
    }
    key_clear(&conjugate);
    fmpz_clear(candidate_degree);
}

static size_t *find_frobenius_collisions_parallel(
    const search_graph_t *graph, const fmpz_t degree_bound,
    unsigned n_threads)
{
    frobenius_collision_tasks_t tasks;
    size_t *matches;
    size_t task_count;

    matches = xrealloc(NULL, graph->length, sizeof(*matches));
    if (graph->length == 0)
        return matches;
    task_count = graph->length < (size_t)n_threads
                     ? graph->length : (size_t)n_threads;
    tasks.graph = graph;
    tasks.degree_bound = degree_bound;
    tasks.matches = matches;
    atomic_init(&tasks.next_index, 0);
    tasks.claim_size =
        adaptive_root_claim_size(graph->length, task_count);
    flint_parallel_do(
        frobenius_collision_task, &tasks, (slong)task_count,
        (int)n_threads, FLINT_PARALLEL_UNIFORM);
    return matches;
}

/*
 * Enumerate the full ball from start, then intersect it with its implicit
 * Frobenius conjugate. If x^p is also in the stored ball, the path to x and
 * the conjugate of the path to x^p join to give a path from start to start^p.
 */
static int ordered_search_frobenius(
    search_graph_t *graph, jkey_t start, const fmpz_t degree_bound,
    const modular_poly_t *polynomials, size_t n_polynomials,
    const field_context_t *field, unsigned n_threads, const char *label,
    size_t *left_meet, size_t *conjugate_meet, fmpz_t path_degree,
    int report, const atomic_int *cancel_requested, int *cancelled_out)
{
    ordered_ball_t ball;
    root_edge_cache_entry_t *root_cache = NULL;
    fmpz_t radius, one, left_degree, conjugate_degree;
    fmpz_t candidate_degree, best_degree;
    jkey_t conjugate;
    size_t polynomial_index, index;
    int found = 0;

    memset(&ball, 0, sizeof(ball));
    ball.graph = graph;
    ball.label = label;
    if (cancelled_out != NULL)
        *cancelled_out = 0;

    fmpz_init(radius);
    fmpz_init(one);
    fmpz_init(left_degree);
    fmpz_init(conjugate_degree);
    fmpz_init(candidate_degree);
    fmpz_init(best_degree);
    key_init(&conjugate);
    fmpz_sqrt(radius, degree_bound);
    fmpz_one(one);
    graph_init(graph);
    graph_add(graph, start, one, NO_INDEX, 0);
    if (n_threads >= 4 && n_polynomials >= 2 &&
        polynomials[n_polynomials - 1].ell + 1 >=
            PARALLEL_FACTOR_MIN_DEGREE) {
        double precompute_started = wall_seconds();

        root_cache = precompute_root_edges(
            start, polynomials, n_polynomials, field, n_threads);
        if (report)
            fprintf(stderr,
                    "Precomputed %lu root-j specializations across "
                    "all primes in %.3f s\n",
                    (unsigned long)n_polynomials,
                    wall_seconds() - precompute_started);
    }

    if (report) {
        fprintf(stderr, "Search radius floor(sqrt(bound)) = ");
        fmpz_fprint(stderr, radius);
        fprintf(stderr,
                "\nFrobenius mode = one explicit ball; "
                "collision scan after complete enumeration\n"
                "Enumeration order = descending canonical prime blocks");
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

    polynomial_index = n_polynomials;
    while (polynomial_index > 0) {
        polynomial_index--;
        if (report)
            fprintf(stderr, "Prime stage ell = %u\n",
                    polynomials[polynomial_index].ell);
        ordered_expand_prime_stage_single(
            &ball, radius, polynomials + polynomial_index, field,
            n_threads, cancel_requested, cancelled_out,
            root_cache == NULL ? NULL
                               : root_cache + polynomial_index);
        if (cancelled_out != NULL && *cancelled_out)
            break;
    }

    if (cancelled_out == NULL || !*cancelled_out) {
        size_t *matches = find_frobenius_collisions_parallel(
            graph, degree_bound, n_threads);

        /*
         * Scan only now: the complete explicit ball is already available.
         * Keep the collision with the smallest product of stored degrees.
         */
        for (index = 0; index < graph->length; index++) {
            size_t match = matches[index];

            if (match == NO_INDEX)
                continue;
            fmpz_mul(candidate_degree, graph->nodes[index].degree,
                     graph->nodes[match].degree);
            if (!found || fmpz_cmp(candidate_degree, best_degree) < 0) {
                *left_meet = index;
                *conjugate_meet = match;
                fmpz_set(best_degree, candidate_degree);
                found = 1;
            }
        }
        free(matches);
    }

    if (found) {
        graph_chain_degree(left_degree, graph, *left_meet);
        graph_chain_degree(conjugate_degree, graph, *conjugate_meet);
        if (fmpz_cmp(left_degree, radius) > 0 ||
            fmpz_cmp(conjugate_degree, radius) > 0)
            die("internal error: Frobenius collision exceeds search radius");
        fmpz_mul(path_degree, left_degree, conjugate_degree);
        if (fmpz_cmp(path_degree, degree_bound) > 0)
            die("internal error: Frobenius path exceeds degree bound");
    }

    if (report) {
        ordered_ball_report(&ball, 0);
        fprintf(stderr,
                "Implicit conjugate ball: %lu j-invariants "
                "(not stored separately)\n",
                (unsigned long)graph->length);
        fprintf(stderr, "Post-enumeration Frobenius collision scan: %s\n",
                found ? "intersection found" : "no intersection");
    }

    key_clear(&conjugate);
    root_edge_cache_clear(root_cache, root_cache == NULL
                                          ? 0 : n_polynomials);
    fmpz_clear(best_degree);
    fmpz_clear(candidate_degree);
    fmpz_clear(conjugate_degree);
    fmpz_clear(left_degree);
    fmpz_clear(one);
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
    const char *right_label;
    int slot;
    int single_ball_frobenius;
    int found;
    int cancelled;
    double elapsed_seconds;
} rerandomization_search_job_t;

static int find_frobenius_collision(
    search_graph_t *graph, const fmpz_t degree_bound,
    const fmpz_t radius, size_t *left_meet,
    size_t *conjugate_meet, fmpz_t path_degree,
    unsigned n_threads)
{
    fmpz_t candidate_degree, best_degree;
    fmpz_t left_degree, conjugate_degree;
    size_t *matches;
    size_t index;
    int found = 0;

    fmpz_init(candidate_degree);
    fmpz_init(best_degree);
    fmpz_init(left_degree);
    fmpz_init(conjugate_degree);
    matches = find_frobenius_collisions_parallel(
        graph, degree_bound, n_threads);

    for (index = 0; index < graph->length; index++) {
        size_t match = matches[index];

        if (match == NO_INDEX)
            continue;
        fmpz_mul(candidate_degree, graph->nodes[index].degree,
                 graph->nodes[match].degree);
        if (!found || fmpz_cmp(candidate_degree, best_degree) < 0) {
            *left_meet = index;
            *conjugate_meet = match;
            fmpz_set(best_degree, candidate_degree);
            found = 1;
        }
    }

    if (found) {
        graph_chain_degree(left_degree, graph, *left_meet);
        graph_chain_degree(conjugate_degree, graph, *conjugate_meet);
        if (fmpz_cmp(left_degree, radius) > 0 ||
            fmpz_cmp(conjugate_degree, radius) > 0)
            die("internal error: shared Frobenius collision exceeds "
                "search radius");
        fmpz_mul(path_degree, left_degree, conjugate_degree);
        if (fmpz_cmp(path_degree, degree_bound) > 0)
            die("internal error: shared Frobenius path exceeds degree "
                "bound");
    }

    free(matches);
    fmpz_clear(conjugate_degree);
    fmpz_clear(left_degree);
    fmpz_clear(best_degree);
    fmpz_clear(candidate_degree);
    return found;
}

/*
 * Run a complete batch of rerandomization candidates under one global FLINT
 * worker budget.  No controller pthread owns a private search.  Instead, all
 * active balls contribute points to the same prime-stage batches.
 */
static int run_rerandomization_search_batch_shared(
    rerandomization_search_job_t *jobs, size_t n_jobs,
    unsigned n_threads)
{
    ordered_shared_lane_t *lanes;
    size_t n_lanes, lane_index = 0, job_index;
    size_t polynomial_index, winner_slot = NO_INDEX;
    size_t batch_marker = 0;
    fmpz_t radius, one, left_degree, right_degree;
    double started, elapsed;
    int single_ball_frobenius, stopped_early = 0;

    if (n_jobs == 0)
        return -1;
    single_ball_frobenius = jobs[0].single_ball_frobenius;
    for (job_index = 1; job_index < n_jobs; job_index++)
        if (jobs[job_index].single_ball_frobenius !=
            single_ball_frobenius)
            die("internal error: mixed rerandomization search modes");
    if (!single_ball_frobenius && n_jobs > SIZE_MAX / 2)
        die("too many shared rerandomization lanes");
    n_lanes = n_jobs * (single_ball_frobenius ? 1 : 2);
    lanes = xrealloc(NULL, n_lanes, sizeof(*lanes));
    memset(lanes, 0, n_lanes * sizeof(*lanes));

    fmpz_init(radius);
    fmpz_init(one);
    fmpz_init(left_degree);
    fmpz_init(right_degree);
    fmpz_sqrt(radius, jobs[0].degree_bound);
    fmpz_one(one);
    started = wall_seconds();

    for (job_index = 0; job_index < n_jobs; job_index++) {
        rerandomization_search_job_t *job = jobs + job_index;
        ordered_shared_lane_t *left_lane = lanes + lane_index++;

        job->found = 0;
        job->cancelled = 0;
        graph_init(&job->left);
        graph_add(&job->left, job->start, one, NO_INDEX, 0);
        left_lane->ball.graph = &job->left;
        left_lane->ball.label = "j'";
        left_lane->this_meet = &job->left_meet;
        left_lane->other_meet = &job->right_meet;
        left_lane->job_found = &job->found;
        left_lane->job_slot = (size_t)job->slot;

        if (single_ball_frobenius)
            continue;

        graph_init(&job->right);
        graph_add(&job->right, job->target, one, NO_INDEX, 0);
        left_lane->other = &job->right;
        {
            ordered_shared_lane_t *right_lane = lanes + lane_index++;

            right_lane->ball.graph = &job->right;
            right_lane->ball.label = job->right_label;
            right_lane->other = &job->left;
            right_lane->this_meet = &job->right_meet;
            right_lane->other_meet = &job->left_meet;
            right_lane->job_found = &job->found;
            right_lane->job_slot = (size_t)job->slot;
        }
        if (key_equal(job->start, job->target) &&
            winner_slot == NO_INDEX) {
            job->found = 1;
            job->left_meet = 0;
            job->right_meet = 0;
            fmpz_one(job->path_degree);
            winner_slot = (size_t)job->slot;
            stopped_early = 1;
        }
    }
    if (lane_index != n_lanes)
        die("internal error constructing shared rerandomization lanes");

    polynomial_index = jobs[0].n_polynomials;
    while (polynomial_index > 0 && winner_slot == NO_INDEX) {
        polynomial_index--;
        if (ordered_expand_prime_stage_shared(
                lanes, n_lanes, radius,
                jobs[0].polynomials + polynomial_index,
                jobs[0].field, n_threads, &winner_slot,
                &batch_marker)) {
            stopped_early = 1;
            break;
        }
    }

    if (single_ball_frobenius) {
        /*
         * Preserve the one-ball Frobenius rule: every explicit ball is fully
         * enumerated before any conjugate collision is tested.
         */
        for (job_index = 0; job_index < n_jobs; job_index++) {
            rerandomization_search_job_t *job = jobs + job_index;

            job->found = find_frobenius_collision(
                &job->left, job->degree_bound, radius,
                &job->left_meet, &job->right_meet,
                job->path_degree, n_threads);
            if (job->found && winner_slot == NO_INDEX)
                winner_slot = (size_t)job->slot;
        }
    } else if (winner_slot != NO_INDEX) {
        rerandomization_search_job_t *job = jobs + winner_slot;

        graph_chain_degree(left_degree, &job->left, job->left_meet);
        graph_chain_degree(right_degree, &job->right, job->right_meet);
        if (fmpz_cmp(left_degree, radius) > 0 ||
            fmpz_cmp(right_degree, radius) > 0)
            die("internal error: shared meeting endpoint exceeds "
                "search radius");
        fmpz_mul(job->path_degree, left_degree, right_degree);
        if (fmpz_cmp(job->path_degree, job->degree_bound) > 0)
            die("internal error: shared meeting path exceeds degree "
                "bound");
    }

    elapsed = wall_seconds() - started;
    for (job_index = 0; job_index < n_jobs; job_index++) {
        rerandomization_search_job_t *job = jobs + job_index;

        job->elapsed_seconds = elapsed;
        if (stopped_early && (size_t)job->slot != winner_slot &&
            !job->found)
            job->cancelled = 1;
    }

    for (lane_index = 0; lane_index < n_lanes; lane_index++) {
        index_vector_clear(&lanes[lane_index].next);
        index_vector_clear(&lanes[lane_index].current);
    }
    free(lanes);
    fmpz_clear(right_degree);
    fmpz_clear(left_degree);
    fmpz_clear(one);
    fmpz_clear(radius);
    return winner_slot == NO_INDEX ? -1 : (int)winner_slot;
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

static void print_frobenius_path(
    const search_graph_t *graph, size_t left_meet,
    size_t conjugate_meet, const fmpz_t total_degree)
{
    size_t *reverse_path = NULL;
    size_t length = 0, conjugate_length = 0, alloc = 0, i, current;
    jkey_t conjugate;

    current = left_meet;
    while (current != NO_INDEX) {
        if (length == alloc) {
            alloc = alloc ? 2 * alloc : 16;
            reverse_path = xrealloc(reverse_path, alloc,
                                    sizeof(*reverse_path));
        }
        reverse_path[length++] = current;
        current = graph->nodes[current].parent;
    }

    current = conjugate_meet;
    while (graph->nodes[current].parent != NO_INDEX) {
        conjugate_length++;
        current = graph->nodes[current].parent;
    }

    printf("Path found. Total degree = ");
    fmpz_print(total_degree);
    printf("\nNumber of steps = %lu\n\n",
           (unsigned long)(length - 1 + conjugate_length));

    print_key(graph->nodes[reverse_path[length - 1]].key);
    printf("\n");
    for (i = length - 1; i > 0; i--) {
        const search_node_t *child = graph->nodes + reverse_path[i - 1];

        printf("  --[%u]--> ", child->parent_ell);
        print_key(child->key);
        printf("\n");
    }

    key_init(&conjugate);
    current = conjugate_meet;
    while (graph->nodes[current].parent != NO_INDEX) {
        printf("  --[%u]--> ", graph->nodes[current].parent_ell);
        current = graph->nodes[current].parent;
        frobenius_conjugate(&conjugate, &graph->nodes[current].key);
        print_key(conjugate);
        printf("\n");
    }
    key_clear(&conjugate);
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
                                    int conjugate_middle_return,
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
    if (conjugate_middle_return) {
        jkey_t conjugate;

        key_init(&conjugate);
        while (right->nodes[current].parent != NO_INDEX) {
            printf("  --[%u]--> ", right->nodes[current].parent_ell);
            current = right->nodes[current].parent;
            frobenius_conjugate(
                &conjugate, &right->nodes[current].key);
            print_key(conjugate);
            printf("\n");
        }
        key_clear(&conjugate);
    } else {
        while (right->nodes[current].parent != NO_INDEX) {
            printf("  --[%u]--> ", right->nodes[current].parent_ell);
            current = right->nodes[current].parent;
            print_key(right->nodes[current].key);
            printf("\n");
        }
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
    if (errno || *end != '\0' || bits < 2 ||
        bits > MAX_CHARACTERISTIC_BITS)
        die("pbits must be a decimal integer between 2 and 127");
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
        big_component_from_fmpz(&key->value.big.re, re);
        big_component_from_fmpz(&key->value.big.im, im);
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
    quadratic_basis_is_minus_one =
        fmpz_equal(quadratic_nonresidue, p) == 0 &&
        fmpz_cmp_ui(p, 2) > 0 &&
        fmpz_cmp_ui(quadratic_nonresidue, 0) > 0;
    if (quadratic_basis_is_minus_one) {
        fmpz_t minus_one;

        fmpz_init(minus_one);
        fmpz_sub_ui(minus_one, p, 1);
        quadratic_basis_is_minus_one =
            fmpz_equal(quadratic_nonresidue, minus_one);
        fmpz_clear(minus_one);
    }

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

static cubic_kernel_kind_t parse_cubic_kernel(const char *s)
{
    if (strcmp(s, "flint") == 0)
        return CUBIC_KERNEL_FLINT;
    if (strcmp(s, "batched") == 0)
        return CUBIC_KERNEL_BATCHED;
    if (strcmp(s, "compare") == 0)
        return CUBIC_KERNEL_COMPARE;
    die("--cubic-kernel must be 'flint', 'batched', or 'compare'");
    return CUBIC_KERNEL_FLINT;
}

static const char *cubic_kernel_name(cubic_kernel_kind_t kernel)
{
    switch (kernel) {
        case CUBIC_KERNEL_FLINT:
            return "flint";
        case CUBIC_KERNEL_BATCHED:
            return "batched";
        case CUBIC_KERNEL_COMPARE:
            return "compare";
    }
    return "unknown";
}

static void print_root_profile_report(void)
{
    uint64_t total_count = 0, total_cpu_ns = 0;
    uint64_t control_wall =
        atomic_load_explicit(&cubic_control_wall_ns,
                             memory_order_relaxed);
    uint64_t batched_wall =
        atomic_load_explicit(&cubic_batched_wall_ns,
                             memory_order_relaxed);
    uint64_t batched_count =
        atomic_load_explicit(&cubic_batched_counts,
                             memory_order_relaxed);
    uint64_t batched_cpu =
        atomic_load_explicit(&cubic_batched_cpu_ns,
                             memory_order_relaxed);
    uint64_t fallback_count =
        atomic_load_explicit(&cubic_batched_fallbacks,
                             memory_order_relaxed);
    unsigned degree;

    if (root_profile_enabled) {
        for (degree = 0; degree < ROOT_PROFILE_DEGREES; degree++) {
            total_count += atomic_load_explicit(
                root_profile_counts + degree, memory_order_relaxed);
            total_cpu_ns += atomic_load_explicit(
                root_profile_cpu_ns + degree, memory_order_relaxed);
        }

        fprintf(stderr,
                "\nResidual-degree root profile "
                "(summed worker thread CPU)\n");
        fprintf(stderr,
                "  degree       calls       CPU s    us/call   CPU share\n");
        for (degree = 0; degree < ROOT_PROFILE_DEGREES; degree++) {
            uint64_t count = atomic_load_explicit(
                root_profile_counts + degree, memory_order_relaxed);
            uint64_t cpu_ns = atomic_load_explicit(
                root_profile_cpu_ns + degree, memory_order_relaxed);

            if (count == 0)
                continue;
            fprintf(stderr,
                    "  %6u  %10llu  %10.6f  %9.3f  %9.2f%%\n",
                    degree, (unsigned long long)count,
                    (double)cpu_ns * 1.0e-9,
                    (double)cpu_ns * 1.0e-3 / (double)count,
                    total_cpu_ns == 0
                        ? 0.0
                        : 100.0 * (double)cpu_ns /
                              (double)total_cpu_ns);
        }
        fprintf(stderr,
                "  total   %10llu  %10.6f\n",
                (unsigned long long)total_count,
                (double)total_cpu_ns * 1.0e-9);
    }

    if (cubic_kernel != CUBIC_KERNEL_FLINT) {
        fprintf(stderr,
                "\nBatched cubic prototype "
                "(%u lanes, exact modular-polynomial specializations)\n",
                CUBIC_BATCH_LANES);
        fprintf(stderr,
                "  batched cubics:         %llu\n"
                "  batched worker CPU:     %.6f s\n"
                "  generic fallbacks:      %llu\n",
                (unsigned long long)batched_count,
                (double)batched_cpu * 1.0e-9,
                (unsigned long long)fallback_count);
        if (cubic_kernel == CUBIC_KERNEL_COMPARE) {
            uint64_t control_cpu = atomic_load_explicit(
                root_profile_cpu_ns + 3, memory_order_relaxed);

            fprintf(stderr,
                    "  FLINT cubic pass wall:  %.6f s\n"
                    "  batched pass wall:      %.6f s\n",
                    (double)control_wall * 1.0e-9,
                    (double)batched_wall * 1.0e-9);
            if (batched_cpu != 0)
                fprintf(stderr,
                        "  CPU speed ratio "
                        "(FLINT/batched): %.3fx\n",
                        (double)control_cpu / (double)batched_cpu);
            if (batched_wall != 0)
                fprintf(stderr,
                        "  wall speed ratio "
                        "(FLINT/batched): %.3fx\n",
                        (double)control_wall / (double)batched_wall);
            fprintf(stderr,
                    "  verification: every sorted cubic root set "
                    "matched FLINT\n");
        }
    }
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
    } else if (result->value.big.im.lo != 0 ||
               result->value.big.im.hi != 0) {
        fmpz_t imaginary;

        fmpz_init(imaginary);
        big_component_to_fmpz(imaginary, result->value.big.im);
        fmpz_sub(imaginary, field_characteristic, imaginary);
        big_component_from_fmpz(&result->value.big.im, imaginary);
        fmpz_clear(imaginary);
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
        key_set_ui(&current, UWORD(1728), 0);
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
    int found, left_graph_live = 0, right_graph_live = 0;
    int attempted_rerandomization = 0;
    int rerandomization_phi_loaded = 0;
    int automatic_j = 0, seed_supplied = 0, threads_supplied = 0;
    int multipoint_batch_supplied = 0, force_multipoint_supplied = 0;
    int evaluation_shards_supplied = 0;
    int parallel_factor_tree_supplied = 0;
    int root_profile_supplied = 0, cubic_kernel_supplied = 0;
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
                          "--evaluation-shards") == 0) {
            if (evaluation_shards_supplied ||
                argument_index + 1 == argc)
                die("--evaluation-shards must occur once and have a "
                    "value");
            evaluation_shard_limit =
                parse_multipoint_batch_limit(
                    argv[++argument_index]);
            evaluation_shards_supplied = 1;
        } else if (strcmp(argv[argument_index],
                          "--root-profile") == 0) {
            if (root_profile_supplied)
                die("--root-profile must occur at most once");
            root_profile_enabled = 1;
            root_profile_supplied = 1;
        } else if (strcmp(argv[argument_index],
                          "--parallel-factor-tree") == 0) {
            if (parallel_factor_tree_supplied ||
                argument_index + 1 == argc)
                die("--parallel-factor-tree must occur once and have "
                    "a value");
            parallel_factor_tree_enabled =
                parse_on_off(argv[++argument_index]);
            parallel_factor_tree_supplied = 1;
        } else if (strcmp(argv[argument_index],
                          "--cubic-kernel") == 0) {
            if (cubic_kernel_supplied ||
                argument_index + 1 == argc)
                die("--cubic-kernel must occur once and have a value");
            cubic_kernel =
                parse_cubic_kernel(argv[++argument_index]);
            cubic_kernel_supplied = 1;
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
                "[--evaluation-shards n] [--force-multipoint]\n"
                "          [--root-profile] "
                "[--parallel-factor-tree on|off]\n"
                "[--cubic-kernel flint|batched|compare]\n"
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
    if (fmpz_bits(field_characteristic) > MAX_CHARACTERISTIC_BITS)
        die("this optimized build requires p < 2^127");
    if (cubic_kernel == CUBIC_KERNEL_COMPARE)
        root_profile_enabled = 1;
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
                : "fq with fixed two-limb graph keys "
                  "(65--127-bit characteristic)");
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
    fprintf(stderr, "Total FLINT worker budget = %u\n", n_threads);
    if (ordered_multipoint_batch_limit == 0)
        fprintf(stderr,
                "Multipoint batch limit = unlimited "
                "(complete eligible prime-power layer)\n");
    else
        fprintf(stderr, "Multipoint batch limit = %lu%s\n",
                (unsigned long)ordered_multipoint_batch_limit,
                multipoint_batch_supplied ? " [explicit]" : "");
    if (evaluation_shard_limit == 0)
        fprintf(stderr,
                "Evaluation shards = automatic (up to worker budget)\n");
    else
        fprintf(stderr, "Evaluation shards = at most %lu%s\n",
                (unsigned long)evaluation_shard_limit,
                evaluation_shards_supplied ? " [explicit]" : "");
    fprintf(stderr, "Force multipoint evaluation = %s%s\n",
            ordered_force_multipoint ? "on" : "off",
            force_multipoint_supplied ? " [explicit]" : " [default]");
    fprintf(stderr, "Cubic root kernel = %s%s\n",
            cubic_kernel_name(cubic_kernel),
            cubic_kernel_supplied ? " [explicit]" : " [default]");
    fprintf(stderr, "Residual-degree root profile = %s%s\n",
            root_profile_enabled ? "on" : "off",
            root_profile_supplied ? " [explicit]" :
            cubic_kernel == CUBIC_KERNEL_COMPARE
                ? " [implied by compare]"
                : " [default]");
    fprintf(stderr, "Parallel high-degree factor tree = %s%s\n",
            parallel_factor_tree_enabled ? "on" : "off",
            parallel_factor_tree_supplied ? " [explicit]"
                                          : " [default]");
    if (allow_rerandomization)
        fprintf(stderr,
                "Concurrent rerandomization balls = %u%s "
                "(shared worker scheduler)\n",
                rerandomization_jobs,
                rerandomization_jobs_supplied ? " [explicit]"
                                              : " [default]");

    initial_start = wall_seconds();
    if (j_mode == 2) {
        found = ordered_search_between(
            &left, &right, j1, target, N,
            polynomials, n_polynomials, &field, n_threads,
            "j1", "j2", &left_meet, &right_meet,
            path_degree, 1, NULL, NULL);
        left_graph_live = 1;
        right_graph_live = 1;
    } else {
        found = ordered_search_frobenius(
            &left, j1, N, polynomials, n_polynomials, &field,
            n_threads, "j1", &left_meet, &right_meet,
            path_degree, 1, NULL, NULL);
        left_graph_live = 1;
    }
    initial_end = wall_seconds();

    if (found) {
        if (j_mode == 2)
            print_path(
                &left, left_meet, &right, right_meet, path_degree);
        else
            print_frobenius_path(
                &left, left_meet, right_meet, path_degree);
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

        if (left_graph_live) {
            graph_clear(&left);
            left_graph_live = 0;
        }
        if (right_graph_live) {
            graph_clear(&right);
            right_graph_live = 0;
        }

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
                int winner;

                if (batch_count > (size_t)rerandomization_jobs)
                    batch_count = (size_t)rerandomization_jobs;
                jobs = xrealloc(NULL, batch_count, sizeof(*jobs));

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
                    job->right_label =
                        j_mode == 2 ? "j2" : "(j')^p";
                    job->single_ball_frobenius = j_mode != 2;
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

                winner = run_rerandomization_search_batch_shared(
                    jobs, batch_count, n_threads);

                attempted_rerandomization = 1;
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
                    if (job->single_ball_frobenius)
                        fprintf(stderr,
                                "Rerandomization attempt %lu %s in "
                                "%.3f s (one explicit ball %lu; "
                                "conjugate ball implicit)\n",
                                (unsigned long)job->attempt_number,
                                status, job->elapsed_seconds,
                                (unsigned long)job->left.length);
                    else
                        fprintf(stderr,
                                "Rerandomization attempt %lu %s in "
                                "%.3f s (left ball %lu, "
                                "right ball %lu)\n",
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
                    left_graph_live = 1;
                    if (!winning_job->single_ball_frobenius) {
                        right = winning_job->right;
                        right_graph_live = 1;
                    }
                    left_meet = winning_job->left_meet;
                    right_meet = winning_job->right_meet;
                    fmpz_set(path_degree,
                             winning_job->path_degree);
                }

                for (job_index = 0; job_index < batch_count;
                     job_index++) {
                    rerandomization_search_job_t *job =
                        jobs + job_index;

                    if ((int)job_index != winner) {
                        graph_clear(&job->left);
                        if (!job->single_ball_frobenius)
                            graph_clear(&job->right);
                    }
                    fmpz_clear(job->path_degree);
                    key_clear(&job->target);
                    key_clear(&job->start);
                }
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
            &left, left_meet,
            j_mode == 2 ? &right : &left, right_meet, successful_walk,
            successful_walk_length, rerandomization_ell, path_degree,
            j_mode != 2,
            j_mode != 2);
    }
    fflush(stdout);
    output_end = wall_seconds();
    print_root_profile_report();

    cleanup_start = wall_seconds();
    if (left_graph_live)
        graph_clear(&left);
    if (right_graph_live)
        graph_clear(&right);
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

    printf("\nTiming (wall clock, total FLINT worker budget %u",
           n_threads);
    if (allow_rerandomization)
        printf(", %u concurrent rerandomization ball%s",
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
