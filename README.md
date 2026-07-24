**WARNING: The files in this repository were created with assistance from ChatGPT Sol 5.6 (i.e. *vibecoded*). It is research/prototype software and may contain serious mathematical or implementation mistakes. Use at your own risk.**

## Programs

### `find_smooth_path.c`

Searches for an isogeny path between supersingular j-invariants using
classical modular polynomials.

The middle portion of the path consists of isogenies of prime degree
ell <= B, and its total degree is bounded by N. The program
enumerates the two degree-bounded balls simultaneously using interleaved
Dijkstra searches and stops as soon as they intersect.

The implementation:

- loads each required modular polynomial once and reduces it modulo p;
- batches equal-distance vertices and uses FLINT multipoint evaluation when
  the batch is sufficiently large;
- evaluates and factors independent polynomials in parallel;
- uses `fq_nmod` when p fits in a machine word;
- automatically switches to the arbitrary-precision `fq` backend for larger
  primes;
- reports wall-clock timings and peak resident memory.

Usage:

```text
./find_smooth_path [N [phi_directory]] [options]
```

Options:

```text
--p p                    field characteristic
--pbits n                largest prime p < 2^n with p = 3 (mod 4)
--B B                    smoothness bound, 2 <= B <= 100
--N N                    degree bound for the middle smooth path
--phi-dir directory      modular-polynomial directory
--j1 re im               starting j-invariant
--j2 re im               target j-invariant; requires --j1
--rerandomization on|off override the mode-dependent default
--threads n              number of FLINT worker threads, 1 <= n <= 256
--seed s                 64-bit seed for automatic j-generation
```

Defaults:

- p = 4294967311;
- N = floor(cuberoot(floor(p/2)));
- `phi_directory=mod_pols`;
- B is selected by an empirical prime-valued step function of
  log(N).

The options `--p p` and `--pbits n` are mutually exclusive. If neither is
given, the finder retains its default p = 4294967311.

There are three input modes:

1. With no supplied j-invariant, the program generates one by taking a
   random non-backtracking 2-isogeny walk from j = 1728, then searches for
   a path from j to j^p.
2. With `--j1 re im`, it searches from the supplied j to j^p.
3. With both `--j1` and `--j2`, it searches from j1 to j2.

Rerandomization is enabled by default in the first two modes and disabled by
default in the third. When enabled after a failed search, the program
enumerates non-backtracking walks of degree ell_prime, where ell_prime is the
smallest prime larger than B, until a rerandomized endpoint
admits a smooth middle path. Use `--rerandomization off` when a bounded
unsuccessful search should terminate immediately.

Examples:

```bash
# Use all principal defaults, with reproducible automatic j-generation.
./find_smooth_path --threads 8 --seed 12345

# Generate the characteristic from its requested bit length.
./find_smooth_path \
    --pbits 100 \
    --threads 8 \
    --seed 12345

# Supply j1 and search for a path to its Frobenius conjugate.
./find_smooth_path \
    --p 2305843009213693951 \
    --j1 1992837699471099977 521720088548403283 \
    --B 41 \
    --N 1048576 \
    --phi-dir mod_pols \
    --threads 8

# Search between two explicitly supplied j-invariants.
./find_smooth_path \
    --p 2305843009213693951 \
    --j1 1992837699471099977 521720088548403283 \
    --j2 1992837699471099977 1784122920665290668 \
    --B 41 \
    --N 1048576 \
    --phi-dir mod_pols \
    --rerandomization off \
    --threads 8
```

The finder exits with status 0 when it finds a path, status 2 when a bounded
search terminates without finding one, and a nonzero error status for invalid
input or missing data.

### `random_walk.c`

Computes a random non-backtracking 2-isogeny walk starting at a supplied
supersingular j-invariant.

At every step, the program specializes the classical modular polynomial
Phi_2(X,Y) at the current j-invariant, finds its roots in
Fp2, removes one copy of the edge just traversed, and chooses
one of the remaining edges.

The program uses FLINT's `fq_nmod` backend when p fits in one machine
word and automatically switches to the arbitrary-precision `fq` backend for
larger primes. The coefficients `j_re` and `j_im` represent j = j_re + j_im*a.

Usage:

```text
./random_walk (--p p | --pbits n) \
    [--j j_re j_im] [--steps length] [--seed seed]
```

For compatibility, the following positional form is also accepted:

```text
./random_walk p [j_re j_im [length]] [--seed seed]
```

The default starting j-invariant is 1728, and the default walk length
is 2*ceil(log2(p)). When the default starting invariant is used,
the program requires p = 3 mod 4, so that 1728 is supersingular.

Example:

```bash
./random_walk \
    --p 2305843009213693951 \
    --j 1992837699471099977 521720088548403283 \
    --steps 20 \
    --seed 12345

# Choose p as the largest prime below 2^100 with p = 3 (mod 4),
# start at j = 1728, and use the default 200-step walk.
./random_walk --pbits 100 --seed 12345
```

The optional 64-bit seed makes the walk reproducible. The internal
SplitMix64 generator is suitable for graph experiments but is not a
cryptographic random-number generator.

## Field representation

Both programs use the field Fp2 = Fp[a]/((a^2 - q)), where q is a
deterministic quadratic nonresidue modulo p. If p = 3 mod 4, the programs
choose q = -1, so the field is represented as Fp[a]/((a^2 + 1)).

In this common case, the Frobenius conjugate of `(re, im)` is
`(re, p-im)` when `im` is nonzero.

Automatic j-generation from 1728 requires p = 3 mod 4, which
ensures that j = 1728 is supersingular. Supplied j-invariants are
assumed to be supersingular; the programs do not test this condition.

## Requirements

- a C11 compiler;
- FLINT 3.x;
- GMP;
- MPFR;
- POSIX threads.

## Compilation

```bash
gcc -O3 -march=native -DNDEBUG -std=gnu11 -Wall -Wextra \
    -Wno-unused-parameter \
    find_smooth_path.c \
    -pthread -lflint -lmpfr -lgmp \
    -o find_smooth_path

gcc -O3 -march=native -DNDEBUG -std=gnu11 -Wall -Wextra \
    random_walk.c \
    -lflint -lmpfr -lgmp \
    -o random_walk
```

For a portable binary, omit `-march=native`.

## Modular-polynomial files

`find_smooth_path` requires one classical modular polynomial
Phi_ell(X,Y) for every prime ell <= B. If rerandomization is
enabled, it also requires Phi_ell_prime, where ell_prime is the smallest
prime larger than B.

By default, the program searches in `mod_pols/` for either of these names:

```text
mod_pols/Phi_<ell>.txt
mod_pols/phi_j_<ell>.txt
```

It accepts two sparse text formats:

```text
x_exponent y_exponent integer_coefficient
[x_exponent,y_exponent] integer_coefficient
```

The bracketed format stores one symmetric half of the modular polynomial; the
swapped monomial is inserted automatically.

### Acknowledgement

The classical modular polynomials used by this software were obtained from Andrew V. Sutherland's [Modular Polynomials database](https://math.mit.edu/~drew/ClassicalModPolys.html).

If you use these modular polynomials in your research, please cite:

> Reinier Bröker, Kristin Lauter, and Andrew V. Sutherland, “Modular polynomials via isogeny volcanoes,” *Mathematics of Computation* 81 (2012), 1201–1231. [doi:10.1090/S0025-5718-2011-02508-1](https://doi.org/10.1090/S0025-5718-2011-02508-1).
