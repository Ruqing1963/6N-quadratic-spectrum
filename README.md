# 6N-quadratic-spectrum

**Part XXXII — The Richness Spectrum of Quadratic Primes on the 6N Skeleton: Mask Rotation in n²+a and the Heegner Corridor of n²+n+A**

Ruqing Chen · GUT Geoservice Inc., Montreal · June 2026

Companion code and data for Part XXXII of *Arithmetic Geodynamics on the 6N Skeleton*. Two families of
quadratics, one principle: **a quadratic is prime-rich exactly insofar as its discriminant dodges the
small primes.** **Everything here is a measured sieve result — no fitted parameters, no fabricated
numbers.**

## One principle

For a quadratic `f` and prime `q`, let `ω(q)` = #roots of `f ≡ 0 (mod q)`. The Bateman–Horn richness
is `C = Π_q g(q)`, with local factor `g(q) = (1 − ω(q)/q)/(1 − 1/q)`:
- `ω(q)=0` (q **cannot** divide f) → `g = q/(q−1) > 1`: a **dodge / enhancement**.
- `ω(q)=2` (q divides on two residues) → `g = (q−2)/(q−1) < 1`: a **suppression**.

Richness is the product of dodged vs hit primes. Discriminants that are non-residues mod many small
primes win.

## The two families

**(i) n²+a — the shift rotates a mod-6 mask.** Which break-zones are open is set by `a mod 6`:

| a mod 6 | open wings | example C(a) |
|---|---|---|
| 1, 4 | both, right:left = 1:2 | a=1:1.37, a=7:1.97, a=−2:1.85 |
| 0, 3 | right only (left sealed) | a=3:1.12, a=6:0.71 |
| 2, 5 | left only (right sealed) | a=2:0.71, a=5:0.53 |

The mask sets *which* zones; the richness C(a) rides independently on top (0.53 → ~2). Bateman–Horn
holds across the family to a fraction of a percent. (a=−1 → n²−1=(n−1)(n+1) factors → C→0, correctly.)

**(ii) n²+n+A — the Heegner corridor.** Discriminant `1−4A`. The class-number-one (Heegner) values
`A = 3, 5, 11, 17, 41` (with `4A−1 = 11,19,43,67,163`) thread the longest corridors — the first prime
that can divide `n²+n+A` is `A` itself — and the richness climbs monotonically:

| A | disc | first prime \| f | C |
|---|---|---|---|
| 3 | −11 | 3 | 1.02 |
| 5 | −19 | 5 | 1.88 |
| 11 | −43 | 11 | 3.26 |
| 17 | −67 | 17 | 4.17 |
| **41** | **−163** | **41** | **6.64** |

For `n²+n+41` we confirm `ω(q)=0` for **every** prime `q ≤ 37` — it dodges all of them — first divisible
at q=41. Non-Heegner `A` (9, 21, 33) are hit at q=3 and stay poor.

**Scope (honest).** None of this is new as theory — the mask is elementary; the local factors, the
constant C, and the Bateman–Horn asymptotic are classical; the class-number-one criterion is
Rabinowitsch's and the Heegner list is the Heegner–Stark–Baker theorem. What this adds is the 6N
geometric reading (mask rotation, wing sealing, the corridor as dodged annihilation channels) and a
high-precision measurement via the square-root sieve. **No claim is made about the infinitude of any
prime-rich polynomial.**

## Reproducing

```bash
pip install -r requirements.txt
cd code
python3 explore_n2a.py    # quick look: n²+a wing masks + richness (console)
python3 explore_euler.py  # quick look: Euler family + corridor (console)
python3 final_n2a.py      # n²+a representative table + C(a) scan -> data/n2a_*.csv   (~30 s)
python3 final_euler.py    # Euler family + corridor -> data/euler_*.csv, quad_summary.csv  (~15 s)
python3 makefigs_quad.py  # reads the CSVs -> figures/p32_fig1.pdf, p32_fig2.pdf
```

Paths resolve relative to the script (outputs → `../data`, `../figures`). The sieves use a Tonelli–
Shanks square root mod q and need ~0.1 GB RAM; single-threaded; NumPy 2.x compliant.

## Files

```
code/    explore_n2a.py  explore_euler.py  final_n2a.py  final_euler.py  makefigs_quad.py
data/    n2a_representative.csv   a, a_mod6, Q, right, left, channels, C_product, C_empirical
         n2a_Cscan.csv           a, a_mod6, C_product, wing_mask
         euler_family.csv        A, disc, class_number_1, first_prime_divisor, Q, right, left, C
         euler_corridor.csv      q, omega_A41, omega_A9
         quad_summary.csv        parameter, value
figures/ p32_fig1.pdf  p32_fig2.pdf
paper/   paper32.tex   paper32.pdf
```

All data files are plain CSV — openable in any text editor or spreadsheet.

## Citation

See `CITATION.cff`. The paper is archived on Zenodo (DOI in the citation file once minted).

## License

MIT (see `LICENSE`).
