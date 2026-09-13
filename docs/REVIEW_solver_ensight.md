# Review dossier: port of `solver.f90` and `ensight_cavity.f90` to C++

Files changed: `src/solver.cpp`, `src/ensight_cavity.cpp` (bodies of the stubs only; the
`extern "C"` hybrid wrappers were left untouched, see section D). Base branch:
`devin/cpp-port-skeleton`. Toolchain used: g++ 11.4.0, gfortran 11.4.0 (Ubuntu 22.04),
flags from `Makefile` (`-std=c++17 -O2 -Wall -Wextra` / `-fdefault-real-8 -O3`).

All line numbers below refer to the files as committed on `devin/port-solver-ensight`.

---

## A. `solver.f90` -> `src/solver.cpp`: storage order and diagonal offsets

Storage conventions (unchanged): `coef(1:n,1:m)` and `L(1:n,1:5)` are column-major
`FArray2<double>` views with Fortran bounds; `coef(i,j)` is the matrix entry
`(i, i+jcoef(j))`, `L(i,j)` is the factor entry `(i, i+Ldiag(j))`. `matgen_cavity.f90`
sets `jcoef = (0, 1, nx)`, hence `Ldiag = (0, 1, nx, nx-1, nx-2)`. Because the views keep the
Fortran bounds, every index expression is transcribed verbatim (no `-1` shifts anywhere).

Notation: F = `solver.f90` line, C = `src/solver.cpp` line.

### A.1 `ICCG2`

| F | Fortran | C | C++ | Justification |
|---|---------|---|-----|---------------|
| 43 | `do i=1,3; Ldiag(i)=jcoef(i); end do` | 128 | `for (i = 1; i <= 3; ++i) Ldiag(i) = jcoef(i);` | same 1-based range, `Ldiag`/`jcoef` are `FArray1<int>(1:5)`/`(1:m)` |
| 44 | `Ldiag(4)=Ldiag(3)-1` | 129 | `Ldiag(4) = Ldiag(3) - 1;` | literal |
| 45 | `Ldiag(5)=Ldiag(3)-2` | 130 | `Ldiag(5) = Ldiag(3) - 2;` | literal |
| 47 | `L=0.` | 132 | `L.fill(0.);` | whole-array assignment over all `n*5` elements |
| 49 | `do i=1,n` | 134 | `for (i = 1; i <= n; ++i)` | inclusive upper bound |
| 51 | `L(i,1)=coef(i,1)` | 135 | `L(i, 1) = coef(i, 1);` | literal |
| 52-54 | `do j=2,5; col=Ldiag(j); if (col<i) L(i,1)=L(i,1)-L(i-col,j)*L(i-col,j)*s(i-col)` | 136-138 | `for (j = 2; j <= 5; ++j) { col = Ldiag(j); if (col < i) L(i, 1) = L(i, 1) - L(i - col, j) * L(i - col, j) * s(i - col); }` | same guard `col<i` keeps `i-col >= 1`; product evaluated left-to-right in both languages (`(a*b)*c`) |
| 56 | `s(i)=1./L(i,1)` | 140 | `s(i) = 1. / L(i, 1);` | literal |
| 58-61 | `L(i,2)=coef(i,2); do j=4,5; col=Ldiag(j); if (col<i) L(i,2)=L(i,2)-L(i-col,j)*L(i-col,j-1)*s(i-col)` | 142-145 | `L(i, 2) = coef(i, 2); for (j = 4; j <= 5; ++j) { col = Ldiag(j); if (col < i) L(i, 2) = L(i, 2) - L(i - col, j) * L(i - col, j - 1) * s(i - col); }` | note the `j-1` second index (columns 3 and 4) kept as is |
| 64 | `L(i,3)=coef(i,3)` | 148 | `L(i, 3) = coef(i, 3);` | literal |
| 66-69 | `if (i>1) then; L(i,4)=-L(i-1,3)*L(i-1,2)*s(i-1); L(i,5)=-L(i-1,4)*L(i-1,2)*s(i-1)` | 150-153 | `if (i > 1) { L(i, 4) = -L(i - 1, 3) * L(i - 1, 2) * s(i - 1); L(i, 5) = -L(i - 1, 4) * L(i - 1, 2) * s(i - 1); }` | `L(i,5)` uses `L(i-1,4)` computed in the previous iteration, same sequential dependence |
| 73 | `do i=1,n; L(i,1)=sqrt(s(i)); end do` | 156 | `for (i = 1; i <= n; ++i) L(i, 1) = std::sqrt(s(i));` | the diagonal of `L` is overwritten by the scaling factor `1/sqrt(D)`; `std::sqrt(double)` = IEEE sqrt as gfortran's |
| 75-80 | `do j=1,3; col=Ldiag(j); do i=1,n-col; coef(i,j)=coef(i,j)*L(i,1)*L(i+col,1)` | 158-163 | `for (j = 1; j <= 3; ++j) { col = Ldiag(j); for (i = 1; i <= n - col; ++i) coef(i, j) = coef(i, j) * L(i, 1) * L(i + col, 1); }` | symmetric scaling of the 3 stored diagonals; upper bound `n-col` keeps `i+col <= n`; only the first 3 columns of `Ldiag` are used here as in Fortran |
| 82-97 | second factorisation pass on the scaled matrix (same statements as 58-69 without the `s(...)` factors) | 165-178 | same, without `s(...)` | literal; `L(i,1)` is *not* recomputed (stays `sqrt(s)`), exactly as in Fortran |
| 99 | `norm0=norm2(b,n)` | 180 | `norm0 = norm2(b, n);` | computed on the **unscaled** `b`, before line 103 |
| 101 | `do i=1,n; s(i)=1./s(i); end do` | 182 | `for (i = 1; i <= n; ++i) s(i) = 1. / s(i);` | `s` becomes the (unscaled) pivot `D` again, used as weight in the norm |
| 103 | `do i=1,n; b(i)=b(i)*L(i,1); end do` | 184 | `for (i = 1; i <= n; ++i) b(i) = b(i) * L(i, 1);` | rhs scaling (the Fortran modifies `b`; so does the port) |
| 109-110 | `call matmul_ell(p,coef,jcoef,x,n,m); r=b-p` | 190-191 | `matmul_ell(p, coef, jcoef, x, n, m); for (...) r(i) = b(i) - p(i);` | literal |
| 112-115 | `call resLU_ell5(r2,L,Ldiag,r,n); p=r2; nu=scal(r,r2,n)` | 193-196 | same | `p=r2` array copy written as loop |
| 117-119 | `norm=0.; do i=1,n; norm=norm+r(i)*r(i)*s(i); end do; norm=sqrt(norm)/norm0` | 198-200 | same | `(r*r)*s` left-to-right; the convergence quantity is `sqrt(sum r_i^2 D_i) / ||b_unscaled||_2`, **not** `norm2(r)/norm2(b)` literally - the port keeps the source's exact expression |
| 121-123 | `j=0; do while (norm>eps.and.j.lt.itmax)` | 202-204 | `j = 0; while (norm > eps && j < itmax)` | strict `>` and `<` kept |
| 126-145 | CG body: `q=A p; alpha=nu/scal(p,q,n); x+=alpha p; r-=alpha q; r2=M^-1 r; mu=scal(r,r2,n); beta=mu/nu; p=r2+beta p; nu=mu; norm=...` | 207-226 | same statements, same order | two separate `do i` loops for `x` and `r` kept separate (no fusion) |
| 153 | `do i=1,n; x(i)=x(i)*L(i,1); end do` | 233 | `for (i = 1; i <= n; ++i) x(i) = x(i) * L(i, 1);` | solution unscaling |
| 155-159 | `if(j.ge.itmax) print*,'non convergence apres =',j,norm; else print*,' Nombre ... =',j` | 237-241 | `if (j >= itmax) printf(" non convergence apres = norm=%24.17G it=%11d\n", norm, j); else printf("  Nombre Iterations ICCG2 ( Fill-In 2 ) = %11d\n", j);` | **deviation, see D.1**: in the non-converged branch the Fortran prints `j` *then* `norm`; the port prints `norm` first so that `j` is the last token as required by the extraction rule |

### A.2 `resLU_ell5` (forward sweep with `L`, backward sweep with `L^T`)

| F | Fortran | C | C++ | Justification |
|---|---------|---|-----|---------------|
| 180-183 | `col1=Ldiag(2); col2=Ldiag(5); col3=Ldiag(4); col4=Ldiag(3)` | 23-26 | identical | offsets in increasing order: `1, nx-2, nx-1, nx`; columns of `L` used with each: `2, 5, 4, 3` |
| 185-187 | `do i=1,col1; x(i)=y(i)` | 28-30 | `for (i = 1; i <= col1; ++i) x(i) = y(i);` | no lower diagonal reaches these rows |
| 189-191 | `do i=col1+1,col2; x(i)=y(i)-L(i-col1,2)*x(i-col1)` | 32-34 | literal | `L(i-col1,2)` is entry `(i-col1, i)` = column `i` of `L^T`, i.e. row `i` of `L` |
| 193-195 | `do i=col2+1,col3; ...-L(i-col2,5)*x(i-col2)` | 36-38 | literal | adds offset `nx-2` (column 5) |
| 197-200 | `do i=col3+1,col4; ...-L(i-col3,4)*x(i-col3)` | 40-43 | literal | adds offset `nx-1` (column 4) |
| 202-205 | `do i=col4+1,n; ...-L(i-col4,3)*x(i-col4)` | 45-48 | literal | adds offset `nx` (column 3); subtraction chain left-to-right in both |
| 207-209 | `do i=n-col1,n-col2+1,-1; x(i)=x(i)-L(i,2)*x(i+col1)` | 50-52 | `for (i = n - col1; i >= n - col2 + 1; --i)` | Fortran `do i=a,b,-1` runs `a` down to `b` inclusive |
| 211-213 | `do i=n-col2,n-col3+1,-1; ...-L(i,5)*x(i+col2)` | 54-56 | literal | |
| 215-217 | `do i=n-col3,n-col4+1,-1; ...-L(i,4)*x(i+col3)` | 58-60 | literal | |
| 219-222 | `do i=n-col4,1,-1; ...-L(i,3)*x(i+col4)` | 62-65 | literal | |

Note: the Fortran declares `col5, col6` but never uses them; they are omitted in C++.

### A.3 `matmul_ell` (symmetric product from the stored upper diagonals)

| F | Fortran | C | C++ | Justification |
|---|---------|---|-----|---------------|
| 244-246 | `do i=1,n; x(i)=coef(i,1)*y(i)` | 73-75 | literal | diagonal `jcoef(1)=0` |
| 248-250 | `do j=2,m; col=jcoef(j); do i=1,n-col` | 77-79 | `for (j = 2; j <= m; ++j) { col = jcoef(j); for (i = 1; i <= n - col; ++i)` | `n-col` keeps `i+col<=n` |
| 251 | `x(i)=x(i)+coef(i,j)*y(i+col)` | 80 | literal | upper part |
| 252 | `x(i+col)=x(i+col)+coef(i,j)*y(i)` | 81 | literal | symmetric (lower) part, same statement order inside the loop so the accumulation order of `x` is identical |

### A.4 `scal`, `norm1`, `norm2`

Lines F 261-310 -> C 86-112: `res=0.; do i=1,n; res=res+x(i)*y(i)` / `x(i)*x(i)`; `norm2` applies
`sqrt` at the end. Same accumulation order (sequential, no pairwise/vectorised reduction:
neither `-O2` g++ nor `-O3` gfortran reassociates FP sums without `-ffast-math`).
`norm1` is ported although it is unused (as in the Fortran).

---

## B. `ensight_cavity.f90` -> `src/ensight_cavity.cpp`: byte layout

Fortran: `open(form='UNFORMATTED', access='direct', recl=reclength)` + one `write(rec=1)`. gfortran
writes the record payload with no markers and pads the record with NUL bytes to `recl`. In C++
(`DirectRecord`, `src/ensight_cavity.cpp` 31-64) the items are appended to a byte buffer in the
order of the Fortran `write` list, the buffer is zero-padded to `reclength` and written with a
single `ofstream::write` (file truncated). 80-byte strings are blank padded (`fstring80`,
l. 24), integers are 4-byte native little-endian (`put_int4`), reals are `float(x)` =
`real(sngl(x),4)` (`put_real4`). Run values: `nx=ny=11`, `nz=1`, `imin=jmin=kmin=1`,
`imax=jmax=11`, `kmax=1`.

### B.1 `EnsightOutput.geo` (`WriteRectEnsightGeo`, F 94-144, C 87-128)

`reclength = 80*8 + 4*(4 + isize + jsize + ksize) = 640 + 4*27 = 748` (F 126 / C 108). The
payload is exactly 748 bytes, so there is no padding; file size 748 in both outputs.

| offset (hex) | size | content | Fortran write item (F 129-140) | C++ (line) |
|---|---|---|---|---|
| 0x000 | 80 | `C Binary` + blanks | `binary_form` | `put_string80(binary_form)` 111 |
| 0x050 | 80 | `Ensight Model Geometry File Created by ` | `file_description1` | 112 |
| 0x0a0 | 80 | `WriteRectEnsightGeo Routine` | `file_description2` | 113 |
| 0x0f0 | 80 | `node id off` | `node_id` | 114 |
| 0x140 | 80 | `element id off` | `element_id` | 115 |
| 0x190 | 80 | `part` | `part` | 116 |
| 0x1e0 | 4 | int `1` | `npart` | `put_int4(npart)` 117 |
| 0x1e4 | 80 | `Sortie Ensight` | `description_part` | 118 |
| 0x234 | 80 | `block rectilinear` | `blocck` | 119 |
| 0x284 | 4 | int `11` | `isize` | 120 |
| 0x288 | 4 | int `11` | `jsize` | 121 |
| 0x28c | 4 | int `1` | `ksize` | 122 |
| 0x290 | 44 | 11 floats `x1(1..11)` | `(real(sngl(x1(i)),4),i=imin,imax)` | loop 123 |
| 0x2bc | 44 | 11 floats `x2(1..11)` | `(real(sngl(x2(j)),4),j=jmin,jmax)` | loop 124 |
| 0x2e8 | 4 | 1 float `0.0` | `(real(sngl(x3(k)),4),k=kmin,kmax)` | loop 125 |
| 0x2ec | - | end of file (748 bytes) | | |

`xxd` head of the reference `reference/upwind/EnsightOutput.geo`:

```
00000000: 4320 4269 6e61 7279 2020 2020 2020 2020  C Binary
00000010: 2020 2020 2020 2020 2020 2020 2020 2020
...
00000050: 456e 7369 6768 7420 4d6f 6465 6c20 4765  Ensight Model Ge
00000060: 6f6d 6574 7279 2046 696c 6520 4372 6561  ometry File Crea
00000070: 7465 6420 6279 2020 2020 2020 2020 2020  ted by
...
00000190: 7061 7274 2020 2020 2020 2020 2020 2020  part
...
000001e0: 0100 0000 536f 7274 6965 2045 6e73 6967  ....Sortie Ensig
000001f0: 6874 2020 2020 2020 2020 2020 2020 2020  ht
...
00000230: 2020 2020 626c 6f63 6b20 7265 6374 696c      block rectil
00000240: 696e 6561 7220 2020 2020 2020 2020 2020  inear
...
00000280: 2020 2020 0b00 0000 0b00 0000 0100 0000      ............
00000290: 8c2e 3a3d e9a2 0b3e 2fba 683e bae8 a23e  ..:=...>/.h>...>
000002a0: 5d74 d13e 0000 003f d145 173f a38b 2e3f  ]t.>...?.E.?...?
000002b0: 74d1 453f 4617 5d3f 175d 743f 8c2e 3a3d  t.E?F.]?.]t?..:=
000002c0: e9a2 0b3e 2fba 683e bae8 a23e 5d74 d13e  ...>/.h>...>]t.>
000002d0: 0000 003f d145 173f a38b 2e3f 74d1 453f  ...?.E.?...?t.E?
000002e0: 4617 5d3f 175d 743f 0000 0000            F.]?.]t?....
```

`xxd` of the C++ output (`build/run_bc/EnsightOutput.geo`, hybrid `PORTED="solver ensight_cavity"`,
`-DNS_BOUNDS_CHECK`, scheme 1) - identical (`cmp` silent, 748 bytes):

```
00000000: 4320 4269 6e61 7279 2020 2020 2020 2020  C Binary
00000010: 2020 2020 2020 2020 2020 2020 2020 2020
...
00000280: 2020 2020 0b00 0000 0b00 0000 0100 0000      ............
00000290: 8c2e 3a3d e9a2 0b3e 2fba 683e bae8 a23e  ..:=...>/.h>...>
000002a0: 5d74 d13e 0000 003f d145 173f a38b 2e3f  ]t.>...?.E.?...?
000002b0: 74d1 453f 4617 5d3f 175d 743f 8c2e 3a3d  t.E?F.]?.]t?..:=
000002c0: e9a2 0b3e 2fba 683e bae8 a23e 5d74 d13e  ...>/.h>...>]t.>
000002d0: 0000 003f d145 173f a38b 2e3f 74d1 453f  ...?.E.?...?t.E?
000002e0: 4617 5d3f 175d 743f 0000 0000            F.]?.]t?....
```

### B.2 `velocityNNNNNN.vec` / `omegaNNNNNN.scl` / `pressNNNNNN.scl` (`WriteEnsightVar`, F 159-198, C 142-173)

`reclength = 80*3 + 4*(1 + isize*jsize*ksize*ndv)` (F 182 / C 158):
`.vec` (ndv=3): `240 + 4*(1+363) = 1696`; `.scl` (ndv=1): `240 + 4*(1+121) = 728`.
Payload = recl in both cases, so no padding; file sizes match (1696 / 728).

| offset (hex) | size | content | Fortran write item (F 188-194) | C++ (line) |
|---|---|---|---|---|
| 0x000 | 80 | `velocity000200.vec` (the **file name**, incl. extension) | `VarFileName` | `put_string80(VarFileName)` 161 |
| 0x050 | 80 | `part` | `part` | 162 |
| 0x0a0 | 4 | int `1` | `npart` | `put_int4(npart)` 163 |
| 0x0a4 | 80 | `block rectilinear` | `blocck` | 164 |
| 0x0f4 | 4*121 | component 1: `var(1,i,j,1)`, i fastest, then j | `((((real(sngl(var(m,i,j,k)),4),i=imin,imax),j=jmin,jmax),k=kmin,kmax),m=1,ndv)` | nested loops `m` (outer) / `k` / `j` / `i` (inner), 166-170 |
| 0x2d8 | 4*121 | component 2 (`.vec` only) | idem, `m=2` | |
| 0x4bc | 4*121 | component 3 = 0.0 (`.vec` only) | idem, `m=3` | |
| 0x6a0 / 0x2d8 | - | end of file (1696 / 728 bytes) | | |

The `var` argument is a `std::vector<double>` with the Fortran `(ndv,nx,ny,nz)` column-major
layout: element `(m,i,j,k)` at `(m-1) + ndv*((i-1) + nx*((j-1) + ny*(k-1)))` (C 170); the
fill loops in `write_result_ensight` use the same formula (`vel`/`sca` lambdas, C 251-254).

`xxd` head of the reference `reference/upwind/velocity000200.vec`:

```
00000000: 7665 6c6f 6369 7479 3030 3032 3030 2e76  velocity000200.v
00000010: 6563 2020 2020 2020 2020 2020 2020 2020  ec
...
00000050: 7061 7274 2020 2020 2020 2020 2020 2020  part
...
000000a0: 0100 0000 626c 6f63 6b20 7265 6374 696c  ....block rectil
000000b0: 696e 6561 7220 2020 2020 2020 2020 2020  inear
...
000000f0: 2020 2020 c2a3 59b9 f9d1 1fba c4da 7cba      ..Y.......|.
00000100: 9c57 a2ba 6ff3 b8ba 42aa c0ba abf6 b8ba  .W..o...B.......
00000110: 165d a2ba c8e6 7cba 2fdb 1fba 7db1 59b9  .]....|./...}.Y.
...
00000690: 0000 0000 0000 0000 0000 0000 0000 0000  ................   (third component = 0)
```

`xxd` of the C++ output `build/run_bc/velocity000200.vec` - identical (`cmp` silent, 1696 bytes):

```
00000000: 7665 6c6f 6369 7479 3030 3032 3030 2e76  velocity000200.v
00000010: 6563 2020 2020 2020 2020 2020 2020 2020  ec
...
000000a0: 0100 0000 626c 6f63 6b20 7265 6374 696c  ....block rectil
000000b0: 696e 6561 7220 2020 2020 2020 2020 2020  inear
...
000000f0: 2020 2020 c2a3 59b9 f9d1 1fba c4da 7cba      ..Y.......|.
00000100: 9c57 a2ba 6ff3 b8ba 42aa c0ba abf6 b8ba  .W..o...B.......
...
00000690: 0000 0000 0000 0000 0000 0000 0000 0000  ................
```

### B.3 `EnsightOutput.case` (`EnsightCase`, F 212-252, C 183-229) - text

| Fortran | output line(s) | C++ (line) |
|---|---|---|
| format 10: `'FORMAT',/,'type: ensight gold',//,'GEOMETRY',/,'model:<TAB>',A,//,'VARIABLE'` | `FORMAT` / `type: ensight gold` / `` / `GEOMETRY` / `model:<TAB>EnsightOutput.geo` / `` / `VARIABLE` | 193-199 (the literal tab of the format string is reproduced as `\t`) |
| format 15/25 (`nfile>1`): `'scalar per node: ',A,'   ',A` with `trim(name)//'******.scl'` | `vector per node: velocity   velocity******.vec`, `scalar per node: omega   omega******.scl`, `scalar per node: press   press******.scl` | 203-215 (the `nfile==1` branch without `******` is ported too) |
| format 45 (`nfile>1`): `/,'TIME            ',/,'time set: 1     ',/,'number of steps:',i4,/,'filename start number:',i10/,'filename increment:',i10/,'time values: '` | `` / `TIME            ` / `time set: 1     ` / `number of steps: 100` / `filename start number:       200` / `filename increment:       200` / `time values: ` | 219-225; `fint(w,v)` (C 67) = Fortran `Iw` (right-justified, `*`-filled on overflow) |
| `write(FileUnit,'(i5)') (j*isto,j=1,nfile)` | one `%5d` value per line, `  200` ... `20000` | 227 |

Trailing blanks of the format literals (`TIME            `, `time set: 1     `, `time values: `)
are kept. The `.case` file is rewritten identically at every `isto` step, as in Fortran.

### B.4 `write_result_ensight` (F 1-82, C 234-310)

* file names: `'(A,I6.6)'` -> `"%06d"` (C 279-282); the `nstep/isto == 1` branch (no step
  suffix, C 274-277) is ported as well.
* `zEnsight` is a local `FArray1<double>(1:1)` set to `0.` (C 245-246, 266).
* The omega call passes `jmin` in the `kmin` slot exactly like the Fortran (C 303); both are 1.
* `nbImages` is computed and unused, as in Fortran.

---

## C. Checks run (all on this branch, from a clean `build/`)

1. `make clean && make check-hybrid PORTED="solver"`
   ```
   == upwind (.../build/run_upwind)
   all files match the reference
   == centred (.../build/run_centred)
   all files match the reference
   ```
   exit 0.
2. `make clean && make check-hybrid PORTED="ensight_cavity"` - same two
   `all files match the reference` lines, exit 0.
3. `make clean && make check-hybrid PORTED="solver ensight_cavity"` - same
   two `all files match the reference` lines, exit 0. Additionally
   `cmp build/run_upwind/EnsightOutput.geo reference/upwind/EnsightOutput.geo` and
   `cmp build/run_upwind/velocity000200.vec reference/upwind/velocity000200.vec` are silent.
4. Bounds-checked build, one scheme:
   `make clean && make hybrid PORTED="solver ensight_cavity" CXXFLAGS="-std=c++17 -O2 -Wall -Wextra -DNS_BOUNDS_CHECK"`,
   then `cd build/run_bc && echo 1 | NS_DUMP_FIELDS=1 ../../NS_lid_driven_cavity_hybrid`
   -> exit 0, empty stderr (no `array bounds error`), 20000 `Nombre` lines;
   `python3 tools/compare_reference.py build/run_bc reference/upwind` -> `all files match the
   reference`, exit 0.
5. `iccg2_iterations.txt` of both reference schemes contains only the values 10 (2 steps) and
   11 (19998 steps); the C++ runs reproduce them exactly (checked by 1-4, exact text compare).
6. `make cpp` (pure C++ executable, remaining units still stubs) compiles both files
   with `-Wall -Wextra` without warnings.
7. No file under `reference/`, `tools/`, `Makefile`, `src/variables.*`, `src/routines.hpp` was
   modified; no tolerance changed.

---

## D. Points for the reviewer / things I was unsure about

1. **Non-convergence message order (deliberate deviation).** Fortran prints
   `'non convergence apres =', j, norm` (norm is the last token). The task/PORTING.md require the
   iteration count to be the LAST token of that line, so the C++ prints
   `non convergence apres = norm=<norm> it=<j>`. It never triggers in the reference runs (all
   20000 steps converge in 10-11 iterations), so the check does not exercise it. If the reviewer
   prefers strict fidelity to the Fortran text instead, swap the two arguments (C 238).
2. **stdout buffering in the hybrid build.** The C++ `printf` lines share stdout with
   libgfortran's own buffer. To guarantee the `Nombre ...` line is emitted whole (so its last
   token is `j` regardless of how the two buffers interleave), `std::fflush(stdout)` is called
   after each ICCG2 print (C 123, 242). This only affects diagnostic output.
3. **Convergence criterion.** It is *not* literally `norm2(r)/norm2(b)`: the source uses
   `sqrt(sum r_i^2 * s_i) / norm2(b_unscaled)` with `s = D` (the pivots recovered at F 101), and
   `norm0` is computed **before** `b` is scaled (F 99 vs 103). The port keeps this exact
   sequence; do not "simplify" it.
4. **`L(i,1)` semantics.** After F 73 the first column of `L` no longer holds the pivot but
   `1/sqrt(D)`; the second factorisation pass (F 82-97) does not recompute it. The port mirrors
   this; the reviewer may want to confirm this is the intended algorithm (it is what the
   reference oracle was produced with).
5. **Direct-access file creation.** `WriteEnsightVar` opens without `status='replace'`. With a
   pre-existing longer file gfortran would keep the extra tail bytes; the C++ always truncates.
   Not observable in the checks (fresh run directories) but a behavioural difference on reruns
   in a dirty directory (the C++ behaviour is arguably the intended one).
6. **Padding.** For the actual sizes (`nx=ny=11`, `nz=1`) every record payload equals `recl`, so
   the NUL padding path of `DirectRecord::write_file` is never exercised by the checks. It
   follows the gfortran behaviour documented in PORTING.md.
7. **Hybrid wrappers.** `iccg2_` and `write_result_ensight_` were checked against the Fortran
   call sites (`main.f90` 36 and 62: `coef(ndim,mdim)`, `jcoef(mdim)`, `L(ndim,5)`, `Ldiag(5)`,
   vectors `(ndim)`; `xx(nx)`, `yy(ny)`, `uc/vc/rot/pre(nx,ny)`, scalars by reference) and found
   correct; **no wrapper was modified**.
8. **Endianness / float format.** `put_int4`/`put_real4` write the native representation
   (x86-64: little-endian, IEEE binary32), which is what gfortran writes on the same platform.
   The Ensight header says `C Binary`, so this matches the file's declared convention.
9. **Local variable `nu` in `ICCG2`** shadows the global `ns::nu` (viscosity) inside the function
   only, same as the Fortran local shadows the module variable. No `-Wshadow` warning is enabled
   in `Makefile`.
