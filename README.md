# BlockMatrices

[![Build Status](https://travis-ci.org/mfalt/BlockMatrices.jl.svg?branch=master)](https://travis-ci.org/mfalt/BlockMatrices.jl)

[![Coverage Status](https://coveralls.io/repos/mfalt/BlockMatrices.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/mfalt/BlockMatrices.jl?branch=master)

[![codecov.io](http://codecov.io/github/mfalt/BlockMatrices.jl/coverage.svg?branch=master)](http://codecov.io/github/mfalt/BlockMatrices.jl?branch=master)

Implements Block Matrices as an efficient alternative to concatenating matrices of different types.

The package supplies the sized `ZeroMatrix` and `IdentityMatrix`.

Supported methods: `A_mul_B!, At_mul_B!, Ac_mul_B!`, i.e. also `A*b, A'b, A.'b'`.

### Examples

Time and allocation comparison for large and sparse blocks:
```julia
A = randn(2000,2000);
B = spzeros(2000,20000);
C = IdentityMatrix{Float64}(2.0,2000);
D = sprandn(2000,20000,0.001);

#Create [A B; C D] BlockMatrix
M1 = BlockMatrix(((A,B),(C,D)));          # 16.821 μs (135 allocations: 176.58 KiB)
M2 = [A B; 2I D];       # 117.897 ms (104 allocations: 215.84 MiB)
M3 = full([A B; 2I D]); # 200.390 ms (106 allocations: 887.23 MiB)

b = randn(size(M1,2));
y = randn(size(M1,1));

A_mul_B!(y, M1, b); #  1.692 ms (9 allocations: 432 bytes)
A_mul_B!(y, M2, b); #  4.049 ms (0 allocations: 0 bytes)
A_mul_B!(y, M3, b); # 32.707 ms (0 allocations: 0 bytes)
```

For smaller blocks:
```julia
A = randn(20,20);
B = spzeros(20,20);
C = IdentityMatrix{Float64}(2.0,20);
D = sprandn(20,20,0.1);

M1 = BlockMatrix(((A,B),(C,D)));  # 17.069 μs (113 allocations: 4.81 KiB)
M2 = [A B; 2I D];    # 17.565 μs (81 allocations: 31.30 KiB)
M3 = full([A B; 2I D]) # 19.735 μs (64 allocations: 45.48 KiB)

b = randn(size(M1,2));
y = randn(size(M1,1));

A_mul_B!(y, M1, b); # 440.182 ns (9 allocations: 432 bytes)
A_mul_B!(y, M2, b); # 464.604 ns (0 allocations: 0 bytes)
A_mul_B!(y, M3, b); # 233.207 ns (0 allocations: 0 bytes)
```
