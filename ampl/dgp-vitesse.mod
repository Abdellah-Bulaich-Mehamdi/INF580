# DEFINE THE NUMBER OF ELEMENTS AND THE NEEDED PARAMETERS
param n integer;
param k integer;
param time integer;
param eps default 1e-3;


# DEFINE THE SET OF INDICES AND THE NEEDED PARAMETERS
set T := 1..time;
set P := 1..(time-1);
set K := 1..k;
set S := 1..n;
set E := {T, S, S};
param d{E};
param vit{P, S};

# DEFINE THE VARIABLES
var x{T, S, K};

# PAGE 185
# DEFINE THE OBJECTIVE FUNCTION
minimize f: sum{(t, u, v) in E: u < v} sum{l in K} (x[t, u, l] - x[t, v, l])^2;

# DEFINE THE CONSTRAINTS
s.t. dist{(t, u, v) in E: u < v}: sum{l in K} (x[t, u, l] - x[t, v, l])^2  >= d[t, u, v]^2;
s.t. ones{t in T, l in K}       : sum{u in S}  x[t, u, l]== 0;
s.t. vitu{t in P, u in S}       : sum{l in K} (x[t, u, l] - x[t+1, u, l])^2  <= vit[t, u]^2 +eps;
s.t. vitl{t in P, u in S}       : sum{l in K} (x[t, u, l] - x[t+1, u, l])^2  >= vit[t, u]^2 -eps;
