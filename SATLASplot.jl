using RDatasets, DataFrames
using PyPlot
include("setupft.jl")
include("oichi2.jl")
include("oiplot.jl")
include("vfunc.jl")
include("hankel.jl")
#include("readoifits.jl")
sat1 = readtable("satlas_rsg/ld_satlas_surface.2t3300g-25m25", header = false, separator = ' ', nastrings=["NA", "na", "n/a", "missing"]);
names!(sat1, [Symbol("col$i") for i in 1:9]);
rename!(sat1, [:col1, :col6], [:mu, :flux])
fI = sat1[:6];
mu = sat1[:1];
rstar = 5;
r = rstar*(1-mu.^2);
f = diff(fI);
ra = diff(r);
a = f./ra;
b = fI[1:999].-(a.*r[1:999]);
a.*(hankel(x)) + b.*(visibility_ud(param, v_r))


#vis = (a.*(hankel(x)- b.*visibility_ud(param, v_r))
