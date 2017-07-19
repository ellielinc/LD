using RDatasets, DataFrames
using PyPlot
include("oiplot.jl")

satfile = readtable("satlas_rsg/ld_satlas_surface.2t3300g-25m25", header = false, separator = ' ', nastrings=["NA", "na", "n/a", "missing"]);
names!(satfile, [Symbol("col$i") for i in 1:9]);
rename!(satfile, [:col1, :col6], [:u, :v]);
fI = satfile[:6];
m = satfile[:1];
rstar = 20;
r = rstar*(1-m.^2);
#t = linspace(0, 6.283185307179586, 190)
t = linspace(0, 1.5707963267948966, 135)
x = zeros(length(t)*length(r));
for i = 1:length(t)
  for u = 1:length(r)
    j = t[i]
  x[u*i] = r[u].*cos(j)
  end
end

y = zeros(length(t)*length(r));
for i = 1:length(t)
  for u = 1:length(r)
  j = t[i]
  y[u*i] = r[u].*sin(j)
  end
end
scatter(x, y)
r = sqrt(x.^2+y.^2);
rstar = 20;
#mu = sqrt(1-(r/rstar).^2)
