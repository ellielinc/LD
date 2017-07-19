using PyPlot
include("oiplot.jl")

inty = zeros(1024,1024)
center = [1025/2, 1025/2]
nx = 1024
#x = reshape([(div(i-1,nx)+1) for i=1:nx*nx], nx, nx)-(nx+1)/2
#y = reshape([(div(i-1,nx)+1) for i=1:nx*nx], nx, nx)-(nx+1)/2
r = linspace(0, 300, 20)
t = linspace(0, 6.283185307179586, 20)
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

for i=1:length(x)
    for j=1:length(y)
    r=sqrt(x[i].^2+y[j].^2)
  if r>=200
  x[i]=0
  y[j]=0
  else
  x[i]=1
  y[j]=1
  end
end
end

#=
mu = sqrt(1-(r./rstar).^2)
alpha = 0.89
f = mu.^alpha
flux = vec(r)
#imdisp(flux)
#theta = atan(x./y)
#scatter(r, theta)
#make into a circle
=#
