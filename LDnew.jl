using PyPlot
include("oiplot.jl")
include("readoifits.jl")
include("setupft.jl")
include("oichi2.jl")
include("vfunc.jl")

nx = 256
#tic();
y = reshape([(div(i-1,nx)+1) for i=1:nx*nx], nx, nx)-(nx+1)/2;
x = reshape([(mod(i-1,nx)+1) for i=1:nx*nx], nx, nx)-(nx+1)/2;
#toc();

r=sqrt(x.^2+y.^2);
rstar= 100;
disc_indx = find(r.<=rstar);
mu = zeros(size(r));
mu[disc_indx] = sqrt(1-(r[disc_indx]/rstar).^2);
alpha = 0.89;
f = mu.^alpha;
flux = vec(f);
imdisp(flux)
#savefig("LDsamp.jpg")
# step 2
image = copy(flux);
pixellation = 0.01; # in mas/pixel
scale_rad = pixellation * (pi / 180.0) / 3600000.0;

#how to save img as oifits
# open data file
oifitsfile = "2004-data1.oifits";
data = read_oifits(oifitsfile);

#SATLAS?

# setup Fourier transform matrix
dft = setup_ft(data, nx);
# compute Fourier transform of the image
cvis_image = image_to_cvis(image, dft)

u=data.uv[1,:];
v=data.uv[2,:];
ruv=sqrt(u.^2+v.^2);
cvis_analytic = visibility_ldpow([2.0, 0.89], ruv)
v2_model = cvis_to_v2(cvis_analytic, data.indx_v2);
v2_model2 = cvis_to_v2(cvis_image, data.indx_v2);
opt = maximum(abs(v2_model-v2_model2)./abs(v2_model))*100.
println("maximum = $opt")
