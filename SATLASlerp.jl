using RDatasets, DataFrames
using PyPlot
include("setupft.jl")
include("oichi2.jl")
include("readoifits.jl")

sat1 = readtable("satlas_rsg/ld_satlas_surface.2t3300g-25m25", header = false, separator = ' ', nastrings=["NA", "na", "n/a", "missing"]);
names!(sat1, [Symbol("col$i") for i in 1:9]);
rename!(sat1, [:col1, :col6], [:mu, :flux]);
f_satlas = convert(Array,sat1[:6]);
push!(f_satlas,1.0);
mu_satlas = convert(Array,sat1[:1]);
push!(mu_satlas,1.0);
data = read_oifits("AlphaCenA.oifits");

function complex_vis_satlas_img(diameter, f_satlas, mu_satlas; npix=256, pixsize=0.05)
x=repmat(collect(1:npix),1,npix)-(npix+1)/2;
y=x';
r=sqrt(x.^2+y.^2);
rstar = 0.5*diameter/pixsize # rstar in pixels
mask=find(r.>rstar);
r[mask]=rstar;
#new mu values to interpolate:
mu_star = sqrt.(1-(r/rstar).^2);
indx_lo = vec(Int.((floor(mu_star.*1000)+1)))
indx_hi = vec(Int.((floor(mu_star.*1000)+2)))
mu_lo = mu_satlas[indx_lo]
mu_hi = mu_satlas[indx_hi]
f_lo = f_satlas[indx_lo]
f_hi = f_satlas[indx_hi]
f = f_lo+(vec(mu_star)-mu_lo).*((f_hi-f_lo)./(mu_hi-mu_lo))
#f = reshape(f, (npix,npix)), imshow(f)
dft = setup_ft(data, npix, pixsize * (pi / 180.0) / 3600000.0)
#ERROR:type #data has no field nuv
cvis = image_to_cvis(f, dft)
return cvis
end


cvis= complex_vis_satlas_img(2.5,f_satlas, mu_satlas)


#loop for different pixel sizes and different SATLAS files
