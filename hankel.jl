using PyCall
sp=pyimport("scipy.special")
#hankel = sp[:struve](order, x)

function phi(x)
 phix = (pi*x./2).*(besselj1(x).*sp[:struve](0, x).-besselj(0,x).*sp[:struve](1, x))
end

function hankel(x)
    hankel= x.^2.*besselj1(x)-phi(x)
end
