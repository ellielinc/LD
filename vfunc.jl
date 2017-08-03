#uniform disk
function visibility_ud(param, v_r)
theta = param[1]/2.0626480624709636e8;
V = 2*(besselj1(pi*theta*v_r))./(pi*theta*v_r)
end

#power law
function visibility_ldpow(param, v_r)
maxk=200;
theta = param[1]/2.0626480624709636e8;
f=k->gamma(param[2]/2+2)/(gamma(param[2]/2+k+2)*gamma(k+1))*(-0.25*(pi*theta*v_r).^2).^k;
V = map(f,collect(0:maxk));
sum(V)
end

#quadratic law
function visibility_ldquad(param,v_r)
theta = param[1]/2.0626480624709636e8;
zeta = (pi*v_r*theta);
V = ((1.0-param[2]-param[3])*(besselj1(zeta)./zeta)+((theta+2*param[3])/sqrt(2/pi))*(
(sqrt(2 ./(pi*zeta)).*((sin(zeta)./zeta)-cos(zeta)))./zeta.^(3/2))-2*param[3]*
(besselj(2, zeta)./zeta.^2))./(0.5-param[2]/6-param[3]/12)
end
