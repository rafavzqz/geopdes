function s = rectangular_plate_analytical_solution (a, b, h, E, nu, q0)
%
% INPUT
%   a:  length of the rectangular section
%   b:  height of the rectangular section
%   h:  thickness of the plate
%   E:  Young modulus
%   nu: Poisson ratio
%   q0: distributed load
%
% OUTPUT
%   displ: maximum displacement

x=a/2;
y=b/2;
%pi = %pi;
D=E*h^3/(12*(1-nu^2));
ks = 16*q0*a^4*b^4/(pi^6*D);
km = 16*q0*a^4*b^4/(pi^4);

s=0;
mx=0;
my=0;
n = 1;
nMax = 1001;
while n<=nMax
    m = 1;
    while m<=nMax
        den = (n^2*a^2 + m^2*b^2)^2;
        sinA = sin(m*pi*x/a);
        sinB = sin(n*pi*y/b);
        knm = sinA*sinB/(m*n*den);
        s = s + ks*knm;
        mx = mx + km*knm*((m/a)^2 + nu*(n/b)^2);
        my = my + km*knm*(nu*(m/a)^2 + (n/b)^2);
        m = m+2;
    end;
    n = n+2;
end

x=0;
y=0;
kmxy = -(1-nu)*16*q0*a^3*b^3/(pi^4);
mxy = 0;
n=1;
while n<=nMax
    m = 1;
    while m<=nMax
        den = (n^2*a^2 + m^2*b^2)^2;
        cosA = cos(m*pi*x/a);
        cosB = cos(n*pi*y/b);
        mxy = mxy + kmxy*cosA*cosB/den;
        m = m+2;
    end;
    n = n+2;
end

