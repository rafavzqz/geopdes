function [bp,wf]=grule(n)
%
% [bp,wf]=grule(n)
%  This function computes Gauss base points and weight factors
%  using the algorithm given by Davis and Rabinowitz in 'Methods
%  of Numerical Integration', page 365, Academic Press, 1975.
%
bp=zeros(n,1); wf=bp; iter=2; m=fix((n+1)/2); e1=n*(n+1);
mm=4*m-1; t=(pi/(4*n+2))*(3:4:mm); nn=(1-(1-1/n)/(8*n*n));
xo=nn*cos(t);
for j=1:iter
   pkm1=1; pk=xo;
   for k=2:n
      t1=xo.*pk; pkp1=t1-pkm1-(t1-pkm1)/k+t1;
      pkm1=pk; pk=pkp1;
   end
   den=1.-xo.*xo; d1=n*(pkm1-xo.*pk); dpn=d1./den;
   d2pn=(2.*xo.*dpn-e1.*pk)./den;
   d3pn=(4*xo.*d2pn+(2-e1).*dpn)./den;
   d4pn=(6*xo.*d3pn+(6-e1).*d2pn)./den;
   u=pk./dpn; v=d2pn./dpn;
   h=-u.*(1+(.5*u).*(v+u.*(v.*v-u.*d3pn./(3*dpn))));
   p=pk+h.*(dpn+(.5*h).*(d2pn+(h/3).*(d3pn+.25*h.*d4pn)));
   dp=dpn+h.*(d2pn+(.5*h).*(d3pn+h.*d4pn/3));
   h=h-p./dp; xo=xo+h;
end
bp=-xo-h;
fx=d1-h.*e1.*(pk+(h/2).*(dpn+(h/3).*(d2pn+(h/4).*(d3pn+(.2*h).*d4pn))));
wf=2*(1-bp.^2)./(fx.*fx);
if ( (m+m) > n )
	bp(m)=0; 
end
if ( ~ ((m+m) == n) )
	m=m-1;
end
jj=1:m; n1j=(n+1-jj); bp(n1j)=-bp(jj); wf(n1j)=wf(jj);
bp = reshape(bp,length(bp),1);
wf = reshape(wf,length(wf),1);

end
