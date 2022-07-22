function[RN_out]=pwr(th,cgw)
%th=5400;
global tmean
tmean=th/3600;
global distr;
%display(tmean);
%cgw=12000;
global xmin;
global inc_gamm;
xmin=0.001;
global alpha 
alpha=1;
global gamma_alpha;
%gamma_alpha=2;
lambda=0.18/24;
gammaa=cgw*lambda;


if distr == 2 %power law
  %find parameter llambda2 for power law with cutoff with fr a defined tmean
  r=0.0747;
  %r=1;
  options = optimoptions('fsolve','Display','off');
  llambda2=fsolve(@myfun,r,options);
  
  %Convolution Integral to find output Rn concentration leaving the HZ
  H1 = @(t) t.^(1-alpha-1).*exp(-t);
  PPWR = @(t) t.^(-alpha).*exp(-llambda2.*t)*llambda2^(1-alpha)./integral(H1,llambda2*xmin,Inf);
  REX = @(t) exp(-lambda.*t).*PPWR(t);
  
  RN_out(1)=gammaa/lambda.*integral(PPWR,xmin,Inf)-gammaa/lambda.*integral(REX,xmin,Inf);
  RN_out(2)=integral(REX,xmin,Inf);

else %==3 gamma 
    
    gamma_beta=tmean/gamma_alpha;
    GAM = @(t) t.^(gamma_alpha-1).*(gamma_beta.^gamma_alpha.*gamma(gamma_alpha))^-1.*exp(-t./gamma_beta);
    REX = @(t) exp(-lambda.*t).*GAM(t);
        
    RN_out(1)=gammaa./lambda.*integral(GAM,0,Inf)-gammaa./lambda.*integral(REX,0,Inf);
    RN_out(2)=integral(REX,0,Inf);
    
end  
  
  
end