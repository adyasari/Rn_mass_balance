function F = myfun(pwr_lambda)
global xmin;
global alpha;
global tmean;
% H1,H2 replacement for igamma 
id='MATLAB:integral:MaxIntervalCountReached';
warning('off',id);
H1 = @(t) t.^(2-alpha-1).*exp(-t);
H2 = @(t) t.^(1-alpha-1).*exp(-t);
F(1)=integral(H1,0,Inf)*(pwr_lambda*integral(H2,pwr_lambda*xmin,Inf))^(-1)-tmean;
end