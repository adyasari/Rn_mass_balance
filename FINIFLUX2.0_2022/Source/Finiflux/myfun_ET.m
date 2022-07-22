function F = myfun_ET(start_tau)
global cutoff_area;
global a;
% H1,H2 replacement for igamma 
id='MATLAB:integral:MaxIntervalCountReached';
warning('off',id);
H1 = @(t) 1/start_tau.*exp(-t./start_tau);
F(1)=integral(H1,a,Inf)-cutoff_area;
end