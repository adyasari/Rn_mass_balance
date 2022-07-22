alpha = 0.5;
beta = 1.75;
gamma_pdf = @(t) t.^(alpha-1).*(beta.^alpha.*gamma(alpha))^-1.*exp(-t./beta);
gamma_t = @(t) t.*gamma_pdf(t);

tmean=integral(gamma_t,0,inf);
area1=integral(gamma_pdf,0,inf);

global cutoff_area;
global a;
a=2.0;
cutoff_area=1-integral(gamma_pdf,a,inf);
Ar=integral(gamma_pdf,a,inf);
r=1;
options = optimoptions('fsolve','Display','off');
tau=fsolve(@myfun_ET,r,options);

cor_exp = @(t) 1/tau.*exp(-t./tau);
area2=integral(cor_exp,a,inf);

ET=@(t) 1/Ar*gamma_pdf(t+a);
ET2=@(t) gamma_pdf(t)+cor_exp(t);
area3=integral(ET,0,inf);

ET3=@(t) t.*ET(t);
ET_mean=integral(ET3,0,inf);

%pdf
figure (1)
hold on
fplot(gamma_pdf,[0 10])
 line([a a], get(gca, 'ylim'));
hold off


figure (2)
hold on
fplot(gamma_pdf,[a 10])
fplot(cor_exp,[a 10])
fplot(ET2,[a 10])
hold off



figure (3)
hold on
%fplot(gamma_pdf,[0 10])
%fplot(cor_exp,[0 10])
fplot(ET,[0 10])
hold off


%-------option 2

ET4=@(t) (Ar^-1)*gamma_pdf(t);
ET5=@(t) (Ar^-1)*gamma_pdf(t+a);
ET5mean=@(t) t.*ET5(t);
ET5_mean=integral(ET5mean,0,inf);
ET5_area=integral(ET5,0,inf);

figure (4)
hold on
%fplot(gamma_pdf,[0 10])
%fplot(cor_exp,[0 10])
fplot(gamma_pdf,[a 10],'b')
fplot(ET4,[a 10],'r')
hold off

figure (5)
hold on
%fplot(gamma_pdf,[0 10])
%fplot(cor_exp,[0 10])
%fplot(gamma_pdf,[a 10],'b')
fplot(ET5,[0 10],'r')
hold off





% xmin=0.01;
% alpha = 0.074697150906270;
% power = @(t) (alpha-1).*xmin^(alpha-1).*t.^(-alpha);
% power_t =@(t) t.*power(t);
% fplot(power,[0 10])
% integral(power,xmin,inf)
% integral(power_t,xmin,inf)