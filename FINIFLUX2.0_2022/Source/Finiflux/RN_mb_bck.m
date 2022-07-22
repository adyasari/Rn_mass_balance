observation_points = importdata('n_observation.dat'); %[-]
n_reaches = observation_points-1;
X=importdata('observation_points.txt');%[m]
width=importdata('river_width.txt');%[m]
Q=importdata('discharge.dat');%[m³/s]
depth_hz=importdata('depth_hz.dat');%[m]
por=importdata('porosity_hz.dat');%[-]
rt_hz=importdata('rtimes_hz.dat');%[sec]
Cin=importdata('nitrate_cin.dat');%[g/m^3]
distr=importdata('distribution.dat'); %1=exponential, 2=powerlaw, 3=gamma
%Monod Parameters
my_max=0.682; %[g/m³/h]
K=3.224; %[g/m³]
%first order
lambda = 0.04; %[1/h]
a=0.5; %lag time between tmean and ETmean for exponential model// Cutoff for GAmma and Power
global tmm;
global alpha;
global xmin;

if distr==3
    global gamma_alpha; 
    gamma_alpha=importdata('gamma_alpha.dat');
else
end

tmean=rt_hz/3600;
if distr == 1    %Exponential Model
    for i=1:observation_points-1;
    dx(i)=X(i+1)-X(i);%[m]
    average_flowriver(i)=Q(i)*3600; %[m³/h]
    Qh(i)=width(i)*depth_hz(i)*por(i)*((rt_hz(i)/3600)^-1)*dx(i); %[m³/h]
    Vh(i)=depth_hz(i)*width(i)*por(i)*dx(i); %[m³]
    
     pdf=@(t) ((tmean(i)-a)^-1)*exp(-t./(tmean(i)-a));
     %Ar(i)=integral(pdf,a,inf);
     %ET=@(t) (1/Ar(i)).*pdf(t+a);
     %ET=@(t) (1/Ar(i)).*pdf(t);
     %NO3 =@(t) -Qh(i)*Cin(i)*(exp(-lambda.*(t))-1).*(lambda)^-1;
     NO3_decay =@(t) exp(-lambda.*t);
     
     REX = @(t) Cin(i).*NO3_decay(t).*pdf(t);
          
     conv=integral(REX,0,Inf);   
     Cout(i)=conv;
%     
    %Cout(i)=(Qh(i)*Cin(i)*(rt_hz(i)/3600)*(lambda*exp(-a/(rt_hz(i)/3600))*(rt_hz(i)/3600)+1)^-1)/Vh(i); %[g/m³]
    Cout_monod(i)=Cin(i)/2-K/2+0.5*(Cin(i)^2+2*Cin(i)*K-2*Cin(i)*my_max*(Vh(i)/Qh(i))+K^2+2*K*my_max*(Vh(i)/Qh(i))+my_max^2*(Vh(i)/Qh(i))^2)^(1/2)-my_max*(Vh(i)/Qh(i))/2;
    dec_mass(i)=(Cin(i)-Cout(i))*Qh(i); %[g/h]
    dec_mass_monod(i)=(Cin(i)-Cout_monod(i))*Qh(i); %[g/h]
    dec_mass_rate(i)=dec_mass(i)/dx(i); %[g/h/m]
    dec_mass_rate_monod(i)=dec_mass_monod(i)/dx(i);%[g/h/m]
    river_mass(i)=average_flowriver(i)*((Cin(i+1)+Cin(i))/2); %[g/h/m]
    end
    
elseif distr == 2 %Power Law 
    tmean=rt_hz/3600;
    alpha =1;
    xmin=0.001;
    
            
    for i=1:observation_points-1;
        dx(i)=X(i+1)-X(i);%[m]
        average_flowriver(i)=Q(i)*3600; %[m³/h]
        Vh(i)=depth_hz(i)*width(i)*por(i)*dx(i); %[m³]
        Qh(i)=width(i)*depth_hz(i)*por(i)*(tmean(i)^-1)*dx(i); %[m³/h]
        tmm=tmean(i)-a;
        r=0.0747;
        options = optimoptions('fsolve','Display','off');
        llambda(i)=fsolve(@myfun_nit,r,options);
                
        NO3_decay =@(t) exp(-lambda.*t);
        %NO3 =@(t) -Qh(i)*Cin(i)*(exp(-lambda.*t)-1).*(lambda)^-1; 
        H1 = @(t) t.^(1-alpha-1).*exp(-t);
        PPWR = @(t) t.^(-alpha).*exp(-llambda(i).*t)*llambda(i)^(1-alpha)./integral(H1,llambda(i)*xmin,Inf);
         
        REX = @(t) Cin(i).*NO3_decay(t).*PPWR(t);        
        conv=integral(REX,xmin,Inf);      
        Cout(i)=conv;
        dec_mass(i)=(Cin(i)-Cout(i))*Qh(i); %[g/h]
        dec_mass_rate(i)=dec_mass(i)/dx(i); %[g/h/m]
        river_mass(i)=average_flowriver(i)*((Cin(i+1)+Cin(i))/2); %[g/h/m]
        
    end
else %Gamma Model
    
    tmean=rt_hz/3600;
    gamma_beta=(tmean)./gamma_alpha;
    
    for i=1:observation_points-1;
        dx(i)=X(i+1)-X(i);%[m]
        average_flowriver(i)=Q(i)*3600; %[m³/h]
        Vh(i)=depth_hz(i)*width(i)*por(i)*dx(i); %[m³]
        Qh(i)=width(i)*depth_hz(i)*por(i)*(tmean(i)^-1)*dx(i); %[m³/h]
        
        gamma_zeta=(tmean-a)./gamma_alpha;       
        NO3_decay =@(t) exp(-lambda.*t);
        %NO3 =@(t) -Qh(i)*Cin(i)*(exp(-lambda.*t)-1).*(lambda)^-1;       
        GAM = @(t) t.^(gamma_alpha-1).*(gamma_zeta(i)^gamma_alpha.*gamma(gamma_alpha))^-1.*exp(-t./gamma_zeta(i));     
        %Ar(i)=integral(GAM,a,inf);
        %ET=@(t) (1/Ar(i)).*GAM(t+a);
        %ET=@(t) (1/Ar(i)).*pdf(t);
                
        REX = @(t) Cin(i).*NO3_decay(t).*GAM(t);  
        conv=integral(REX,0,Inf);
        Cout(i)=conv;
        dec_mass(i)=(Cin(i)-Cout(i))*Qh(i); %[g/h]
        dec_mass_rate(i)=dec_mass(i)/dx(i); %[g/h/m]
        river_mass(i)=average_flowriver(i)*((Cin(i+1)+Cin(i))/2); %[g/h/m]
    end
end
       

figure(1)
ax1=subplot(2,1,1);


for i=1:X(observation_points)
    XX(i,1)=i/1000;
    for j=1:observation_points-1
        if i >= X(j)&& i<= X(j+1)
        XX(i,2)=dec_mass_rate(j)/1000;
        break;
        else
        end
    end
end

for i=1:X(observation_points)  
    if i==1
        XX(i,3)=0;
    else
        XX(i,3)=XX(i-1,3)+XX(i,2);
    end 
end

plot(XX(1:end,1),XX(1:end,2));
ylabel(ax1,'MASS DECAY [kg*h^-1*m^-1]');
xlabel(ax1,'X [km]');

ax2=subplot(2,1,2);
plot(XX(1:end,1),XX(1:end,3));
ylabel(ax2,'CUMM MASS DECAY [kg*h^-1]');
xlabel(ax2,'X [km]');

figure(2)
if distr == 1    
    boxplot(tmean,'notch','on')
    xlabel('-');
    ylabel('Exposure Times [h]');
    title('PDF Exposure Times - EXPONENTIAL FIT');
    
    
    for i=1:X(observation_points)
    YY(i,1)=i/1000;
    for j=1:observation_points-1
        if i >= X(j)&& i<= X(j+1)
        YY(i,2)=dec_mass_rate_monod(j)/1000;
        break;
        else
        end
    end
end

for i=1:X(observation_points)  
    if i==1
        YY(i,3)=0;
    else
        YY(i,3)=YY(i-1,3)+YY(i,2);
    end 
end

figure (3)
ax1=subplot(2,1,1);
    
plot(YY(1:end,1),YY(1:end,2));
ylabel(ax1,'MASS DECAY Monod [kg*h^-1*m^-1]');
xlabel(ax1,'X [km]');

ax2=subplot(2,1,2);
plot(YY(1:end,1),YY(1:end,3));
ylabel(ax2,'CUMM MASS DECAY Monod [kg*h^-1]');
xlabel(ax2,'X [km]');    
     
    
elseif distr == 2  
    boxplot(tmean,'notch','on')
    ylabel('RT [h]');
    xlabel('-');
    title('PDF Residence Times - TRUNCATED POWER LAW');
else
    boxplot(tmean,'notch','on')
    ylabel('RT [h]');
    xlabel('-');
    title('PDF Residence Times - GAMMA');   
end

figure(4)
xlabel('RT [h]');
ylabel('prop [1/h]');    
   
if distr == 1  
   title('PDF Residence Times - EXPONENTIAL');        
    for i=1:observation_points-1;
    hold on
    dist=@(t) 1/tmean(i).*exp(-t./tmean(i));
    fplot(dist,[0 10]);
    if a~=0
        line([a a], get(gca, 'ylim'));
    else
    end
    hold off
    end
elseif distr == 2
    title('PDF Residence Times - TRUNCATED POWER LAW');        
    for i=1:observation_points-1;
    hold on
    dist=@(t) t.^(-alpha).*exp(-llambda(i).*t)*llambda(i)^(1-alpha)./integral(H1,llambda(i)*xmin,Inf);
    fplot(dist,[0 10]);
    if a~=0
       line([a a], get(gca, 'ylim'));
    else
    end
    hold off
    end
else 
    title('PDF Residence Times - GAMMA');        
    for i=1:observation_points-1;
    hold on
    dist=@(t) t.^(gamma_alpha-1).*(gamma_beta(i)^gamma_alpha.*gamma(gamma_alpha))^-1.*exp(-t./gamma_beta(i));
    fplot(dist,[0 10]);
    if a~=0
       line([a a], get(gca, 'ylim'));
    else
    end
    
    end
end
   
