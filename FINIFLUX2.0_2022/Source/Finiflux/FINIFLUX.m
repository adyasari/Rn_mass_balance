%---------------------MODEL DESCRIPTION---------------------------------
%FINIFLUX is an implicit Finite Element model that numerically solves the
%steady state mass balance equation for Radon (Rn) in lotic systems.
%Degasing is represented using two very popular models which originally are
%from O'Connor and Dobbins (1958) and Negulescu and Rojanski (1969) with
%modifications from Cartwright et al. (2011). Alternativley, FINIFLUX also can work
%with degasing fluxes specified by the user. FINIFLUX is intended to
%estimate groundwater fluxes into river systems as well as hyporheic
%exchnage characteristics based on measured Rn concentrations. The model is
%coupled to the optimasation software package PEST (Doherty 2010) for
%inverse parameter estimation. FINIFLUX is intended to help scientists and
%authorities that use the Rn tecnique to estimate surface-groundwater
%exchange for river systems at the reach scale.
%Authors: S. Frei (sven.frei@uni-bayreuth.de)
%         B.S. Gilfedder (benjamin-silas.gilfedder@uni-bayreuth.de)
%June 2015
%June 2016: We included the posssiblity to account for surficial Rn inflow
%via tributaries or drainage.
%Mai 2017: We added uppwinding Finite Element correction to reduce the effect of numerical dispersion.
% Hyporheic areas can now be represented using 1) exponential, 2) power law
% and 3) gamma distribted RT. The entire code can be compiled using
% MATLAb's application compiler. 
%February 2022: fixed a bug in the degassing model 2
%------------DATA INPUT SECTION-------------------------------------------- 
%number of observation points
global observation_points;
observation_points = importdata('n_observation.dat');

%number of reaches
global n_reaches;
n_reaches = observation_points-1;

%radioactive decay constant for radon [s^-1]
global lambda;
lambda = 0.18/86400;

%X-location [m] of observation points
X=importdata('observation_points.txt');

%river width [m] for each reach; 
width=importdata('river_width.txt'); 

%river depth [m] for each reach; 
depth=importdata('river_depth.txt'); 

%GW inflow rate [m^3/(m*s)];
I=importdata('gw_inflow.dat');

%Rn endmember concentration of GW [Bq/m�];
cgw=importdata('gw_concentration.txt');

%river discharge [m�/s]
Q=importdata('discharge.dat');

%Rn concentration at x=0 [Bq/m�]
initial_C = importdata('initial_c.dat');

%switch with/without hyporheic exchange 1 = with hz; 0 = without hz
sw_hz=importdata('hz_on_off.dat');

%switch degasing model 1= O'Connor and Dobbins (1958) model, 2=Negulescu and Rojanski (1969), 3=user defined degasing values 
deg_k=importdata('k_degas.dat');

%switch for hyporheic one RT distribution model, 1=exponential, 2=powerlaw, 3=gamma
global distr;
distr=importdata('distribution.dat'); 

%inflow from tributaries [m�/s], zero for no river inflow reaches
Qin=importdata('river_inflow.dat'); 

%inflow_length for tributaries [m], zero for no inflow reaches
inflow_length=importdata('inflow_length.dat'); 

%Rn concentration for inflowing tributaries [Bq/m�], zero for no inflow reaches
cin=importdata('inflow_conc.dat'); 

%output of mass balance on/off [0] off [1] on, For PEST optimizzation turn
%off mass blanace
switch_mass=0; 

%Petrov Galerkin Finite Element for upstream weighting; off beta = 0; on%beta >=0
beta = importdata('upstream_weight.dat');
alp=0.75;

%for gamma distribution model define alpha parameter
if distr==3
    global gamma_alpha; 
    gamma_alpha=importdata('gamma_alpha.dat');
else
end



%------------MAIN ROUTINE--------------------------------------------------

for i=1:numel(inflow_length)
    if inflow_length(i)==0
        Rl(i)=1;
    else
        Rl(i)=inflow_length(i);
    end
end
        
j=1;
for i=1:n_reaches
    if Qin(i)==0
        flag(j)=1;
        j=j+1;
    else
        flag(j)=1;
        j=j+1;
    end
end




%Element discretization [m]; dx
for i=1:observation_points-1;
dx(i)=X(i+1)-X(i);
end
dx=dx.';

%degasing
%cross-section and velocity
A=depth.*width;
v=(Q./A)*86400;
%degasing k1 [m/s]

if deg_k == 1    
    k=((9.301*10^-3*((v.^0.5)./(depth.^1.5))).*depth)/86400;
elseif deg_k == 2
    k=((4.87*10^-4*(v./depth).^0.85)).*depth/86400;
elseif deg_k == 3
    k=importdata('k_values.dat'); 
end

    
if sw_hz == 1 

    %depth of hyporheic zone (m)
    h=importdata('depth_hz.dat');
    %porosity hyporheic zone (-)
    porosity=importdata('porosity_hz.dat');
    %residence time hz (s)
    th=importdata('rtimes_hz.dat');
    %gamma hz (Bq/m�/s)
    %gamma=importdata('gamma_hz.dat');
    gamma=cgw.*lambda;
    
    %hyporheic exchange coeff. 
    if distr==1 %exponential
        alpha1=(gamma.*h.*width.*porosity)./(1+lambda.*th);
        alpha2=(lambda.*h.*width.*porosity)./(1+lambda.*th);
    else
        for i=1:n_reaches %power law with cutoff, gamma
        RN_HZ{i}=pwr(th(i),cgw(i));
        alpha1(i)=RN_HZ{i}(1)*h(i)*width(i)*porosity(i)/th(i);
        alpha2(i)=-(RN_HZ{i}(2)*h(i)*width(i)*porosity(i)/th(i)-h(i)*width(i)*porosity(i)/th(i));
        end 
    end
    
    entry0_0=[];
    entry0_1=[];
    entry1_0=[];
    entry1_1=[];
    vector0=[];
    vector1=[];
    upstream_corr0_0=[];
    upstream_corr0_1=[];
    upstream_corr1_1=[];
    upstream_corr1_0=[];
    corr_vector0=[];
    corr_vector1=[];

        for i=1:n_reaches
            entry0_0(i)=1/3*(-I(i)*flag(i)-alpha2(i)*flag(i)-k(i)*width(i)*flag(i)-depth(i)*width(i)*lambda*flag(i)-Qin(i)/Rl(i))*X(i)+1/3*(I(i)*flag(i)+alpha2(i)*flag(i)+k(i)*width(i)*flag(i)+depth(i)*width(i)*lambda*flag(i)+Qin(i)/Rl(i))*X(i+1)-Q(i)*flag(i)/2;
             entry0_1(i)=1/6*(-I(i)*flag(i)-alpha2(i)*flag(i)-k(i)*width(i)*flag(i)-depth(i)*width(i)*lambda*flag(i)-Qin(i)/Rl(i))*X(i)+1/6*(I(i)*flag(i)+alpha2(i)*flag(i)+k(i)*width(i)*flag(i)+depth(i)*width(i)*lambda*flag(i)+Qin(i)/Rl(i))*X(i+1)+Q(i)*flag(i)/2;
             entry1_0(i)=1/6*(-I(i)*flag(i)-alpha2(i)*flag(i)-k(i)*width(i)*flag(i)-depth(i)*width(i)*lambda*flag(i)-Qin(i)/Rl(i))*X(i)+1/6*(I(i)*flag(i)+alpha2(i)*flag(i)+k(i)*width(i)*flag(i)+depth(i)*width(i)*lambda*flag(i)+Qin(i)/Rl(i))*X(i+1)-Q(i)*flag(i)/2;
             entry1_1(i)=1/3*(-I(i)*flag(i)-alpha2(i)*flag(i)-k(i)*width(i)*flag(i)-depth(i)*width(i)*lambda*flag(i)-Qin(i)/Rl(i))*X(i)+1/3*(I(i)*flag(i)+alpha2(i)*flag(i)+k(i)*width(i)*flag(i)+depth(i)*width(i)*lambda*flag(i)+Qin(i)/Rl(i))*X(i+1)+Q(i)*flag(i)/2;
             upstream_corr0_0(i)=1/3*beta*alp*(2*Q(i)-(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-alpha2(i)-Qin(i)/Rl(i))*X(i)+(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-alpha2(i)-Qin(i)/Rl(i))*X(i+1));
             upstream_corr0_1(i)=-1/3*beta*alp*(2*Q(i)+(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-alpha2(i)-Qin(i)/Rl(i))*X(i)-(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-alpha2(i)-Qin(i)/Rl(i))*X(i+1));
             upstream_corr1_0(i)=-1/3*beta*alp*(2*Q(i)-(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-alpha2(i)-Qin(i)/Rl(i))*X(i)+(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-alpha2(i)-Qin(i)/Rl(i))*X(i+1));
             upstream_corr1_1(i)=1/3*beta*alp*(2*Q(i)+(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-alpha2(i)-Qin(i)/Rl(i))*X(i)-(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-alpha2(i)-Qin(i)/Rl(i))*X(i+1));
             vector0(i)=-0.5*((I(i)*cgw(i)*flag(i)+alpha1(i)*flag(i)+(Qin(i)/Rl(i))*cin(i))*(X(i)-X(i+1)));
             vector1(i)=-0.5*((I(i)*cgw(i)*flag(i)+alpha1(i)*flag(i)+(Qin(i)/Rl(i))*cin(i))*(X(i)-X(i+1)));
             corr_vector0(i)=1/3*beta*2*(I(i)*cgw(i)+alpha1(i)+Qin(i)*cin(i)/Rl(i))*(X(i)-X(i+1)); %check this term (I*cgw+alpha1*dx)
             corr_vector1(i)=-1/3*beta*2*(I(i)*cgw(i)+alpha1(i)+Qin(i)*cin(i)/Rl(i))*(X(i)-X(i+1)); %check this term
             AM(:,:,i) = zeros(observation_points); 
             CO(:,:,i) = zeros(observation_points);
            R(observation_points,1,i) = 0;
            R_cor(observation_points,1,i) = 0;
        end 

        for i=1:n_reaches   
        AM(i,i,i) = entry0_0(i);
        AM(i,i+1,i) = entry0_1(i);
        AM(i+1,i,i) = entry1_0(i);
        AM(i+1,i+1,i) = entry1_1(i);
        CO(i,i,i) = upstream_corr0_0(i);
        CO(i,i+1,i) = upstream_corr0_1(i);
        CO(i+1,i,i) = upstream_corr1_0(i);
        CO(i+1,i+1,i) = upstream_corr1_1(i);
        R(i,1,i) = vector0(i);
        R(i+1,1,i) = vector1(i);
        R_cor(i,1,i) = corr_vector0(i);
        R_cor(i+1,1,i) = corr_vector1(i);
        end

        %global Matrix and Vector assembly
        Global= zeros(observation_points);
        R_global (observation_points,1)= 0;

        for i=1:n_reaches
            Global=Global+AM(:,:,i)+CO(:,:,i);
            R_global=R_global+R(:,:,i)+R_cor(:,:,i);
        end

        %modify_matrix and vector for IC
        Global(1,1)=Global(1,1)*10^6;
        R_global(1)=Global(1,1)*initial_C;
        Global=sparse(Global);

        %solve system
        Cx=Global\R_global;

        %scatter(X,Cx);
        %plot(X,Cx,'-o');
        II=[X,Cx];
        dlmwrite('Rn_modeled.csv',II,'delimiter','\t','precision', 16)

else 
    
    for i=1:n_reaches 
        alpha1(i)=0;
        alpha2(i)=0;
        porosity(i)=0;
        h(i)=0;
        th(i)=0;
    end 

    entry0_0=[];
    entry0_1=[];
    entry1_0=[];
    entry1_1=[];
    vector0=[];
    vector1=[];
    upstream_corr0_0=[];
    upstream_corr0_1=[];
    upstream_corr1_1=[];
    upstream_corr1_0=[];
    corr_vector0=[];
    corr_vector1=[];

        for i=1:n_reaches
            entry0_0(i)=1/3*(-I(i)*flag(i)-k(i)*width(i)*flag(i)-depth(i)*width(i)*lambda*flag(i)-Qin(i)/Rl(i))*X(i)+1/3*(I(i)*flag(i)+k(i)*width(i)*flag(i)+depth(i)*width(i)*lambda*flag(i)+Qin(i)/Rl(i))*X(i+1)-Q(i)*flag(i)/2;
             entry0_1(i)=1/6*(-I(i)*flag(i)-k(i)*width(i)*flag(i)-depth(i)*width(i)*lambda*flag(i)-Qin(i)/Rl(i))*X(i)+1/6*(I(i)*flag(i)+k(i)*width(i)*flag(i)+depth(i)*width(i)*lambda*flag(i)+Qin(i)/Rl(i))*X(i+1)+Q(i)*flag(i)/2;
             entry1_0(i)=1/6*(-I(i)*flag(i)-k(i)*width(i)*flag(i)-depth(i)*width(i)*lambda*flag(i)-Qin(i)/Rl(i))*X(i)+1/6*(I(i)*flag(i)+k(i)*width(i)*flag(i)+depth(i)*width(i)*lambda*flag(i)+Qin(i)/Rl(i))*X(i+1)-Q(i)*flag(i)/2;
             entry1_1(i)=1/3*(-I(i)*flag(i)-k(i)*width(i)*flag(i)-depth(i)*width(i)*lambda*flag(i)-Qin(i)/Rl(i))*X(i)+1/3*(I(i)*flag(i)+k(i)*width(i)*flag(i)+depth(i)*width(i)*lambda*flag(i)+Qin(i)/Rl(i))*X(i+1)+Q(i)*flag(i)/2;
             upstream_corr0_0(i)=1/3*beta*alp*(2*Q(i)-(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-Qin(i)/Rl(i))*X(i)+(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-Qin(i)/Rl(i))*X(i+1));
             upstream_corr0_1(i)=-1/3*beta*alp*(2*Q(i)+(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-Qin(i)/Rl(i))*X(i)-(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-Qin(i)/Rl(i))*X(i+1));
             upstream_corr1_0(i)=-1/3*beta*alp*(2*Q(i)-(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-Qin(i)/Rl(i))*X(i)+(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-Qin(i)/Rl(i))*X(i+1));
             upstream_corr1_1(i)=1/3*beta*alp*(2*Q(i)+(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-Qin(i)/Rl(i))*X(i)-(-I(i)-k(i)*width(i)-depth(i)*width(i)*lambda-Qin(i)/Rl(i))*X(i+1));            
             vector0(i)=-0.5*((I(i)*cgw(i)*flag(i)+(Qin(i)/Rl(i))*cin(i))*(X(i)-X(i+1)));
             vector1(i)=-0.5*((I(i)*cgw(i)*flag(i)+(Qin(i)/Rl(i))*cin(i))*(X(i)-X(i+1)));
             corr_vector0(i)=1/3*beta*2*(I(i)*cgw(i)+Qin(i)*cin(i)/Rl(i))*(X(i)-X(i+1)); %check this term (I*cgw+alpha1*dx)
             corr_vector1(i)=-1/3*beta*2*(I(i)*cgw(i)+Qin(i)*cin(i)/Rl(i))*(X(i)-X(i+1)); %check this term         
             AM(:,:,i) = zeros(observation_points); 
             CO(:,:,i) = zeros(observation_points);
            R(observation_points,1,i) = 0;
             R_cor(observation_points,1,i) = 0;
        end 

        for i=1:n_reaches   
        AM(i,i,i) = entry0_0(i);
        AM(i,i+1,i) = entry0_1(i);
        AM(i+1,i,i) = entry1_0(i);
        AM(i+1,i+1,i) = entry1_1(i);
        CO(i,i,i) = upstream_corr0_0(i);
        CO(i,i+1,i) = upstream_corr0_1(i);
        CO(i+1,i,i) = upstream_corr1_0(i);
        CO(i+1,i+1,i) = upstream_corr1_1(i);
        R(i,1,i) = vector0(i);
        R(i+1,1,i) = vector1(i);
        R_cor(i,1,i) = corr_vector0(i);
        R_cor(i+1,1,i) = corr_vector1(i);
        end

        %global Matrix and Vector assembly
        Global= zeros(observation_points);
        R_global (observation_points,1)= 0;

        for i=1:n_reaches
            Global=Global+AM(:,:,i)+CO(:,:,i);
            R_global=R_global+R(:,:,i)+R_cor(:,:,i);
        end

        %modify_matrix and vector for IC
        Global(1,1)=Global(1,1)*10^6;
        R_global(1)=Global(1,1)*initial_C;
        Global=sparse(Global);

        %solve system
        Cx=Global\R_global;

        
        %scatter(X,Cx);
        %plot(X,Cx,'-o');
        II=[X,Cx];
        dlmwrite('Rn_modeled.csv',II,'delimiter','\t','precision', 16)
       
  
end


if switch_mass == 1
%massbalance for all reaches-----
%fluxes
for i=1:n_reaches
    %Rn Hyporheic Flux [Bq/s]
    H(i)=alpha1(i)*dx(i)-alpha2(i)*(Cx(i)+Cx(i+1))*dx(i)/2;
    H1(i)=-(porosity(i)*h(i)*width(i))/(lambda*th(i)+1)*dx(i)*(lambda/2*(Cx(i)+Cx(i+1))-gamma(i));
    %Rn Groundwater Flux [Bq/s]
    GW(i)=I(i)*dx(i)*cgw(i)-I(i)*(Cx(i)+Cx(i+1)*dx(i)/2);
    %Rn Advective Flux [Bq/s]
    Ad(i)=Q(i)*(Cx(i+1)-Cx(i));
    %Rn Degasing Flux [Bq/s]
    Deg(i)=-k(i)*width(i)*(dx(i)/2)*(Cx(i)+Cx(i+1));
    %Rn Decay Flux [Bq/s]
    Dec(i)=-depth(i)*width(i)*lambda*dx(i)/2*(Cx(i)+Cx(i+1));
    %Rn Inflow/Outflow [Bq/s] maybe wrong
    IO(i)=Qin(i)*cin(i)-Qin(i)*(Cx(i+1)+Cx(i))/2;
    
end

%mass balance
for i=1:n_reaches
    %Rn Hyporheic Flux [Bq/s] always positive
    H(i)=alpha1(i)*dx(i)-alpha2(i)*(Cx(i)+Cx(i+1))*dx(i)/2;
    H1(i)=-(porosity(i)*h(i)*width(i))/(lambda*th(i)+1)*dx(i)*(lambda/2*(Cx(i)+Cx(i+1))-gamma(i));
    %Rn Groundwater Flux [Bq/s] always positive
    GW(i)=I(i)*dx(i)*cgw(i)-I(i)*(Cx(i)+Cx(i+1)*dx(i)/2);
    %Rn Advective Flux [Bq/s] positive or negative
    Ad(i)=-Q(i)*(Cx(i+1)-Cx(i));
    %Rn Degasing Flux [Bq/s] always negative
    Deg(i)=-k(i)*width(i)*(dx(i)/2)*(Cx(i)+Cx(i+1));
    %Rn Decay Flux [Bq/s] always negative
    Dec(i)=-depth(i)*width(i)*lambda*dx(i)/2*(Cx(i)+Cx(i+1));
    %Rn Inflow/Outflow [Bq/s] positive or negative
    IO(i)=Qin(i)*cin(i)-Qin(i)*(Cx(i+1)+Cx(i))/2;
    
    %Balance for reach
    In(i)=H(i)+GW(i);
    Out(i)=Deg(i)+Dec(i);
    if Ad(i) > 0
        In(i)=In(i)+Ad(i);
    else
        Out(i)=Out(i)+Ad(i);
    end
    if IO(i) > 0
        In(i)=In(i)+IO(i);
    else
        Out(i)=Out(i)+IO(i);
    end
    
    Res(i)=In(i)-Out(i);
    Res_rel(i)=Res(i)/In(i);
    
end


total_Res=Q(1)*initial_C+sum(H)+sum(GW)+sum(Deg)+sum(Dec)+sum(IO)-Q(n_reaches)*Cx(n_reaches); 
total_Res_rel=total_Res/(Q(1)*initial_C+sum(H)+sum(GW));

text1 = ['Residual Rn Mass Balance [Bq/s]: ',num2str(total_Res)];
text2 = ['Relative Error (Residual/Inflow) [%]: ',num2str(total_Res_rel*100)];
disp(text1);
disp(text2);

figure (2);
 ax1=subplot(4,1,1);
 plot(H,'r-o');
 ylabel(ax1,'Rn Hyporheic Flux [Bq/s]');
 %xlabel(ax1,'Reach Number [-]');
 
 ax2=subplot(4,1,2);
 plot(GW,'r-o');
 ylabel(ax2,'Rn Groundwater Flux [Bq/s]');
 %xlabel(ax2,'Reach Number [-]');
 
  ax3=subplot(4,1,3);
 plot(Ad,'r-o');
 ylabel(ax3,'Rn Advective Flux [Bq/s]');
 %xlabel(ax3,'Reach Number [-]');
 
  ax4=subplot(4,1,4);
 plot(Deg,'r-o');
 ylabel(ax4,'Rn Degasing Flux [Bq/s]');
 xlabel(ax4,'Reach Number [-]');
 
 figure(3);
  ax5=subplot(4,1,1);
 plot(Dec,'r-o');
 ylabel(ax5,'Rn Decay Flux [Bq/s]');
 %xlabel(ax5,'Reach Number [-]');
 
 ax6=subplot(4,1,2);
 plot(IO,'r-o');
 ylabel(ax6,'Rn Inflow/Outflow Flux [Bq/s]');
 %xlabel(ax5,'Reach Number [-]');
 
 ax7=subplot(4,1,3);
 plot(Res_rel*100,'r-o');
 ylabel(ax7,'Residual Error [%]');
 %xlabel(ax7,'Reach Number [-]');
 
  ax8=subplot(4,1,4);
 plot(dx,'r-o');
 ylabel(ax8,'Reach Length [m]');
 xlabel(ax8,'Reach Number [-]');
else
    exit
end
 
 
 



