clear all;
clc;
%This code is run with example data
%First load multiple experiments

%In this case: 4 experiments at 4 different temperatures
%I: number of experiments (4)



I = [1,2,3,4];

%Generates an overall input matrix Dall by combining all input matrices D(spectra over time)
%input_1: name of first experiment
%change i for the numbers of experiments (4)

for i=1:4
    load(['input_',num2str(I(i))]);
    Dall{i} = D(1:end,:);
    len(i) = size(D,1);
end

%mMrges all timevectors to an overall time vector Tall
%t1: time vector of first experiment in seconds
%tall: overall time vector

tall{1} = t1(2:end)';
tall{2} = t2(2:end)';
tall{3} = t3(2:end)';
tall{4} = t4(2:end)';

%Enter temperatures for each experiment in  °C
%T: temperature vector
T = [15,25,35,45];

% Do not change the next two lines
s = 2;
z = 2;

%Change j for the numbers of experiments (6)

for j=1:4
    [U{j},S{j},V{j}] = svds(Dall{j},z);
    for i=1:z
        if min(V{j}(:,i)) < -max(V{j}(:,i))
            V{j}(:,i) = -V{j}(:,i);
            U{j}(:,i) = -U{j}(:,i);
        end
    end
    L{j} = U{j}*S{j};
    pL{j} = pinv(L{j});
    R{j} = V{j}';
end

%Set constraints
% 1 for applying
% 0 for applying
%First number for kinetic model as constraint
%Second number for only positive concentration profiles
%Third number for only positive values in spectrum

w = [1;1;0];

%Enter initial educt concentrations for each experiment in mol/L

c0 = [1 1 0 0; 1 1 0 0; 1 1 0 0;1 1 0 0];

%If required enter initial catalyst concentrations for each experiment in the following line in mol/L
%cat= [0.005;0.005;0.005;0.0025;0.0025;0.0025];


p = 1;
%Initial estimation of reaction rate coefficient K0(1) and activation energy K0(2)
K0 = [0.2; 80];

%Nonlinear least-squares solver (matlab function lsqnonlin)
%If required enter cat in optimGES(K,tall,T,L,pL,R,w,c0,cat,len)
[Kopt,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(@(K) optimGES(K,tall,T,L,pL,R,w,c0,len), K0, zeros(length(len)+2,p),[],optimset('Display','Iter'));

%calculation of Confidence Intervalls (CI)
conf = nlparci(Kopt,residual,'jacobian',jacobian);
%calculation of Confidence Intervalls (CI) in percent
conf_pro=(Kopt-conf(:,1))./Kopt*100

%calculation of Covariance Matrix vK0
vK0 = resnorm*inv(jacobian'*jacobian)/length(K0);
%calculation of Correlation Matrix  corr_matr
corr_matr = corrcov(vK0); 

% K is final result (reaction rate coefficient K(1) and activation energy K(2))
K = Kopt



for i=1:length(len)
%If required enter cat in rhs(t,c,K,cat(i),T(:,i)),[0; tall{i}],c0(i,1:4));
[t,Code{i}] = ode15s(@(t,c) rhs(t,c,K,T(:,i)),[0; tall{i}],c0(i,1:4));
 
 ct{i} = pL{i}*Code{i}(1:end,[1 3]); %four components are in the initial concentration, but only two are can be modelled using MCR [1 3]
 Csvd{i} = L{i}*ct{i};  
 Ssvd{i} = pinv(ct{i})*R{i};
      
    
    figure(1)
    hold on;
    plot(t,Code{i}(:,[1 3]),'k-.'); %four components are in the initial concentration, but only two are can be modelled using MCR [1 3]
    scatter([0; tall{i}],Csvd{i},'linewidth',2);
    hold off;
    
end





for i=1:length(len)
E{i}=Dall{i}-Csvd{i}*Ssvd{i};           %Calculation of Residual Matrices
u{i}=sum(sum(E{i}.*E{i}));      
v{i}=sum(sum(Dall{i}.*Dall{i}));        %Summe der Einträge der Inputmatrix
R2{i}=100*sqrt((v{i}-u{i})/v{i});       %Calculation of explained data variance
LoF{i}=100*sqrt(u{i}/v{i});             %Calculation ofLack of Fit
end

save Result

%If required enter cat in optimGES(K,tall,T,L,pL,R,w,c0,cat,len)
function R = optimGES(K,tall,T,L,pL,R,w,c0,len)
i0 = 0;
Code = [];
R1 = [];
R2 = [];
R3 = [];
for i=1:length(len)
    %If required enter cat in rhs(t,c,K,cat(i),T(:,i))
    [t,code] = ode15s(@(t,c) rhs(t,c,K,T(:,i)),[0; tall{i}],c0(i,1:4));
    Code{i} = code(1:end,:);
    ct = pL{i}*Code{i}(:,[1 3]); %four components are in the initial concentration, but only two are can be modelled using MCR [1 3]
    Csvd = L{i}*ct;
    Ssvd = pinv(ct)*R{i};
    r1 = Code{i}(:,[1 3])-Csvd; %four components are in the initial concentration, but only two are can be modelled using MCR [1 3]
    r2 = min(Csvd,0);
    r3 = min(Ssvd,0);
    R1 = [R1; r1(:)];
    R2 = [R2; r2(:)];
    R3 = [R3; r3(:)];
end

R = [w(1)*R1; w(2)*R2; w(3)*R3]; 
end



%Selection of kinetic model, correlation with reparameterized Arrhenius equation
%in this example:
%kinetic follows second order:
%d(1)/dT=-k*(c(1)^1)*(c(2)^1)
%change kinetic model 
%if neccesarry you can also add cat in rhs(t,c,K,cat,T)
%concentration here!
%set reference temperature in Kelvin (here 298.15)
function dc = rhs(t,c,K,T)
dc = c;
dc(1) = -(K(1)*exp((-K(2)/0.008314)*(1./(T+273.15)-1./298.1500)))*(c(1)^1)*(c(2)^1);%*(cat^1);
dc(2) = -(K(1)*exp((-K(2)/0.008314)*(1./(T+273.15)-1./298.1500)))*(c(1)^1)*(c(2)^1);%*(cat^1);
dc(3) = (K(1)*exp((-K(2)/0.008314)*(1./(T+273.15)-1./298.1500)))*(c(1)^1)*(c(2)^1);%*(cat^1);
dc(4) = (K(1)*exp((-K(2)/0.008314)*(1./(T+273.15)-1./298.1500)))*(c(1)^1)*(c(2)^1);%*(cat^1);
end