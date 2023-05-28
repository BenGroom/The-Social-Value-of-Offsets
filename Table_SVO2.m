%Make table on adjustment factors
clear; close all;
%constant parameters
    g=0.02;
    eta=1.35;
    delta=0.005;
    discount=delta+eta*g;
    GDP0=85000;
    T0=1.2; %°C
    zeta=0.0006;
    tau1=3;
%parameters which will vary
    Endtime=[25 50 100 395];
    tau2=Endtime+tau1;
    phivarphi=[0 0.005 0.01];
    varphi2=[1000 0.5 0.25];
    kappa=[0.0077 0.0025]; % output loss for 1°C
    gamma=2*kappa; %total damages=exp(-gamma/2^*T²)
%Load Temperature path and interpolate.
    [TEMP] = xlsread('iamc_db_Temp.xlsx',1,'H1:P8','basic');%starts in 2020
    TEMP=TEMP([3 4 7 8],:); %select scenario SSP1-26, SSP2-45, SSP4-60 and SSP5-85(Baseline)
    for i=1:4;    Tpath(i,:)=interp1([0:10:80],TEMP(i,:),[0:80], 'spline');end
    Tpath(:,[82:401])=Tpath(:,81)*ones(1,320); %constant temperature after 2100
    Tpath=Tpath';
%Paths
    year=[0:400]';
    GDP=GDP0*(exp(g*year));
    Discount=exp(-discount*year);
    Table=zeros(36,6);


for i=1:2 %gamma
    for j=1:4 %RCP
        for k=1:3 %varphi2
            for l=1:3 %phi+varphi
                for m=1:4 %lifetime
    RiskFactor=exp(-phivarphi(l)*year).*(1-exp(-varphi2(k).*year));
    SCC=sum(gamma(i).*zeta.*Tpath(tau1:end,j).*GDP(tau1:end).*Discount(tau1:end)); %gamma*zeta²*S=gamma*zeta*Temp
    SCO=sum(gamma(i).*zeta.*Tpath(tau1:tau2(m),j).*GDP(tau1:tau2(m)).*Discount(tau1:tau2(m)).*RiskFactor(tau1:tau2(m)));
    AdjFactor=transpose(SCO./SCC);
    Table(9*(j-1)+3*(k-1)+l,m)=AdjFactor;
    Table(9*(j-1)+3*(k-1)+l,5+(i-1))=SCC;
                end
            end
        end
    end
end
%risk adjusted same likelihood for RCP2.6, RCP3.4 and RCP6
meanT=mean(Tpath(tau1:80,1:3), 'all') % mean temperature is 2°C, so 

for i=1:2 %gamma
    %for j=1:3 %RCP
        for k=1:3 %varphi2
            for l=1:3 %phi+varphi
                for m=1:4 %lifetime
    RiskFactor=exp(-phivarphi(l).*(0.5+0.5/meanT*Tpath(tau1:end,1:3)).*year(tau1:end)).*(1-exp(-varphi2(k).*year(tau1:end)));
    SCC=sum(gamma(i).*zeta.*Tpath(tau1:end,1:3).*GDP(tau1:end).*Discount(tau1:end), 'all')/3; %gamma*zeta²*S=gamma*zeta*Temp
    SCO=sum(gamma(i).*zeta.*Tpath(tau1:tau2(m),1:3).*GDP(tau1:tau2(m)).*Discount(tau1:tau2(m)).*RiskFactor(tau1:tau2(m),:),'all')/3;
    AdjFactor=transpose(SCO./SCC);
    TableRisk(3*(k-1)+l,m)=AdjFactor;
    TableRisk(3*(k-1)+l,5+(i-1))=SCC;
                end
            end
        end
    %end
end
