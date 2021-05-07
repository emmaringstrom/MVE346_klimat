
%% Uppgift 1.1
clc, clear, clf

%Parametrar
F = [0 60 0; 15 0 45; 45 0 0]; %GtC
NPP0 = F(1,2);
beta = 0.35; %Kan variera mellan 0.1 och 0.8
B0 = [600 600 1500]; %GtC
global CO2Emissions;
global CO2ConcRCP45;

%Nettoprimärproduktion av biomassa (fotosyntes)
NPP = @(B1) NPP0*(1+beta*log(B1/B0(1)));

%Flödeskoefficienter
alpha = @(i,j) F(i,j)/B0(i);

%Differentialekvationer
dB1 = @(B,t) alpha(3,1)*B(3,:) + alpha(2,1)*B(2,:) - NPP(B(1,:)) + CO2Emissions(t);
dB2 = @(B,t) NPP(B(1,:)) - alpha(2,3)*B(2,:) - alpha(2,1)*B(2,:);
dB3 = @(B,t) alpha(2,3)*B(2,:) - alpha(3,1)*B(3,:);

dB = @(B,t) [dB1(B,t); dB2(B,t); dB3(B,t)];

%Lösning
B = eulerforward(B0, dB);

%Koncentrationerna
conc = 0.469 * B(1,:); %ppm C02

%Plot
t = 1765:2500;
plot(t,conc,'m')
hold on
plot(t,CO2ConcRCP45,'g')
legend('Our estimation','RCP 4.5','Location','northwest')
xlabel('Year'), ylabel('CO_2 concentration (ppm)')


%% Uppgift 1.2

clc, clear, clf

%Parametrar
F = [0 60 0; 15 0 45; 45 0 0]; %GtC
NPP0 = F(1,2);
beta = linspace(0.1,0.8,5);
B0 = [600 600 1500]; %GtC
global CO2Emissions;
global CO2ConcRCP45;

col = ['m','g','b','r','c'];

for b = 1:length(beta)

    %Nettoprimärproduktion av biomassa (fotosyntes)
    NPP = @(B1) NPP0*(1+beta(b)*log(B1/B0(1)));

    %Flödeskoefficienter
    alpha = @(i,j) F(i,j)/B0(i);

    %Differentialekvationer
    dB1 = @(B,t) alpha(3,1)*B(3,:) + alpha(2,1)*B(2,:) - NPP(B(1,:)) + CO2Emissions(t);
    dB2 = @(B,t) NPP(B(1,:)) - alpha(2,3)*B(2,:) - alpha(2,1)*B(2,:);
    dB3 = @(B,t) alpha(2,3)*B(2,:) - alpha(3,1)*B(3,:);

    dB = @(B,t) [dB1(B,t); dB2(B,t); dB3(B,t)];

    %Lösning
    B = eulerforward(B0, dB);

    %Koncentrationerna
    conc = 0.469 * B(2,:); %ppm C02

    %Plot
    t = 1765:2500;
    txt = ['beta = ',num2str(beta(b))];
    plot(t,conc,col(b),'DisplayName',txt)
    hold on

end

%plot(t,CO2ConcRCP45,'k--','DisplayName','Actual','LineWidth',2)

hold off
legend('Location','northwest')


%% Uppgift 1.3

clc, clear, clf

%Parametrar
A = [0.113 0.213 0.258 0.273 0.1430];
tao_0 = [2.0 12.2 50.4 243.3 Inf];
k = 3.06*1e-3;
n = length(A);

%Tidskonstant
tao = @(U) tao_0*(1 + k*U);

%Impulssvaret
imp = @(t, tao) sum(A' .* exp(-t./tao'));

t = 0:500;
U = [0, 140, 560, 1680];

for u = U
    I = imp(t, tao(u));
    plot(t, I), hold on
end

legend('U = 0','U = 140','U = 560','U = 1680')

%% Uppgift 1.4
clc, clear, clf

%Parametrar
A = [0.113 0.213 0.258 0.273 0.1430];
tao_0 = [2.0 12.2 50.4 243.3 Inf];
k = 3.06*1e-3;
M0 = 630;
global CO2Emissions;
global CO2ConcRCP45;

U = CO2Emissions;

t = 1:length(U);

tao = @(t_hat) tao_0 .* (1 + k * cumsum(U))';
I = @(t, t_hat) sum(A' .* exp(-t./tao(t_hat)'));
M = @(t) M0 + cumsum(I(flip(t), t) .* U);

conc = 0.469 * M(t);

years = 1765:2500;
plot(years,conc,'DisplayName','Our model');
hold on
plot(years,CO2ConcRCP45(1:length(t)),'DisplayName','RCP4.5');
legend show
legend('Location','southeast')
xlabel('Year'), ylabel('CO_2 concentration (ppm)')






