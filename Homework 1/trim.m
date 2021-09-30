clc, clear all, close all

addpath('package');
load('package/AeroDB.mat');

Refl = 0.1524;
Refa = 0.01824;
Mass = 56.2951;
Iyy = 44.3618;

tblMachTrim = [ 0.6 0.85 0.95 1.05 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.6 ];
tblAoATrim  = [ 0   2    4    6    8   10  12  14  18  20  22  24  28  30  ];
tblAltTrim  = [ 0 5000 10000 15000 20000 25000 30000 ];

nMachTrim = length(tblMachTrim);
nAoATrim  = length(tblAoATrim);
nAltTrim = length(tblAltTrim);

tblFlight  = zeros(nAltTrim, nMachTrim, nAoATrim);
tblqTrim   = zeros(nAltTrim, nMachTrim, nAoATrim);
tbldelTrim = zeros(nAltTrim, nMachTrim, nAoATrim);
tblazTrim  = zeros(nAltTrim, nMachTrim, nAoATrim);

figure(1);

for i = 1:nAltTrim
    for j = 1:nMachTrim
        for k = 1:nAoATrim
            cAlt = tblAltTrim(i);
            cMach = tblMachTrim(j);
            cAoA = tblAoATrim(k);
            
            Cz0 = interp2(tblMach, tblAoA, tblCz0', cMach, cAoA);
            Czd = interp2(tblMach, tblAoA, tblCzd', cMach, cAoA);
            CM0 = interp2(tblMach, tblAoA, tblCM0', cMach, cAoA);
            CMd = interp2(tblMach, tblAoA, tblCMd', cMach, cAoA);
            CMq = interp1(tblMach, tblCMq, cMach);
            
            [Tmp, SoS, ~, Rho] = atmosisa(cAlt);
            
            V = cMach * SoS;
            Q = 0.5 * Rho * V^2;
            
            qTrim = (CM0-(CMd/Czd)*Cz0)/((CMd/Czd)*(Mass*V/(Q*Refa))-CMq*Refl/(2*V));
            delTrim = (1/CMd)*(-CM0-CMq*qTrim*(Refl/(2*V)));
            azTrim = (Q*Refa/Mass)*(Cz0+Czd*delTrim);
            
            tblqTrim(i,j,k) = qTrim*180/pi;
            tbldelTrim(i,j,k) = delTrim*180/pi;
            tblazTrim(i,j,k) = azTrim/9.81;
            
            if abs(tbldelTrim(i,j,k)) <= 15
                tblFlight(i,j,k) = 1;
                
                plot3(cMach, cAlt/1000, cAoA, 'r*', 'linewidth', 2.0);
                hold on;
            else
                plot3(cMach, cAlt/1000, cAoA, 'b*', 'linewidth', 2.0);
                hold on;
            end
        end
    end
end

xlabel('Mach');
ylabel('H [km]');
zlabel('Az [g]');
grid on;