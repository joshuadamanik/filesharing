clc, clear all, close all

addpath('package');
load('package/AeroDB.mat');

Refl = 0.1524;
Refa = 0.01824;
Mass = 56.2951;
Iyy = 44.3618;

MachTrim = 2.0;
AoATrim = 0.0;
AltTrim = 0.0;
delAoA = 0.05;
            
Cz0  = interp2(tblMach, tblAoA, tblCz0', MachTrim, AoATrim);
Czd  = interp2(tblMach, tblAoA, tblCzd', MachTrim, AoATrim);
CM0  = interp2(tblMach, tblAoA, tblCM0', MachTrim, AoATrim);
CMd  = interp2(tblMach, tblAoA, tblCMd', MachTrim, AoATrim);
CMq  = interp1(tblMach, tblCMq, MachTrim);
Cz0p = interp2(tblMach, tblAoA, tblCz0', MachTrim, AoATrim+delAoA);
CM0p  = interp2(tblMach, tblAoA, tblCM0', MachTrim, AoATrim+delAoA);
            
[Tmp, SoS, ~, Rho] = atmosisa(AltTrim);
            
V = MachTrim * SoS;
Q = 0.5 * Rho * V^2;
            
qTrim = (CM0-(CMd/Czd)*Cz0)/((CMd/Czd)*(Mass*V/(Q*Refa))-CMq*Refl/(2*V));
delTrim = (1/CMd)*(-CM0-CMq*qTrim*(Refl/(2*V)));
azTrim = (Q*Refa/Mass)*(Cz0+Czd*delTrim);

Cz0a = (Cz0p-Cz0)/(delAoA*pi/180);
CM0a = (CM0p-CM0)/(delAoA*pi/180);

Za = (Q*Refa/(Mass*V))*Cz0a;
Zd = (Q*Refa/(Mass*V))*Czd;
Ma = (Q*Refa*Refl/Iyy)*CM0a;
Md = (Q*Refa*Refl/Iyy)*CMd;
Mq = (Q*Refa*Refl^2/(2*V*Iyy))*CMq;
Aa = V*Za;
Ad = V*Zd;

% Model without IMU shift
q_d = tf([Md, Ma*Zd-Md*Za], [1, -(Za+Mq), Za*Mq-Ma]);
figure, pzplot(q_d), title('Center of Gravity');
saveas(gcf, 'fig\q_cg.png');
close all;

az_d = tf([Ad, Aa*Zd-(Za+Mq)*Ad, Aa*Md-Ad*Ma+(Ad*Za-Aa*Zd)*Mq], [1, -(Za+Mq), Za*Mq-Ma]);
figure, pzplot(az_d), title('Center of Gravity');
saveas(gcf, 'fig\az_cg.png');
close all;

% Model with IMU shift
for Xi = -2:0.2:2
    az_ds = tf([Ad-Md*Xi, Aa*Zd-(Za+Mq)*Ad+(Md*Za-Ma*Zd)*Xi, Aa*Md-Ad*Ma+(Ad*Za-Aa*Zd)*Mq], [1, -(Za+Mq), Za*Mq-Ma]);
    figure, pzplot(az_ds), title(strcat('Shifted, X_{IMU} = ', num2str(Xi)));
    saveas(gcf, strcat('fig\az_', num2str(Xi), '.png'));
    close all;
end
