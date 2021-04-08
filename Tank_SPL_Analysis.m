% Alexandre Otávio Guedes Drummond         
% 05/04/2021
%%
clc, clear all, close all
% 

%% System's Constants

rho0= 1.2922;                          %air density [kg/m³]
c0 = 343.0;                            %air sound speed[m/s]
ref = 20e-6;                           %pressure reference [Pa] 
rho1 = 1000.0;                         %water density [kg/m³]
c1 = 1480.0;                           % water sound speed[m/s]

% Subsystems dimensions
r1=5.65; r2=4.1; r3=3.2 
t1=3.962; t2=3.048 ;t3=2.59 ;
Vr = r1*r2*r3;                                %[m³] room's volume  
Vt = t1*t2*t3;                                     %[m³] tank's volume 

Sroom = 2*(r1*r2) + 2*(r1*r3) + 2*(r2*r3) - (t1*t3);      % [m²] Room's surface area
Lroom =  4*r1+4*r2+4*r3;                           % [m] 
Stank =  2*(t1*t2)+2*(t2*t3) + 2*(t1*t3); 
Ltank =  4*t1+4*t2+4*t3; 
Atank = t1*t3;                                    % Tank's surface area

%% Reverberation Time

Paper_T60_r = [8.61 7.54 7.76 8.69 10.43 10.65 8.92 8.59 7.98 6.74 6.07 5.34 4.96 4.41 3.73 3.05 2.68 2.31];
Paper_T60_t = [0 0 0 0 0 0 0 0 .35 .42 .43 .33 .32 .25 .26 .24 .25 .25];

%% Vacuum Sound Power
We = [1.067701213366e-7 3.258485333e-6 2.396535e-7 4.50449e-7 1.349024e-6 2.190076e-6 3.531653e-6 5.11264e-6 0.00024915652 3.21072299e-5 7.0694245e-6 7.5555327e-6 ...
      1.1512413e-5 1.489275e-5 1.7656767e-5 8.76629240e-6 6.7690471e-6 6.984186e-6];

%% Getting pump's sound power

W0 = 10e-12;                               %Reference power
Power = sum(We);                           %Total Power 
Original_SWL = 10*log10(abs((Power)/W0));

We = We*42; %Power Compensation for the power to be 92 dB

Adapted_Power = sum(We);
Adapted_SWL = 10*log10(abs((Adapted_Power)/W0))

%Since Original_SWL is 151.56 dB and I know the required power to be 92 dB,
% the power must be 1000x smaller, corresponding to 60 dB less

%% Schroeder's Frequency

Fsr = 2000*sqrt(mean(Paper_T60_r)/Vr);                 %Rooms Schroeder's frequency
Fst = 2000*sqrt(mean(Paper_T60_t)/Vt);                 %Tans Schroeder's frequency

%% Array of Frequencies

freq_bands = [100 125 160 200 250 315 400 500 630 800 1000 1250 1600 ...
     2000 2500 3150 4000 5000] ;        %Frequencies (one-third octave bands)

for i = drange(1, length(freq_bands))
    if freq_bands(i) >= Fsr
       break
    end 
end
f = freq_bands(i-1:end)                           %Use only frequencies from Schroeder's frequency
w = 2*pi*f;                                       %From linear to angular freq.

T60_r = Paper_T60_r(i-1:end);
T60_t = Paper_T60_t(i-1:end);
%% Energies Calculation
for i = drange(1,length(f))
    
    %Losses Factors
    
    etar= (2.2/(f(i)*T60_r(i))) ;
    m_t = (rho1*Vt)/(t1*t3);        
    a_t = (m_t*t1*w(i))/(2*rho1*c1) ;
    tau_rt = log(1+(a_t^2))./(a_t^2);
    etart= ((c1*Atank)/(4*w(i)*Vr))*tau_rt ;
    etat= (2.2/(f(i)*T60_t(i)));
    m_r = (rho0*Vr)/(t1*t3); 
    a_r = (m_r*t3*w(i))/(2*rho0*c0);
    tau_tr = log(1+(a_r^2))./(a_r^2);
    etatr= ((c0*Atank)/(4*w(i)*Vt))*tau_tr ;



    ETA = [-etart (etat + etatr) ; ...
            (etar+etart)   -etatr ]   %coefficients matrix

    W = [0 We(i)];                     %injected power matrix
    
   Energy_Mat(:,i) = (w(i)*ETA)\W';  %Solve the system (one iteration per frequency band) 
end

% Energies matrix stores each subsystems energy, per frequency band.

 
%% Spatial Average sound pressure
Pref1 = 2e-5;                                   %Air Sound pressure reference
Pref2 = 1e-6;                                   %Underwater sound pressure reference
Pr = sqrt((Energy_Mat(1,:)*rho0*(c0^2))/Vr);   %Average pressure inside the room
Pt = sqrt((Energy_Mat(2,:)*rho1*(c1^2))/Vt);   %Average pressure inside the tank

%To dB
Pr_dB = 20*log10(abs(Pr)/Pref1);
Pt_dB = 20*log10(abs(Pt)/Pref2);        



%% Plots


%Rooms sound pressure
figure()
mean_SPL = mean(Pr_dB);
mean_SPL = mean_SPL*ones(length(f),1);
total_SPL = 20*log10(abs(sum(Pr)/Pref1))
semilogx(f,Pr_dB,'k')
txt = ['Spatial Average SPL inside the room (TOTAL =' num2str(total_SPL, '%0.0f') 'dB)']
title(txt)
xlabel('Frequency Bands [Hz]')
ylabel('SPL [dB](Ref:20e-6)')
xticks(f)
xlim([630 5000])
xticklabels({'500','630','800','1k','1.25k','1.6k','2k','2.5k','3.15k','4k','5k'})
yticks([70, 75, 80, mean_SPL(1), 85, 90, 95, 100, 105]);
yticklabels({'70','75', '80', num2str(mean_SPL(1),'\\bf%0.0f'), '85', '90', '95', '100', '105'})
hold on;
semilogx(f,mean_SPL,'b--','LineWidth',2.5)
legend('SPL', 'Mean per Frequency')
grid minor


figure()
mean_SPL = mean(Pt_dB(2:end));
mean_SPL = mean_SPL*ones(length(f),1);
total_SPL = 20*log10(abs(sum(Pt)/Pref2))
semilogx(f,Pt_dB,'k')
txt = ['Spatial Average SPL inside the tank (TOTAL = ' num2str(total_SPL,'%0.0f') 'dB)']
title(txt)
xlabel('Frequency Bands [Hz]')
ylabel('SPL [dB](Ref:1e-6)')
xticks(f)
xlim([630 5000])
xticklabels({'500','630','800','1k','1.25k','1.6k','2k','2.5k','3.15k','4k','5k'})
yticks([115, 120, 125, 130, mean_SPL(1), 135, 140]);
yticklabels({'115','120', '125', '130', num2str(mean_SPL(1),'\\bf%0.0f'), '135', '140'})
hold on;
semilogx(f,mean_SPL,'b--','LineWidth',2.5)
legend('SPL', 'Mean per Frequency')
text(2500,125,txt)
grid



%%

%Source Sound Power

W0 = 10e-12;
We_NWS = 10*log10(abs((We)/W0));

figure()
semilogx(freq_bands, We_NWS,'m','LineWidth',2.0)
title("Pump's Sound Power")
xlabel('Frequency Bands [Hz]')
ylabel('Sound Power [W]')
xticks(freq_bands)
xticklabels({'100','','160','','','315','','','630','','1k','','','2k','','','4k','5k'})
xlim([100,5000])
grid on

%%

%Room's Reverberation Time

figure()
semilogx(freq_bands, Paper_T60_r,'k','LineWidth',1.0)
title("Room's Reverberation Time")
xlabel('Frequency Bands [Hz]')
ylabel('T60[s]')
xticks(freq_bands)
xticklabels({'100','','160','','','315','','','630','','1k','','','2k','','','4k','5k'})
xlim([100,5000])
grid on 
%Tank's Reverberation Time

figure()
semilogx(freq_bands, Paper_T60_t,'k','LineWidth',1.0)
title("Tank's Reverberation Time")
xlabel('Frequency Bands [Hz]')
ylabel('T60 [s]')
xticks(freq_bands(9:end))
xticklabels({'630','','1k','','','2k','','','4k','5k'})
xlim([630,5000])
grid on 