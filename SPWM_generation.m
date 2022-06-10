% #########################################################################
%
% Creating PWL Files for Sinusoidal PWM to use it in LTspice
%
% This script generates PWM signals for high side and low side to source
% in LTspice. Its inputs consist of rise time, fall time, fundamental
% frequency, switching frequency, deadtime, modulation limit, Vdd, Vcc.
%
% Its features are symmetric deadtime, multi-phase generation with
% different phase shifts, adjustable rise time, fall time, deadband,
% resolution (time step).
%
% Written by ozgur gulsuna
%            samet yakut
%            ısık emir altunkol
%                                                           29.10.2021
% #########################################################################


clc;
clear;
close all;

%% INPUTS
f=200000;                        % switching Frequency
n=350;                              % total number of cycles
Tr=10e-9;                          % rise time
Tf=10e-9;                          % fall time
timeStep =5e-9;                 % minimum time step

waveType = 1 ;                    % 0 = Constant duty
                                  % 1 = Sine

duty = 0.3 ;                      % for the constant duty waveType
fundamental = 2000 ;              % fundamental frequency for the sine
deadBand = 25e-9 ;                % deadtime 

Vdd = 3;                        % 12V supply
Vcc = 0;                        % reference of the Vdd
modulation = 0.95;               % modulation index  
range = 1-2*deadBand*f ;          % duty range to limit duty due deadtime
phaseA = 180 ;                    % phase difference in degree
phaseB = 210 ;
phaseC = 0 ;
% phase difference in degree


%% Calculations
totalTime=n/fundamental;                            % Total simulation time
time = 0:timeStep:totalTime-timeStep;               % time array in seconds

trigRef= 0.5*sawtooth(2*pi*f*time,1/2)+0.5;         % Triangle reference
sineRefA = modulation*0.5*range*sin(2*pi*fundamental*time+phaseC*pi/180)+modulation*0.5+deadBand*f;                % 0 degrees
sineRefB = modulation*0.5*range*sin(2*pi*fundamental*time+phaseA*pi/180)+modulation*0.5+deadBand*f;  % 120 degrees
sineRefC = modulation*0.5*range*sin(2*pi*fundamental*time+phaseB*pi/180)+modulation*0.5+deadBand*f;  % 240 degrees

arrayLength = n*f*4/fundamental +2;


%% Initilization

time1  = zeros(1,arrayLength);
value1 = zeros(1,arrayLength);
time2  = zeros(1,arrayLength);
value2 = zeros(1,arrayLength);

Ref = sineRefA;

t=3;                                % time counter

%% PWL Generation

for p=1:3
    for i=2:length(time)
        if ((trigRef(i)-Ref(i)>=0) && (trigRef(i-1)-Ref(i-1)<=0)) % positive slope triangle part
            time1(t)=time(i)+deadBand/2;
            time1(t-1)=time(i)-Tr+deadBand/2;        %high side - rising edge
            value1(t)=Vdd;
            value1(t-1)=Vcc;
            time2(t-1)=time(i)-deadBand/2;
            time2(t)=time(i)+Tf-deadBand/2;          %low side - falling edge
            value2(t)=Vcc;
            value2(t-1)=Vdd;
            t=t+2;

        end
        if((trigRef(i)-Ref(i)<0) && (trigRef(i-1)-Ref(i-1)>0))  % negative slope triangle part
            time1(t-1)=time(i)-deadBand/2;
            time1(t)=time(i)+Tf-deadBand/2;         %high side - falling edge
            value1(t)=Vcc;
            value1(t-1)=Vdd;
            time2(t)=time(i)+deadBand/2;
            time2(t-1)=time(i)-Tr+deadBand/2;       %low side - rising edge
            value2(t)=Vdd;
            value2(t-1)=Vcc;
            t=t+2;
        end
    end
    t=3;   % reset the time counter

    time1(arrayLength) = totalTime;                  % deal with last point
    value1(arrayLength) = value1(arrayLength-1);    

    time2(arrayLength) = totalTime;
    value2(arrayLength) = value2(arrayLength-1);

    if p == 1
        PhaseA_H  = [time1.' value1.'];
        PhaseA_L  = [time2.' value2.'];

        figure;
        plot(time1,value1);
        hold on
        plot(time,Ref);
        plot(time2,value2);

        Ref = sineRefB;
    end
    if p == 2
        PhaseB_H = [time1.' value1.'];
        PhaseB_L = [time2.' value2.'];


        figure;
     plot(time1,value1);
        hold on
        plot(time,Ref);
        plot(time2,value2);

        Ref = sineRefC;
    end
    if p == 3

        figure;
        plot(time1,value1);
        hold on
        plot(time,Ref);
        plot(time2,value2);
        PhaseC_H = [time1.' value1.'];
        PhaseC_L = [time2.' value2.'];

        Ref = sineRefA;
    end

end

%% Writing PWL files
writematrix(PhaseA_H,'PhaseA-H.txt');
writematrix(PhaseA_L,'PhaseA-L.txt');
writematrix(PhaseB_H,'PhaseB-H.txt');
writematrix(PhaseB_L,'PhaseB-L.txt');
writematrix(PhaseC_H,'PhaseC-H.txt');
writematrix(PhaseC_L,'PhaseC-L.txt');
