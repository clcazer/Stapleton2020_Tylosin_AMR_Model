PK_Monte=1; %1=Monte Carlo of PK parameters, 0=deterministic
PD_Monte=1; %1=Monte Carlo of PD parameters, 0=deterministic
Pop_Monte=1; %1=Monte Carlo of Pop parameters, 0=deterministic
Growth=1; %bw growth is on

burn=50*24; %50 day burn-in
sim_time=221*24; %221 day simulation period, same as Schmidt 2020

%Neg control run with no TYL, no intervention
Treatment=5; %1=withdrawal, 2=antimicrobial-free pen, 
    %3=probiotic, 4=All, 5=no interventions
Method=0; %0=no TYL, 1=yes TYL
filename='NoTYL_validation_';
run('TYL_Model_cc.m');
save(strcat(filename,'cc.mat'));

%Pos control run with TYL, no intervention
Treatment=5; %1=withdrawal, 2=antimicrobial-free pen, 
    %3=probiotic, 4=All, 5=no interventions
Method=1; %0=no TYL, 1=yes TYL
filename='TYL_validation_';
run('TYL_Model_cc.m');
save(strcat(filename,'cc.mat'));