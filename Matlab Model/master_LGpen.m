%run this after burn-in calculation and validation

PK_Monte=1; %1=Monte Carlo of PK parameters, 0=deterministic
PD_Monte=1; %1=Monte Carlo of PD parameters, 0=deterministic
Pop_Monte=1; %1=Monte Carlo of Pop parameters, 0=deterministic
Growth=1; %bw growth is on

burn=50*24; %50 day burn_in
sim_time=143*24; %143 day simulation period

%Neg control run with no TYL, no intervention
Treatment=5; %1=withdrawal, 2=antimicrobial-free pen, 
    %3=probiotic, 4=All, 5=no interventions
Method=0; %0=no TYL, 1=yes TYL
filename='NoTYL_NI_LGpen_';
run('TYL_Model_cc_LGpen.m');
%save(strcat(filename,'cc.mat'));

%Pos control run with TYL, no intervention
Treatment=5; %1=withdrawal, 2=antimicrobial-free pen, 
    %3=probiotic, 4=All, 5=no interventions
Method=1; %0=no TYL, 1=yes TYL
filename='TYL_NI_LGpen_';
run('TYL_Model_cc_LGpen.m');
%save(strcat(filename,'cc.mat'));

%for withdrawal, no neg control because there can be no withdrawal without TYL
%TYL, withdrawal
Treatment=1; %1=withdrawal, 2=antimicrobial-free pen, 
    %3=probiotic, 4=All, 5=no interventions
Method=1; %0=no TYL, 1=yes TYL
filename='TYL_RWT_LGpen_';
run('TYL_Model_cc_LGpen.m');
%save(strcat(filename,'cc.mat'));

%No TYL, AFTP
Treatment=2; %1=withdrawal, 2=antimicrobial-free pen, 
    %3=probiotic, 4=All, 5=no interventions
Method=0; %0=no TYL, 1=yes TYL
filename='NoTYL_AFTP_LGpen_';
run('TYL_Model_cc_LGpen.m');
%save(strcat(filename,'cc.mat'));

%TYL, AFTP
Treatment=2; %1=withdrawal, 2=antimicrobial-free pen, 
    %3=probiotic, 4=All, 5=no interventions
Method=1; %0=no TYL, 1=yes TYL
filename='TYL_AFTP_LGpen_';
run('TYL_Model_cc_LGpen.m');
%save(strcat(filename,'cc.mat'));

%No TYL, probiotic
Treatment=3; %1=withdrawal, 2=antimicrobial-free pen, 
    %3=probiotic, 4=All, 5=no interventions
Method=0; %0=no TYL, 1=yes TYL
filename='NoTYL_DFM_LGpen_';
run('TYL_Model_cc_LGpen.m');
%save(strcat(filename,'cc.mat'));

%TYL, probiotic
Treatment=3; %1=withdrawal, 2=antimicrobial-free pen, 
    %3=probiotic, 4=All, 5=no interventions
Method=1; %0=no TYL, 1=yes TYL
filename='TYL_DFM_LGpen_';
run('TYL_Model_cc_LGpen.m');
%save(strcat(filename,'cc.mat'));
