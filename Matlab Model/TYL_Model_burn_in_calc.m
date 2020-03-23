PK_Monte=1; %1=Monte Carlo of PK parameters, 0=deterministic
PD_Monte=1; %1=Monte Carlo of PD parameters, 0=deterministic
Pop_Monte=1; %1=Monte Carlo of Pop parameters, 0=deterministic

%establish burn-in period. turned off BW growth
%Neg control run with no TYL, no intervention
Treatment=5; %1=withdrawal, 2=antimicrobial-free pen, 
    %3=probiotic, 4=All, 5=no interventions
Method=0; %0=no TYL, 1=yes TYL
Growth=0;
burn=143*24; %entire simulation of burn-time so that feed/water compartments will not be cleaned
sim_time=0;
run('TYL_Model_cc.m'); %write matrices at end off


%confirm that model behaves
%final Ent concentrations (after 145d) below carrying capacity in cattle
sum(Cow_total_conc(end,:)>K_c, 'all') %ok--0 sims with cattle ent conc > carrying capacity
%but Ent conc are not always below carrying capacity
sum(Cow_total_conc(240,:)>K_c, 'all') %362 iterations out of 100 are above cc at day 1

%note that carrying capacity is not used to limit Ent population in other compartments (F/W/P), just to set initial values
%so check that maximum concentration in those compartments is not above plausible values
max(Feed_total_conc, [], 'all')<max(K_f, [], 'all') %ok
max(Water_total_conc, [], 'all')<max(K_w, [], 'all') %not all below K_w
max(Water_total_conc, [], 'all') %max is 5.7e4. plausible for highly contaminated water
max(Pen_total_conc, [], 'all')<max(K_p, [], 'all') %not all below K_p
max(Pen_total_conc, [], 'all') %max is 1.4e7. Plausible for manure

%how many iterations have Ent population hit 0?
sum(sum(Cow_total==0)~=0) %1 iterations have population hit 0
sum(Cow_total(end,:)==0) %but not at 0 at end of simulation

%burn in period should result in all Cow_total_conc below carrying capacity
Below_CC=Cow_total_conc<K_c;
num_iter_below_cc=sum(Below_CC,2); %note that first timepoint has all iterations below cc because of starting point
all_iter_below_cc=find(num_iter_below_cc(2:end)==1000);
min(all_iter_below_cc) %11574
%minimum days to burn-in
ceil(min(all_iter_below_cc)/240) %49 days

%gradient should be low (little change over time)
%accept 90% of iterations with slope less than 1 CFU/g change per hour (?) = 0.1 CFU/g per timestep
slope=abs(diff(Cow_total_conc,1,1));
slope_less=slope<0.1;
num_iter_slope_less=sum(slope_less,2);
pct90_iter_slope_less=find(num_iter_slope_less>899);
ceil(min(pct90_iter_slope_less)/240) %42 days

%50 day burn-in will suffice
%at this point, 93.6% of iterations have low slope
num_iter_slope_less(50*240)





