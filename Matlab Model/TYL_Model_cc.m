
%{
Abbreviations
PK = pharmacokinetic
PD = pharmacodynamic
Pop = population
TYL = tylosin phosphate
USI = upper small intestine
LSI = lower small intestine
LI = large intestine
environ = environment
Ent = Enterococcus
res = resistant
sus = susceptible
conc = concentration
amt = amount
prop = proportion
cc = carrying capacity
T = time
%}

%must define outside of this script:
%burn: hours of burn-in period
%sim_time: hours of simulation after burn-in period
%PD_Monte, PK_Monte, Pop_Monte: logical for 1 simulation (=0) or monte carlo 1000 sims (=1)
%Growth: logical for no body weight increase (=0) or body weight increase (=1)
%Treatment:  1=withdrawal, 2=antimicrobial-free pen, 3=probiotic, 4=All, 5=no interventions
%Method:  0=no TYL, 1=yes TYL

rng(1); %seed set for model development and validation

time=burn+sim_time+0.1; %length of simulation is burn-in plus experiment sim time plus 0.1 because time start at 0: hours

dt=0.1; %time step of 0.1h: hours
    
TYL_start_time=burn; %time when TYL will be given in the feed: hours

if Treatment==1 || Treatment==2 || Treatment==4 %for RWT and AFTP, stop TYL 30 days earlier
    TYL_stop_time=TYL_start_time+sim_time-(30*24); %hours
else
    TYL_stop_time=TYL_start_time+sim_time; %hours
end

if Treatment==2 || Treatment==4 %for AFTP
    AFTP_time=TYL_start_time+sim_time-(30*24); %move to AFTP 30 days prior to end: hours
end

if PK_Monte==0 && PD_Monte==0 && Pop_Monte==0
    n=1; %deterministic model number of simulations
else
    n=1000; %Monte Carlo model number of simulations
end

%Storage arrays for TYL LI concentration & compartment TYL amounts
TYL_li_conc=zeros(time/dt,n); %Store conc of TYL in L
TYL_stomach_amt=zeros(time/dt,n); %Store amt of TYL in stomach
TYL_USI_amt=zeros(time/dt,n); %Store amt of TYL in USI
TYL_LSI_amt=zeros(time/dt,n); %Store amt of TYL in LSI
TYL_LI_amt=zeros(time/dt,n); %Store amt of TYL in LI

%Storage arrays for TYL PD effect 
Er=zeros(time/dt,n); %Store PD effect on res Ent
Es=zeros(time/dt,n); %Store PD effect on sus Ent

%Storage arrays for number of Ent (R vs S) in the different populations
Cow_res=zeros(time/dt,n); %Store res Ent in cow 
Cow_sus=zeros(time/dt,n); %Store sus Ent in cow
Feed_res=zeros(time/dt,n); %Store res Ent in feed
Feed_sus=zeros(time/dt,n); %Store sus Ent in feed
Water_res=zeros(time/dt,n); %Store res Ent in water
Water_sus=zeros(time/dt,n); %Store sus Ent in water
Pen_res=zeros(time/dt,n); %Store res Ent in pen
Pen_sus=zeros(time/dt,n); %Store sus Ent in pen

%Storage arrays of total number Ent in each compartment
Cow_total=zeros(time/dt,n); %Store total Ent in cattle 
Feed_total=zeros(time/dt,n); %Store total Ent in feed
Water_total=zeros(time/dt,n); %Store total Ent in water
Pen_total=zeros(time/dt,n); %Store total Ent in pen

%Storage arrays for concentration of Ent in each compartment
Cow_total_conc=zeros(time/dt,n); %store CFU/g-Ingesta
Feed_total_conc=zeros(time/dt,n); %store CFU/g-Feed
Water_total_conc=zeros(time/dt,n); %store CFU/g-Water
Pen_total_conc=zeros(time/dt,n); %store CFU/g-Pen

%PK model parameters
    delta=random(makedist('Lognormal', 'mu', -5.7, 'sigma', 0.7), 1, n); %degradation rate: 1/hr
    lambda_s=unifrnd(0.05,0.09,1,n); %movement rate out stom: 1/hr
    lambda_usi=unifrnd(0.2,0.4,1,n); %movement rate out of USI: 1/hr
    lambda_lsi=unifrnd(0.1,0.2,1,n); %movement rate out of LSI: 1/hr
    lambda_li=unifrnd(0.1,0.2,1,n); %movement rate out of LI: 1/hr
    mu=random(truncate(makedist('Normal', 'mu', 0.7, 'sigma', 0.1)...
        , 0, 1), 1, n); %sorption of TYL to digesta as a proportion: unitless
    free=1-mu; %calculate proportion of TYL free given proportion sorbed: unitless
    initial_bw=300; %initial body weight of steers: kg

%PD model parameters
    Emax=1; %E_Max constant at 1: unitless
    E_0=1; %E_0 constant at 1: unitless
    aero=abs(normrnd(0,1.2,1,n)); %accounting for anaerobicity in LI: unitless
    log2_MIC_sus=random(truncate(makedist('Normal', 'mu', 1, 'sigma', 0.8), -3, 4), 1, n); %sample MIC from log2MIC distribution
    MIC_sus=repmat(2.^(log2_MIC_sus + aero),1);  %MIC of sus Ent to TYL: ug/mL
    log2_MIC_res=unifrnd(4,7,1,n); %sample MIC from log2MIC distribution
    MIC_res=repmat(2.^(log2_MIC_res + aero),1); %MIC of res Ent to TYL: ug/mL
    Hill_sus=unifrnd(1.3,2.1,1,n); %Hill Coefficient for sus Ent: unitless
    Hill_res=unifrnd(2.6,4.3,1,n); %Hill Coefficient for res Ent: unitless

%Pop model parameters
    N_c=repmat(11,1,n); %number of cattle in the pen: head
    W_c=random(truncate(makedist('Lognormal', 'mu', 7.5, 'sigma', 0.2), 500, 4000), 1, n); %amt cattle drink: g-Water/(head hr)
    F_c=unifrnd(458,517,1,n); %amt of feed cattle eat: g-Feed/(head hr)
    P_c=unifrnd(10,30,1,n); %amt of pen environ cattle ingest: g-Pen/(head hr)
    DFM=unifrnd(0,0.4,1,n); %effect of DFM on K_c in cattle: log10CFU/g
    log10_K_c=random(truncate(makedist('Weibull', 'a', 5.7, 'b', 5.9), 2,7),1,n); %log carrying capacity, a is scale; b is shape: log10CFU/g
    if Treatment==3 || Treatment==4 %for DFM treatments
        %log carrying capacity is reduced by DFM.
        K_c=repmat((10.^(log10_K_c-DFM)),1); %carrying capacity of Ent in cattle LI: CFU/g-Ingesta
    else
         K_c=repmat((10.^log10_K_c),1); %carrying capacity of Ent in cattle LI: CFU/g-Ingesta
    end
    plasmid=random(truncate(makedist('Lognormal', 'mu', -12.2, 'sigma', 4.7), 0, 0.01), 1, n); %plasmid transfer rate from res to sus: 1/hr
    R_c=random(truncate(makedist('Lognormal', 'mu', -2.1, 'sigma'...
        , 1.1),0, 1), 1,n); %rate of growth of Ent in cattle: 1/hr
    Trough=unifrnd(26500,530000,1,n); %size of water trough: g-Water
    log10_K_w=unifrnd(1.2,2.1,1,n); %sample water cc in log10
    K_w=repmat((10.^log10_K_w),1); %cc of Ent in water: CFU/mL
    R_w=unifrnd(0.03,0.05,1,n); %rate of death of Ent in water: 1/hr
    Bunk=unifrnd(100000,170000,1,n); %size of feed bunk: g-Feed
    log10_K_f=random(makedist('Triangular', 'a', 2, 'b', 4.4, 'c', 6.7), 1, n); %sample feed cc in log10
    K_f=repmat((10.^log10_K_f),1); %cc of Ent in feed: CFU/g-Feed
    R_f=unifrnd(0.01,0.02,1,n); %rate of death of Ent in feed: 1/hr
    Ingested_death=unifrnd(0.5,0.9,1,n); %death penalty to ingested bacteria moving from feed, water, and 
        %pen to cattle. Proportion of ingested bacteria that die:
        %unitless
    Pen=unifrnd(30000,50000,1,n); %amt of ingestible pen material: g-Pen
    log10_K_p=random(truncate(makedist('Weibull', 'a', 5.9, 'b'...
        , 4), 2, 7), 1, n); %scale=a; shape=b. log10CFU/g
    K_p=repmat((10.^log10_K_p),1); %cc of Ent in pen: CFU/g-Pen
    R_p=random(truncate(makedist('Lognormal', 'mu', -5.5, 'sigma', 0.8), 0,1), 1,n); %rate of death of Ent in pen environment: 1/hr
    Fit=random(truncate(makedist('Weibull', 'a', 0.03, 'b', 1.7),0,0.1),1,n);%fitness cost for bearing res genes as proportion: unitless
    W_to_P=unifrnd(0.003,0.005,1,n); %movement rate of Ent from water to pen: 1/hr
    P_to_W=unifrnd(0.001,0.003,1,n); %movement rate of Ent from pen to water: 1/hr
    F_to_P=unifrnd(0.003,0.005,1,n); %movement rate of Ent from feed to pen: 1/hr
    P_to_F=unifrnd(0.0001,0.0003,1,n); %movement rate of Ent from pen to feed: 1/hr
    C_to_P=unifrnd(0.01,0.02,1,n); %movement rate of Ent from cattle to pen (defecation): 1/hr
    W_to_C=N_c.*W_c./Trough;  %movement rate of water to cattle (drinking) 1/hr
    F_to_C=N_c.*F_c./Bunk; %movement rate of feed to cattle (eating): 1/hr
    P_to_C=N_c.*P_c./Pen; %movement rate of pen to cattle (ingestion): 1/hr
    Zc=unifrnd(0.1,0.9,1,n);
    Zw=unifrnd(0.1,0.9,1,n);
    Zf=unifrnd(0.1,0.9,1,n);
    Zp=unifrnd(0.1,0.9,1,n);
        %Z is starting amount of Ent expressed as proportion of carrying
        %capacity: unitless
    Yc=random(makedist('Beta', 'a', 0.4, 'b', 4.4),1,n); %a: first shape; b: second shape
    Yw=random(makedist('Beta', 'a', 0.4, 'b', 4.4),1,n);
    Yf=random(makedist('Beta', 'a', 0.4, 'b', 4.4),1,n);
    Yp=random(makedist('Beta', 'a', 0.4, 'b', 4.4),1,n);
        %Y is proportion of starting Ent that are resistant: unitless

%Time, body weight and LI volume are deterministic and the same for each simulation
T=(0:dt:(time-dt))'; %T (time) variable is in increments of 0.1: hours

%bw is only used to calculate V_li
if Growth==0 %if growth is off
    bw=repmat(initial_bw,time/dt,1); %no bw increase
else %if growth is on
    bw1=repmat(initial_bw, burn/dt, 1); %no bw increase during burn-in
    bw2=initial_bw+...
        (0.192676.*(T((burn/dt)+1:end,1)/(24)))+(0.00514.*((T((burn/dt)+1:end,1)/(24)).^2))...
        -((6.38*(10^-6)).*((T((burn/dt)+1:end,1)/(24)).^3)); %body weight increase: kg
    bw=vertcat(bw1, bw2); %kg
end

LI=10^(-0.936).*(bw.^1.032)*.079*1000; %LI content: g-Ingesta/head
V_li=(LI/1000)*1.99*1.5; %store volume of cattle LI: Liters/head
    
for j=1:n %for each simulation
%PK compartments that are overwritten with each model simulation
TYL_feed=zeros(time/dt,1); %storage array for TYL dosage in the feed: mg/hr

%initialize PK compartments. PK is modeled in one animal (head = 1)
TYL_feed(1,1)=0; %starting rate of TYL in feed: mg/hr
TYL_stomach_amt(1,j)=0; %starting amt of TYL in stomach assigned to array: mg
TYL_USI_amt(1,j)=0; %starting amt of TYL in USI assigned to array: mg
TYL_LSI_amt(1,j)=0; %starting amt of TYL in LSI assigned to array: mg
TYL_LI_amt(1,j)=0; %starting amt of TYL in LI assigned to array: mg
TYL_li_conc(1,j)=0; %starting conc of TYL in LI: mg/L == ug/mL

%initialize PD compartments
Er(1,j)=1; %starting PD effect for res Ent: unitless
Es(1,j)=1; %starting PD effect for sus Ent: unitless

%initialize Monte Carlo Pop compartments
%each compartment is the number of Ent CFU (not concentration /g)
Cow_total(1,j)=Zc(1,j)*K_c(1,j)*LI(1,1)*N_c(1,j); %starting total Ent in cattle: CFU
Feed_total(1,j)=Zf(1,j)*K_f(1,j)*Bunk(1,j); %starting total Ent in feed: CFU
Water_total(1,j)=Zw(1,j)*K_w(1,j)*Trough(1,j); %starting total Ent in water: CFU
Pen_total(1,j)=Zp(1,j)*K_p(1,j)*Pen(1,j); %starting total Ent in pen: CFU

Cow_res(1,j)=Cow_total(1,j)*Yc(1,j); %starting res Ent in each cattle: CFU
Cow_sus(1,j)=Cow_total(1,j)*(1-Yc(1,j)); %starting sus Ent in each cattle: CFU
Feed_res(1,j)=Feed_total(1,j)*Yf(1,j); %starting res Ent in feed: CFU
Feed_sus(1,j)=Feed_total(1,j)*(1-Yf(1,j)); %starting sus Ent in feed: CFU
Water_res(1,j)=Water_total(1,j)*Yw(1,j); %starting res Ent in water: CFU
Water_sus(1,j)=Water_total(1,j)*(1-Yw(1,j)); %starting sus Ent in water: CFU
Pen_res(1,j)=Pen_total(1,j)*Yp(1,j); %starting res Ent in pen: CFU
Pen_sus(1,j)=Pen_total(1,j)*(1-Yp(1,j)); %starting sus Ent in pen: CFU

Cow_total_conc(1,j)=Zc(1,j)*K_c(1,j); %starting concentration Ent in cattle: CFU/g-Ingesta
Feed_total_conc(1,j)=Zf(1,j)*K_f(1,j); %starting concentration Ent in feed: CFU/g-Feed
Water_total_conc(1,j)=Zw(1,j)*K_w(1,j); %starting concentration Ent in water: CFU/g-Water
Pen_total_conc(1,j)=Zp(1,j)*K_p(1,j); %starting concentration Ent in pen: CFU/g-Pen

i=2;

while i<((time/dt)+1)


%PK dosage information for single animal
if T(i,1)>TYL_start_time && T(i,1)<TYL_stop_time && Method==1
    TYL_feed(i,1)=90/24; %mg/hr
else
    TYL_feed(i,1)=0;
end

%PK Flow Parameters for single animal
stomach_deg=(delta(1,j))*TYL_stomach_amt(i-1,j); %amount of TYL degrading in the stomach: mg/hr
To_USI=lambda_s(1,j)*TYL_stomach_amt(i-1,j); %movement from stomach to usi: mg/hr

USI_deg=(delta(1,j))*TYL_USI_amt(i-1,j); %amount of TYL degraded in the USI: mg/hr
To_LSI=lambda_usi(1,j)*TYL_USI_amt(i-1,j); %movement of TYL from USI to LSI: mg/hr

LSI_deg=(delta(1,j))*TYL_LSI_amt(i-1,j); %amount of TYL degraded in the LSI: mg/hr
To_LI=lambda_lsi(1,j)*TYL_LSI_amt(i-1,j); %movement of TYL from LSI to LI: mg/hr

LI_deg=(delta(1,j))*TYL_LI_amt(i-1,j); %amount of TYL degraded in the LI: mg/hr
To_pen=lambda_li(1,j)*TYL_LI_amt(i-1,j); %movement of TYL from LI to pen: mg/hr


%PK compartment equations for single animal. must be positive
TYL_stomach_amt(i,j)=max(0,TYL_stomach_amt(i-1,j)+dt*(TYL_feed(i,1)-stomach_deg-To_USI)); 
    %calculate TYL in stomach at each time step: mg

TYL_USI_amt(i,j)=max(0,TYL_USI_amt(i-1,j)+dt*(To_USI-USI_deg-To_LSI)); 
    %calculate TYL in USI at each time step: mg

TYL_LSI_amt(i,j)=max(0,TYL_LSI_amt(i-1,j)+dt*(To_LSI-LSI_deg-To_LI));
    %calculate TYL in LSI at each time step: mg
 
TYL_LI_amt(i,j)=max(0,TYL_LI_amt(i-1,j)+dt*(To_LI-LI_deg-To_pen));
    %calculate TYL in LI at each time step: mg


%PK storage in arrays
TYL_li_conc(i,j)=(TYL_LI_amt(i,j)/V_li(i,1))*free(1,j);
    %store conc of TYL in LI in array: mg/L == ug/mL
    
%PD equations. assumes TYL concentration in all cattle is equal to that found for single animal above
Es(i,j)=E_0-((Emax*(TYL_li_conc(i-1,j).^Hill_sus(1,j)))/...
    ((MIC_sus(1,j).^Hill_sus(1,j))+(TYL_li_conc(i-1,j).^Hill_sus(1,j))));
    %calculated effect of TYL on sus Ent at each time step: unitless
    
Er(i,j)=E_0-((Emax*(TYL_li_conc(i-1,j).^Hill_res(1,j)))/...
    ((MIC_res(1,j).^Hill_res(1,j))+(TYL_li_conc(i-1,j).^Hill_res(1,j))));
    %calculated effect of TYL on res Ent at each time step: unitless

    
%Pop equations. Broken into chunks: Plasmid transfer, movement between compartments, growth/death
%plasmid transfer, only occurs if there are Ent
if Cow_total(i-1,j)>0 %no transfer in cattle if there are no Ent
    C_plasmid=plasmid(1,j)*Cow_res(i-1,j)*(Cow_sus(i-1,j)/Cow_total(i-1,j)); %CFU/hr
else
    C_plasmid=0;
end

if Feed_total(i-1,j)>0 %no transfer in feed if there are no Ent
    F_plasmid=plasmid(1,j)*Feed_res(i-1,j)*(Feed_sus(i-1,j)/Feed_total(i-1,j)); %CFU/hr
else
    F_plasmid=0;
end

if Water_total(i-1,j)>0 %no transfer in water if there are no Ent
    W_plasmid=plasmid(1,j)*Water_res(i-1,j)*(Water_sus(i-1,j)/Water_total(i-1,j)); %CFU/hr
else
    W_plasmid=0;
end

if Pen_total(i-1,j)>0 %no transfer in environ if there are no Ent
    P_plasmid=plasmid(1,j)*Pen_res(i-1,j)*(Pen_sus(i-1,j)/Pen_total(i-1,j)); %CFU/hr
else
    P_plasmid=0;
end


%movement values. must be positive
%Feed, pen environ, and water compartments each have 11 cattle drinking/eating from them
Sus_F_to_C=max(0,(Feed_sus(i-1,j)*F_to_C(1,j))); %sus Ent from feed to cattle (ingestion): CFU/hr
Res_F_to_C=max(0,(Feed_res(i-1,j)*F_to_C(1,j))); %res Ent from feed to cattle (ingestion): CFU/hr

Sus_F_to_P=max(0,Feed_sus(i-1,j)*F_to_P(1,j)); %sus Ent from feed to pen: CFU/hr
Res_F_to_P=max(0,Feed_res(i-1,j)*F_to_P(1,j)); %res Ent from feed to pen: CFU/hr

Sus_P_to_F=max(0,Pen_sus(i-1,j)*P_to_F(1,j)); %sus Ent from pen to feed: CFU/hr
Res_P_to_F=max(0,Pen_res(i-1,j)*P_to_F(1,j)); %res Ent from pen to feed: CFU/hr
    
Sus_W_to_C=max(0,(Water_sus(i-1,j)*W_to_C(1,j))); %sus Ent from water to cattle (ingestion): CFU/hr
Res_W_to_C=max(0,(Water_res(i-1,j)*W_to_C(1,j))); %res Ent from water to cattle (ingestion): CFU/hr

Sus_W_to_P=max(0,Water_sus(i-1,j)*W_to_P(1,j)); %sus Ent from water to pen: CFU/hr
Res_W_to_P=max(0,Water_res(i-1,j)*W_to_P(1,j)); %res Ent from water to pen: CFU/hr

Sus_P_to_W=max(0,Pen_sus(i-1,j)*P_to_W(1,j)); %sus Ent from pen to water: CFU/hr
Res_P_to_W=max(0,Pen_res(i-1,j)*P_to_W(1,j)); %res Ent from pen to water: CFU/hr

Sus_C_to_P=max(0,Cow_sus(i-1,j)*C_to_P(1,j)); %sus Ent from cattle to pen (defecation): CFU/hr
Res_C_to_P=max(0,Cow_res(i-1,j)*C_to_P(1,j)); %res Ent from cattle to pen (defecation): CFU/hr

Sus_P_to_C=max(0,(Pen_sus(i-1,j)*P_to_C(1,j))); %sus Ent from pen to cattle (ingestion): CFU/hr
Res_P_to_C=max(0,(Pen_res(i-1,j)*P_to_C(1,j))); %res Ent from pen to cattle (ingestion): CFU/hr


%calculate growth rates of res and sus Ent in cattle
Carry_Capacity_Distance=1-(Cow_total_conc(i-1,j)/K_c(1,j));
%calculate penalty for growth based on distance of total population from carrying capacity: unitless

C_sus_growth=Cow_sus(i-1,j)*R_c(1,j)*Carry_Capacity_Distance*Es(i,j); %CFU/hr
C_res_growth=Cow_res(i-1,j)*R_c(1,j)*Carry_Capacity_Distance*Er(i,j)*(1-Fit(1,j)); %CFU/hr

%calculate death rate in feed, water, and pen for sus and res Ent
F_res_death=Feed_res(i-1,j)*R_f(1,j)*(1+Fit(1,j)); %CFU/hr
F_sus_death=Feed_sus(i-1,j)*R_f(1,j); %CFU/hr
W_res_death=Water_res(i-1,j)*R_w(1,j)*(1+Fit(1,j)); %CFU/hr
W_sus_death=Water_sus(i-1,j)*R_w(1,j); %CFU/hr
P_res_death=Pen_res(i-1,j)*R_p(1,j)*(1+Fit(1,j)); %CFU/hr
P_sus_death=Pen_sus(i-1,j)*R_p(1,j); %CFU/hr

%put the chunks together
%calculate sus Ent in cattle at each time step. must be positive
Cow_sus(i,j)=max(0,Cow_sus(i-1,j)+...
    dt*(Ingested_death(1,j)*(Sus_F_to_C+Sus_W_to_C+Sus_P_to_C)...
    +C_sus_growth-Sus_C_to_P-C_plasmid)); %CFU

%calculate res Ent in cattle at each time step
Cow_res(i,j)=max(0,Cow_res(i-1,j)+...
    dt*(Ingested_death(1,j)*(Res_F_to_C+Res_W_to_C+Res_P_to_C)...
    +C_plasmid+C_res_growth-Res_C_to_P)); %CFU


%calculate res and sus Ent in feed. 
%after burn-in, Feed is replaced every 24 hours (reset to starting point after burn-in)
%AFTP_time occurs at day 113, so feed will also be replaced at that time
if T(i,1)>burn && mod(T(i,1),24)==0
    Feed_sus(i,j)=Feed_sus(burn/dt,j); %CFU
    Feed_res(i,j)=Feed_res(burn/dt,j); %CFU
else
    Feed_res(i,j)=max(0,Feed_res(i-1,j)+...
        dt*(Res_P_to_F-Res_F_to_C-Res_F_to_P+F_plasmid-F_res_death)); %CFU
    
    Feed_sus(i,j)=max(0,Feed_sus(i-1,j)+...
        dt*(Sus_P_to_F-Sus_F_to_C-Sus_F_to_P-F_plasmid-F_sus_death)); %CFU
end

%res and sus Ent in water. Water trough is cleaned every two weeks (reset to starting point after burn-in)
%AFTP_time does not occur on an even two week point, but new clean water trough will be used at this time
if T(i,1)>burn && (mod(T(i,1),24*14)==0 || ((Treatment==2 || Treatment==4) && T(i,1)==AFTP_time))
    Water_sus(i,j)=Water_sus(burn/dt,j); %CFU
    Water_res(i,j)=Water_res(burn/dt,j); %CFU
else
   Water_sus(i,j)=max(0,Water_sus(i-1,j)+...
        dt*(Sus_P_to_W-Sus_W_to_C-Sus_W_to_P-W_plasmid-W_sus_death)); %CFU 
    
   Water_res(i,j)=max(0,Water_res(i-1,j)+...
        dt*(Res_P_to_W-Res_W_to_C-Res_W_to_P+W_plasmid-W_res_death)); %CFU
end

%calculate res and sus Ent in environment at each time stop    
%for AFTP (Treatment 2 and 4), eliminate resistant Ent in Pen at AFTP_time. Sus Ent reset to starting point after burn-in
%AFTP necessarily occurs after burn time
if (Treatment==2 || Treatment==4) && T(i,1)==AFTP_time
    Pen_res(i,j)=0; %CFU
    Pen_sus(i,j)=Pen_sus(burn/dt,j); %CFU
else
    Pen_sus(i,j)=max(0,Pen_sus(i-1,j)+...
        dt*(Sus_W_to_P+Sus_C_to_P+Sus_F_to_P-Sus_P_to_W-Sus_P_to_C-Sus_P_to_F...
        -P_plasmid-P_sus_death)); %CFU   
    
    Pen_res(i,j)=max(0,Pen_res(i-1,j)+...
        dt*(Res_W_to_P+Res_C_to_P+Res_F_to_P-Res_P_to_W-Res_P_to_C-Res_P_to_F...
        +P_plasmid-P_res_death)); %CFU
end

%calculate total Ent in each compartment at each time step
Cow_total(i,j)=Cow_sus(i,j)+Cow_res(i,j);
Feed_total(i,j)=Feed_sus(i,j)+Feed_res(i,j);
Pen_total(i,j)=Pen_sus(i,j)+Pen_res(i,j);
Water_total(i,j)=Water_sus(i,j)+Water_res(i,j);

%calculate concentration of Ent in each compartment
Cow_total_conc(i,j)=Cow_total(i,j)/(LI(i,1)*N_c(1,j));
Feed_total_conc(i,j)=Feed_total(i,j)/Bunk(1,j);
Pen_total_conc(i,j)=Pen_total(i,j)/Pen(1,j);
Water_total_conc(i,j)=Water_total(i,j)/Trough(1,j);

i=i+1;
end
end

%proportion resistant and susceptible
Prop_cow_res=Cow_res./Cow_total;
Prop_cow_sus=Cow_sus./Cow_total;
Prop_feed_res=Feed_res./Feed_total;
Prop_feed_sus=Feed_sus./Feed_total;
Prop_water_res=Water_res./Water_total;
Prop_water_sus=Water_sus./Water_total;
Prop_pen_res=Pen_res./Pen_total;
Prop_pen_sus=Pen_sus./Pen_total;

%bind parameter values together and export for sensitivity analysis
MC_parameters=transpose(vertcat(delta, lambda_s, lambda_usi, lambda_lsi, lambda_li, mu,...
	aero, log2_MIC_sus, log2_MIC_res, Hill_sus, Hill_res, W_c, F_c, P_c, DFM, log10_K_c, plasmid, R_c, Trough, log10_K_w,...
    R_w, Bunk, K_f, R_f, Ingested_death, Pen, log10_K_p, R_p, Fit, W_to_P, P_to_W, F_to_P, P_to_F, C_to_P, Zc, Zw, Zf,...
    Zp, Yc, Yw, Yf, Yp));

%Days for plotting
Days=T/24;

%export the proportion resistant, number of Ent, TYL_li_conc for analysis in R
%(where I can extract mean, median, etc)
%define base filename outside of script: e.g. 'TYL_NI_'
%round proportions of resistant to 4 decimal places, total Ent to 0 decimal places, Ent conc to 2 decimal places
%don't round MC parameters--some (e.g. plasmid) will be 0
%only export data for every hour after burn-in (every 10 time steps) to reduce file size
writematrix(round(Prop_cow_res(burn/dt+1:10:end,:),4), strcat(filename,'Prop_cow_res.txt'));
writematrix(round(Cow_total(burn/dt+1:10:end,:),0), strcat(filename,'Cow_total.txt'));
writematrix(round(Cow_total_conc(burn/dt+1:10:end,:),2), strcat(filename,'Cow_total_conc.txt'));

writematrix(round(Prop_feed_res(burn/dt+1:10:end,:),4), strcat(filename,'Prop_feed_res.txt'));
writematrix(round(Feed_total(burn/dt+1:10:end,:),0), strcat(filename,'Feed_total.txt'));
writematrix(round(Feed_total_conc(burn/dt+1:10:end,:),2), strcat(filename,'Feed_total_conc.txt'));

writematrix(round(Prop_water_res(burn/dt+1:10:end,:),4), strcat(filename,'Prop_water_res.txt'));
writematrix(round(Water_total(burn/dt+1:10:end,:),0), strcat(filename,'Water_total.txt'));
writematrix(round(Water_total_conc(burn/dt+1:10:end,:),2), strcat(filename,'Water_total_conc.txt'));

writematrix(round(Prop_pen_res(burn/dt+1:10:end,:),4), strcat(filename,'Prop_pen_res.txt'));
writematrix(round(Pen_total(burn/dt+1:10:end,:),0), strcat(filename,'Pen_total.txt'));
writematrix(round(Pen_total_conc(burn/dt+1:10:end,:),2), strcat(filename,'Pen_total_conc.txt'));

writematrix(round(TYL_li_conc(burn/dt+1:10:end,:),4), strcat(filename,'TYL_li_conc.txt'));
writematrix(MC_parameters, strcat(filename,'MC_parameters.txt'));
writematrix(Days(burn/dt+1:10:end,:), strcat(filename,'Days.txt'));
