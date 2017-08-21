%WingHalevyModel  This M-file runs the model for sulfur isotope fraction-
% ation by sulfate-reducing bacteria from Wing and Halevy (2014).
% 
% The goal is to test whether observed effects of phosphate limitation of 
% SRB are adequately captured by the model.  This exercise was conducted as
% part of the research described in Zaarur et al. (in prep.).
%
% [FLAG] represents possible issues with the parameterization of the model.
%
% *** MODEL DESCRIPTION ***
%
% The reaction network that is modeled is the following:
%
%  0        1            2         3        4         5         6        7
% inf -> SO4-2out <-> SO4-2in <-> APS <-> SO3-2 <-> H2Sin <-> H2Sout -> eff
%                  A           B       C         D       
%
% Here, inf = influent, eff = effluent, in = inside the cell, out = outside
% the cell.
%
% Reactions are as follows (step, substrates, products)
%        S#1      +   S#2     +   S#3   =   P#1       +   P#2       +   P#3     +   P#4
% A    SO4(2-)out + n H+out             =   SO4(2-)in + n H+in       
% B    SO4(2-)in  +   ATP(4-) + 2 H+in  =   APS(2-)   +   PPi(3-)
% C    APS(2-)    +   MKred             =   SO3(2-)   +   MKox      +   AMP(2-) + 2 H+in
% D    SO3(2-)    + 3 MKred   + 2 H+in  =   H2S       + 3 MKox      + 3 H2O
% 
% Isotopically-selective reactions are labeled A thru D following Wing and 
% Halevy (2014).  Sulfur pools are numbered 0 thru 7.  
%
% Steady-state is assumed, such that 0 and 7 are constant sources/sinks of
% SO4-2 and H2S, respectively.  The processes 0->1 and 6->7 (hereafter, 01
% and 67 for short) proceed at rate equal to the net sulfate reduction rate
% and do not fractionate sulfur isotopes (like in a chemostat).
%
% *** NOTATION AND UNIT CONVENTIONS ***
% Base units are:
%  time             seconds     s
%                   days        day
%  temperature      Kelvin      K
%  energy           Joule       J 
%  amount           mole        mol
%  volume           liter       L
%  voltage          volt        V
%  cell density     cell        cell
%  isotope ratio    unitless    fraction (i.e., not *1000 permil)
%
% *** REFERENCES ***
%
% Kreke, B., and Cypionka, H. (1992) Protonmotive force in freshwater
%  sulfate-reducing bacteria, and its role in sulfate accumulation in
%  Desulfobulbus propionicus. Arch. Microbiol. 158, 183-187.
% 
% Wing, B. A., and Halevy, I. (2014) Intracellular metabolite levels shape 
%  sulfur isotope fractionation during microbial sulfate respiration. Proc. 
%  Natl. Acad. Sci. U.S.A. 115, 18116–18125.
%
% Zaarur, S., Wang, D. T., Ono, S., and Bosak, T. (in prep.) Effects of 
%  organic substrate, growth rate, and phosphate limitation on the frac-
%  tionation of sulfur isotopes during microbial reduction of sulfate.
%
% See also XX.

% Last modified 02 Nov 2016 // MATLAB R2012b (8.0.0.783) 64-bit Windows 7
% David T. Wang (dtw@mit.edu)

clear all;
% close all;

CRed1 = [0.80, 0.00, 0.00];     % Define custom colors
CBlu1 = [0.28, 0.57, 0.81];
CGrn1 = [0.00, 0.46, 0.37];
COrg1 = [1.00, 0.70, 0.28];
CGry1 = [0.81, 0.81, 0.77];
CBlk1 = [0.28, 0.24, 0.20];
CYlw1 = [1.00, 0.85, 0.00];

R = 8.314;      % ideal gas constant,   J (mol K)^-1
e = 1.602e-19;  % elementary charge,    C
Navo = 6.022e23;% Avogadro's number,    mol^-1
F = e*Navo;     % Faraday constant,     C mol^-1 [or J V^-1 mol^-1]

sperday = 86400;% seconds per day (60*60*24)

T_C = 30;       % K&C expt temperature, degC
T = 273.15+T_C; % K&C expt temperature, K

% Data from Table 2 of Kreke and Cypionka
KCTab2Data = zeros(7,4);
KCTab2Data(:,1) = [  2.5,     5,    10,    20,   100,  1000,  4000]';   % SO4 added (uM)
KCTab2Data(:,2) = [  529,  1047,  2014,  3498,  5325,  8550, 11850]';   % SO4 intracellular observed (uM)
KCTab2Data(:,3) = [ 1676,  4100, 10577, 33813, 410e5,   1e7,   1e7]';   % SO4 intracellular calculated (uM)
KCTab2Data(:,4) = [  214,   209,   202,   174,    53,     9,     3]';   % SO4 Accumulation as SO4in/SO4added


%% Calculate DeltaG_Ao as a function of SO4out
% 
% The following calculations are from [Eqns. S13-S15] in Wing and Halevy, 
% to calculate DeltaG_Ao as a function of outside SO4-2 concentration.  The
% calculations are based on data from Fig. 3 of Kreke and Cypionka. 
% 
% However, DeltaG_Ao should not vary with SO4-2 concentration, unless this 
% represents contributions of both membrane potential and concentration 
% gradient in the thermodynamic term.  
%
% [FLAG] The way they did this is potentially erroneous.  In particular, 
% the x-axis of Fig. 3 in K&C is not actually [SO4]out at equilibrium; it
% is [SO4]added, which is >10x [SO4]out at equilibrium.  This may lead to 
% errors of an order of magnitude or greater in X vs. SO4out.  

SO4out = logspace(-6,-1,100)';   
                % SO4 outside concentration, M

nProtonsUptake = -0.3200*log10(SO4out*1e6)+2.68;    
                % num H+ symported per SO4 vs. SO4out (M)   [K&C Fig. 3]
logAcc = -0.6650*log10(SO4out*1e6)+3.0541;
                % log(SO4in/SO4out) vs. SO4out (M)          [K&C Fig. 3]              
                
PMF = -132e-3;  % proton-motive force,  V
MembranePotential = nProtonsUptake*PMF/2 + logAcc*2.3*R*T/F/2;
                % membrane potential (\Delta\psi), V        [W&H Eqn. S13]
DeltaG_Ao = MembranePotential.*(nProtonsUptake-2)*R*T/(2.3*R*T/F);
                % Gibbs energy of reaction A (SO4 uptake) at standard state
                % (i.e., 1 M H+)                            [W&H Eqn. S14]
ProtonRatio = 10.^((PMF-MembranePotential)/(2.3*R*T/F));
                % H+in/H+out                                [W&H Eqn. S15]

%% K&C Fig. 3 - H+ per sulfate symported and SO4in/SO4out vs. SO4out (uM)

figure(53); clf;
[ha53, ha53l, ha53r] = plotyy(SO4out, nProtonsUptake, SO4out, 10.^logAcc, @semilogx, @loglog, 'LineWidth', 2);

legend([ha53l ha53r], 'H^+/SO_4^{2-} symported', 'sulfate accumulation', 'Location', 'SouthWest');
legend boxoff;

axes(ha53(1));  % left plot (nProtonsUptake)
    set(ha53l, 'LineStyle', '--', 'LineWidth', 1.5, 'Color', CBlu1);
    set(gca(), 'XColor', 'k', 'YColor', 'k');
    set(gca(), 'FontSize', 10, 'TickLength',2.5*get(gca(),'TickLength'));
    set(gca(), 'Position', get(gca(),'Position').*[1 1 0.93 0.95]+[0 0.05 0 0]);    % [left bottom width height]
    set(gca(), 'XTickLabel',[], 'XAxisLocation', 'top');
    ylim([0 3]);
    ylabel('H^+ per SO_4^{2-} symported', 'FontSize', 10);

axes(ha53(2));  % right plot (Accumulation)
    set(ha53r, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', CRed1);
    set(gca(), 'XColor', 'k', 'YColor', 'k');
    set(gca(), 'FontSize', 10, 'TickLength',2.5*get(gca(),'TickLength'));
    ylabel('[SO_4^{2-}]_{in}/[SO_4^{2-}]_{out}', 'FontSize', 10);

    xlabel('[SO_4^{2-}]_{out} / (mol L^{-1})', 'FontSize', 10);

% set(gcf(), 'PaperPositionMode', 'auto');
% print(gcf(), '-depsc2', '-loose', 'K&C Figure 3.eps');

% [FLAG] The Kreke and Cypionka paper has inconsistencies in the data
% displayed in Fig. 3 compared to Table 2.  The table shows accumulation
% ratio as SO4in/SO4added, but the figure shows more data points than in
% the table, and some plotted points do not match data in the table.

%% W&H Fig. S1 - DeltaG_Ao (kJ/mol) vs. logSO4out (M)

figure(21); clf;    
semilogx(SO4out,DeltaG_Ao/1000, '-', 'Color', CGrn1, 'LineWidth', 1.5);

ylabel('\Delta_rG\circ for sulfate uptake / (kJ mol^{-1})', 'FontSize', 10);
xlabel('[SO_4^{2-}]_{out} / (mol L^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 10, 'TickLength', 2.5*get(gca(), 'TickLength'));
grid on;

% set(gcf(), 'paperpositionmode','auto');
% print(gcf(), '-depsc2', '-loose', 'W&H Figure S1.eps');

%% Calculate intracellular metabolite concentrations as a function of csSRR

% In vitro enzyme kinetic parameters [W&H Table S2].
Ks=[%#1     #2      #3      % KM for substrates, mM
     0.01   0       0       % A
    10.0    0.10    0;      % B
     0.02   0.10    0;      % C
     0.05   0.02    0];     % D
Kp=[%#1     #2      #3      #4      % KM for products, mM
     0.01   0       0       0;      % A
     0.17   0.13    0       0;      % B
     0.40   0.10    0.30    0;      % C
     0.01   0.02    0       0];     % D
Vf=[ 3.98e-20; 3.24e-19; 3.49e-19; 4.28e-19 ];  % V+ (in vitro), mol/cell/s

Ks = Ks / 1e3;      % convert to mol/L
Kp = Kp / 1e3;  

DeltaGo = [ 0; 55.9e3; 5.4e3; 31.2e3 ]; % DeltaGo, J mol^-1 [W&H Table S1]

% Tunable parameters [W&H Table S4]
ATP = 2.6 / 1e3;
AMP = 0.3 / 1e3;
MK  = 0.6 / 1e3;    % total menaquinone
MKratio = 100;      % MKred/MKox
uvvF = @(cssrr) 0.9*cssrr + 13.9;   % u_{vivo-vitro} = V+_vivo / V+_vitro
                                    % here, cssrr is in fmol cell^-1 day^-1

% Model parameters or data from experiments [Zaarur et al.]
T_C = 25;           % expt temperature,     degC
T = 273.15+T_C;     % expt temperature,     K

SO4out = 20.0 / 1e3;
H2Sout =  1.0 / 1e3;

csSRRf = logspace(-1.5,+3.5,51);    % csSRR [1x41], in fmol cell^-1 day^-1 
uvv = uvvF(csSRRf);                 % calculate u_{vivo-vitro}
csSRR = csSRRf ./ sperday * 1e-15;  % convert csSRR to mol cell^-1 s^-1
Vfv = Vf * uvv;                     % V+ (in vivo)

% Calculate standard Gibbs energy & related values for sulfate uptake (A).
% [FLAG] Regressions from possible mis-interpreted K&C Fig. 3, see above.
nProtonsUptake = -0.3200*log10(SO4out*1e6)+2.68;    
logAcc = -0.6650*log10(SO4out*1e6)+3.0541;
MembranePotential = nProtonsUptake*PMF/2 + logAcc*2.3*R*T/F/2;
DeltaG_Ao = MembranePotential.*(nProtonsUptake-2)*R*T/(2.3*R*T/F);
ProtonRatio = 10.^((PMF-MembranePotential)/(2.3*R*T/F));

DeltaGo(1) = DeltaG_Ao;     % add DeltaG_Ao to DeltaGo mat (= [A;B;C;D]);

% Calculate intracellular concentrations following W&H Eqns. S22-S25
MKox  = MK/(MKratio+1);
MKred = MKox*MKratio;
H2Sin = H2Sout;
SO4in = (SO4out/Ks(1,1) - csSRR./Vfv(1,:).*(1 + SO4out/Ks(1,1)) ) ...
    ./ (csSRR./(Vfv(1,:)*Kp(1,1)) + (ProtonRatio).^nProtonsUptake .* exp(DeltaG_Ao/R/T)/Ks(1,1) );
SO3 = ( (csSRR./Vfv(4,:))*Ks(4,1)*(Ks(4,2)^3)*(1+H2Sin*MKox.^3/(Kp(4,1)*(Kp(4,2)^3))) + H2Sin*(MKox^3)*exp(DeltaGo(4)/R/T) ) ...
    ./ ( (MKred^3)*(1-csSRR./Vfv(4,:)) );
APS = ( (csSRR./Vfv(3,:))*Ks(3,1)*Ks(3,2).*(1+SO3.*MKox.*AMP./prod(Kp(3,1:3))) + SO3*MKox*AMP*exp(DeltaGo(3)/R/T) ) ...
    ./ (MKred .* (1 - csSRR./Vfv(3,:)));
PPi = ( SO4in.*ATP.*(1-csSRR./Vfv(2,:)) - (csSRR./Vfv(2,:))*Kp(2,1)*Kp(2,2) ) ...
    ./ ( APS .* ( (csSRR./Vfv(2,:)).*Ks(2,1)*Ks(2,2)/(Kp(2,1)*Kp(2,2)) + exp(DeltaGo(2)/R/T) ) );

% Calculate Gibbs energy of reactions according to W&H Eqns. S26-S29
DeltaG_A = DeltaGo(1) + R*T*log( (SO4in./SO4out).*(ProtonRatio.^nProtonsUptake) );
DeltaG_B = DeltaGo(2) + R*T*log( (APS.*PPi)./(SO4in.*ATP) );
DeltaG_C = DeltaGo(3) + R*T*log( (SO3.*AMP.*MKox)./(APS.*MKred) );
DeltaG_D = DeltaGo(4) + R*T*log( (H2Sin./SO3).*((MKox./MKred).^3) );

DeltaG = [DeltaG_A; DeltaG_B; DeltaG_C; DeltaG_D];  % [4x41], J mol^-1

%% Calculate isotopic fractionation

% Calculate reversibilities f_{p,r}, as a function of csSRR [W&H Eqn. S12]
rev = exp(DeltaG/R/T);                          % [4x41], unitless

% Sulfur isotope fractionation factors* (3x/32)             [W&H Table S3]
%  ------------------------------------------------------------------
%   Var     Description                             Rows    Columns
%  ------------------------------------------------------------------
%   aEq     Equilibrium fractionation**             [A:D]   [34, 33]
%   aKf     Forward kinetic fractionation**         [A:D]   [34, 33]
%   lamEq   Mass-dependent exponent, equilibrium    [A:D]   [33]
%   lamKf   Mass-dependent exponent, forward kin.   [A:D]   [33]
%  -------------------------------------------------------------------
%  * See Wing and Halevy (2014) for references.
%  **Note that alpha values here take their 'normal' forms, i.e.:
%     - Equilibrium fractionation factor (aEq) for A <-> B is (RBeq/RAeq).
%     - Forward kinetic fractionation factor (aKf) for A -> B is (RAB/RA).
%    In W&H, alphas were defined 'upside down' (i.e., reciprocals).

lamEq = [0.5150, 0.5150, 0.5167, 0.5147]';      % Equilibrium (minor)
lamKf = [0.5146, 0.5146, 0.5146, 0.5146]';      % Forward kinetic (minor)

aEq = ones(4,2);                                % Equilibrium 
aEq(:,1) = [1.0000, 1.0000, 1.0060, 1.0650]';   % 34S, (B/A)^-1
aEq(:,2) = aEq(:,1).^lamEq;                     % Minor isotope, (B/A)^-1
aEq = aEq.^-1;                                  % Flip to (B/A)

aKf = ones(4,2);                                % Forward kinetic
aKf(:,1) = [0.9970, 1.0000, 1.0220, 1.0250]';   % 34S (A/AB)
aKf(:,2) = aKf(:,1).^lamKf;                     % Minor isotope (A/AB)
aKf = aKf.^-1;                                  % Flip to (AB/A)

% Calculate 34/32 of pools 1-4 relative to that of pool 5 for each csSRR.
% The 34/32 ratio of pool 5 divided by the 34/32 ratio of pool 1 (R5/R1)
% = R(H2Sin)/R(SO4-2out).  Because diffusion of H2S through cell membranes
% is assumed to be rapid (Wing and Halevy (2014) SI p. 4), H2Sin = H2Sout.

% First, write system of linear equations for 34/32 fractionation @ each  
% step.  This can be derived starting from W&H Eqns. S3-S7 [see DTW's
% handwritten notes from 28.Oct.2016 for the system of equations].

Adiag = bsxfun(@times,aKf(:,1),eye(4));
Adiagplusonecol = bsxfun(@times, permute(bsxfun(@times,-rev,aKf(:,1)./aEq(:,1)), [1 3 2]), eye(4));
Adiagplusonecol = [zeros(size(rev,1),1,size(rev,2)), Adiagplusonecol];

B = permute(1-rev, [1 3 2]);
B = B - Adiagplusonecol(:,5,:);     % use a little trick to get B(4,:,:)

A = bsxfun(@plus, Adiag, Adiagplusonecol(:,1:end-1,:));

% Now that we have written out the system of linear equations, we can solve
% it using some matrix math.  

X = zeros(size(B));                 % create matrix of solutions [4x1x41]

for i = 1:size(csSRR,2),            
    X(:,:,i) = A(:,:,i)\B(:,:,i);   % solve @ each csSRR, X = mldivide(A,B) 
                                    % X(:,1,i) = [R1; R2; R3; R4]./R5
end

X = squeeze(X);                     % [4x41] (drop extra 2nd dimension)

eps34 = 1000*(X(1,:)-1)             % enrichment factor (net), permil

%% W&H Fig. 3 - Model calibration to experimental data, for DMSS-1
figure(13); clf;
plot(csSRRf, eps34, '-o', 'Color', CRed1, 'LineWidth', 1.5);

ylabel(['^{34}\epsilon / ' char(8240) ''], 'FontSize', 10);
xlabel('csSRR / (fmol cell^{-1} day^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 10, 'TickLength', 2.5*get(gca(), 'TickLength'));
grid on;

% set(gcf(), 'paperpositionmode','auto');
% print(gcf(), '-depsc2', '-loose', 'W&H Figure 3.eps');

%% Figure X3
figure(3); clf;
semilogx(csSRRf, eps34, '-', 'Color', CBlk1, 'LineWidth', 1.5);

ylabel(['^{34}\epsilon / ' char(8240) ''], 'FontSize', 10);
xlabel('csSRR / (fmol cell^{-1} day^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 10, 'TickLength', 2.5*get(gca(), 'TickLength'));
grid on;

ee = 1000*(X-1);
% hold on;
% hee = plot(csSRRf, ee(2:end,:), '-');
% hold off;
% hleg = legend(hee, 'SO_4^{2-}_{in}', 'APS', 'SO_3^{2-}', 'Location', 'Best');

aa = csvread('forZaarurEtAl_ScratchBasicModels.csv',1)
hold on;
plot(aa(:,1), aa(:,2), '--', 'Color', CBlk1, 'LineWidth', 0.5);  % plot one-phase decay model
plot(aa(:,1), aa(:,3), '-', 'Color', CBlk1, 'LineWidth', 0.5);   % plot MM model
hold off;

ss = csvread('forZaarurEtAl_eps34vscsSRR_DMSS-1.csv',1,1)
ss(ss==0) = NaN;            % filter out zeros
hold on;
plot(ss(:,1), ss(:,3), 'o', 'MarkerEdgeColor', CRed1, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.0);
hold off;

set(gcf(), 'paperpositionmode','auto');
print(gcf(), '-depsc2', '-loose', 'Figure X3.eps');

% grid off

%% Figure X4

% Calculate Pi and ADP @ eq based on Thauer++ p. 101
Pieq = (PPi.*exp(-(-21.92e3)/R/T)).^(1/2);
ADPeq = (ATP.*AMP.*exp(-(0.00e3)/R/T)).^(1/2);

figure(4); clf;
set(gcf(), 'Position', [348, 49, 694, 652]);
hts = tight_subplot(3,2, [0.03 0.12], [0.10 0.05], [0.15 0.05])

axes(hts(1))
loglog(csSRRf, SO4in, 'o-', 'Color', CRed1, 'LineWidth', 1.5);

ylabel('[SO_4^{2-}]_{in} / (mol L^{-1})', 'FontSize', 10);
% xlabel('csSRR / (fmol cell^{-1} day^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 10, 'TickLength', 2.5*get(gca(), 'TickLength'));
grid on;
set(gca(),'XTickLabel',[])

axes(hts(3))
loglog(csSRRf, APS, 'o-', 'Color', CBlu1, 'LineWidth', 1.5);

ylabel('[APS] / (mol L^{-1})', 'FontSize', 10);
% xlabel('csSRR / (fmol cell^{-1} day^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 10, 'TickLength', 2.5*get(gca(), 'TickLength'));
grid on;
set(gca(),'XTickLabel',[])

axes(hts(5))
loglog(csSRRf, SO3, 'o-', 'Color', CGrn1, 'LineWidth', 1.5);

ylabel('[SO_3^{2-}] / (mol L^{-1})', 'FontSize', 10);
xlabel('csSRR / (fmol cell^{-1} day^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 10, 'TickLength', 2.5*get(gca(), 'TickLength'));
grid on;
ylim(10.^[-5 -2]);

axes(hts(2))
loglog(csSRRf, Pieq, 'o-', 'Color', CBlk1, 'LineWidth', 1.5);
% hold on; loglog(csSRRf, PPi, 'o-', 'Color', COrg1, 'LineWidth', 1.5); hold off;

ylabel('[Pi]_{in} / (mol L^{-1})', 'FontSize', 10);
% xlabel('csSRR / (fmol cell^{-1} day^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 10, 'TickLength', 2.5*get(gca(), 'TickLength'));
grid on;
set(gca(),'XTickLabel',[])

axes(hts(4))
loglog(csSRRf, PPi, 'o-', 'Color', COrg1, 'LineWidth', 1.5);

ylabel('[PPi] / (mol L^{-1})', 'FontSize', 10);
xlabel('csSRR / (fmol cell^{-1} day^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 10, 'TickLength', 2.5*get(gca(), 'TickLength'));
grid on;
% set(gca(),'XTickLabel',[])

axes(hts(6))
% set(gca(), 'Color','None',...
%            'XColor','k','YColor','k');
axis off;

set(gcf(), 'paperpositionmode','auto');
print(gcf(), '-depsc2', '-loose', 'Figure X4.eps');

%% Figure X5

figure(5); clf; 

loglog(Pieq, eps34, 'o-', 'Color', CBlk1, 'LineWidth', 1.5);
% hold on; loglog(csSRRf, PPi, 'o-', 'Color', COrg1, 'LineWidth', 1.5); hold off;

ylabel(['^{34}\epsilon / ' char(8240) ''], 'FontSize', 10);
xlabel('[Pi]_{in} / (mol L^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 10, 'TickLength', 2.5*get(gca(), 'TickLength'));
grid on;


%% Figure 7 - Intracellular Phosphorus vs. csSRR

% Total metabolite (transferable) P (i.e., not counting P in DNA or lipids)
TotPin = Pieq + PPi + APS + ATP + ADPeq + AMP;  

figure(7); clf;
hts = tight_subplot(2,1, [0.03 0.15], [0.05 0.05], [0.20 0.05])

axes(hts(1));
loglog(csSRRf, TotPin, '-', 'Color', [0.6700, 0.3100, 0.3200], 'LineWidth', 1.5);

ylabel('\SigmaP_{in} (mol L^{-1})', 'FontSize', 10);
% xlabel('csSRR / (fmol cell^{-1} day^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 9, 'TickLength', 4*get(gca(), 'TickLength'));

xlim(10.^[-1.9,+3.9])
ylim([3e-4 3e0])

set(gca(),'YTick',10.^[-3:1:+0])
set(gca(),'YMinorTick','on')

set(gca(),'XTick',10.^[-1:1:+3])
set(gca(),'XTickLabel',[])
grid on;
grid minor

% ylim([0 0.8])
% set(gca(),'YTick',[0:0.2:0.8])

% Panel Label
xl = xlim
yl = ylim
text(10.^(log10(xl(2)/xl(1))*0.07+log10(xl(1))),10.^(log10(yl(2)/yl(1))*0.90+log10(yl(1))),'a', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

axpos = get(gca(),'Position');
set(gca(),'Position', [axpos(1), axpos(2)+axpos(4)*(0.5), axpos(3), axpos(4)-axpos(4)*(0.5)])

axes(hts(2));
hPhos(1) = loglog(csSRRf, Pieq, '-', 'Color', CBlk1, 'LineWidth', 1.5); hold on
hPhos(2) = loglog(csSRRf, PPi, '-', 'Color', COrg1, 'LineWidth', 1.5); hold on
hPhos(3) = loglog(csSRRf, APS, '-', 'Color', CBlu1, 'LineWidth', 1.5); hold on
hPhos(4) = loglog(csSRRf, repmat(ATP, size(csSRRf)), '-', 'Color', CRed1, 'LineWidth', 1.5); hold on
hPhos(5) = loglog(csSRRf, repmat(ADPeq, size(csSRRf)), '-', 'Color', CGry1, 'LineWidth', 1.5); hold on
hPhos(6) = loglog(csSRRf, repmat(AMP, size(csSRRf)), '-', 'Color', CGrn1, 'LineWidth', 1.5); hold off

ylabel('intracellular concentration (mol L^{-1})', 'FontSize', 10);
xlabel('csSRR / (fmol cell^{-1} day^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 9, 'TickLength', 4*get(gca(), 'TickLength'));
grid on;
grid minor

ylim([3e-9 2e0])
set(gca(),'YTick',10.^[-9:1:+0])
% set(gca(),'YTickLabel',10.^[-9:1:+0])
set(gca(),'YMinorTick','on')

xlim(10.^[-1.9,+3.9])
% set(gca(),'XMinorTick','off')
set(gca(),'XTick',10.^[-1:1:+3])
set(gca(),'XTickLabel',10.^[-1:1:+3])

axpos = get(gca(),'Position');
set(gca(),'Position', [axpos(1), axpos(2), axpos(3), axpos(4)+axpos(4)*(0.5)])

legend(hPhos, {'Pi(eq)', 'PPi', 'APS', 'ATP', 'ADP(eq)', 'AMP'}, ...
    'Location', 'SouthOutside', ...
    'FontSize', 8);
% legend boxoff;

% axis square

% Panel Label
xl = xlim
yl = ylim
text(10.^(log10(xl(2)/xl(1))*0.07+log10(xl(1))),10.^(log10(yl(2)/yl(1))*0.95+log10(yl(1))),'b', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

% set(gcf(), 'Position', [305   297   579   401]);
set(gcf, 'position',  [  478   217   306   581])

set(gcf(), 'paperpositionmode','auto');
print(gcf(), '-depsc2', '-loose', 'Figure 7.eps');



%% Figure A3 - Reversibility & d34S vs. csSRR

figure(83); clf;
set(gcf(),'Position',[ 778   174   357   537])
hts = tight_subplot(2,1, [0.03 0.15], [0.1 0.05], [0.20 0.05])

% panel A, Reversibility vs csSRR
axes(hts(1));

thecm = [38  53  79    
        21 108 213
%         230 223 207
        217  49  43
        140  18  14]./255;
% colormap(gcf(), thecm);

for i = 1:size(rev,1),
    semilogx(csSRRf, rev(i,:).*100, 'color', thecm(i,:), 'linewidth', 1.25); hold on;
end
hold off;

hleg = legend('A','B','C','D', 'location', 'SW')
legend boxoff
set(hleg, 'fontsize', 8)

ylabel('reversibility (%)', 'FontSize', 10);
% xlabel('csSRR / (fmol cell^{-1} day^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 10, 'TickLength', 3*get(gca(), 'TickLength'));

xlim(10.^[-1.9,+3.9])
% ylim([-0.05 1.05])

% set(gca(),'YTick',10.^[-3:1:+0])
% set(gca(),'YMinorTick','on')

set(gca(),'XTick',10.^[-1:1:+3])
set(gca(),'XTickLabel',[])

grid on
% grid minor
set(gca, 'XMinorGrid', 'off')

% Panel Label
xl = xlim
yl = ylim
text(10.^(log10(xl(2)/xl(1))*0.07+log10(xl(1))),(yl(2)-yl(1))*0.88+yl(1),'a', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

axpos = get(gca(),'Position');
set(gca(),'Position', [axpos(1), axpos(2)+axpos(4)*(0.3), axpos(3), axpos(4)-axpos(4)*(0.3)])


% panel B, Plot d34S intracellular
axes(hts(2))

% Calculate intracellular d34S values w.r.t. SO4out
RvsSO4out = [bsxfun(@rdivide,X(:,:),X(1,:)); 1./X(1,:)]     % 34/32 of pool 1-5 vs. SO4out (pool 1)
dvsSO4out = RvsSO4out - 1

% Plot
% thecm =[228,26,28
%         55,126,184
%         77,175,74
%         152,78,163
%         255,127,0]./255;

thecm = [ 0 0 0; thecm];
% thecm = thecm(:,[3 2 1]);

for i = 1:size(dvsSO4out,1),
    if i==1, LS='--', else LS='.-', end
    semilogx(csSRRf, dvsSO4out(i,:).*1000, LS, 'color', thecm(i,:), 'linewidth', 1); hold on;
end
hold off;


hleg = legend('SO_4^{2-}(out)','SO_4^{2-}(in)','APS','SO_3^{2-}', 'H_2S', 'location', 'SE')
legend boxoff
set(hleg,'FontSize',8)

ylabel(['\delta^{34}S_{vs. SO_4^{2-} (out)} (' char(8240) ')'], 'FontSize', 10);
xlabel('csSRR / (fmol cell^{-1} day^{-1})', 'FontSize', 10);
set(gca(), 'FontSize', 10, 'TickLength', 3*get(gca(), 'TickLength'));

xlim(10.^[-1.9,+3.9])
ylim([-70 70])

% set(gca(),'YTick',10.^[-3:1:+0])
set(gca(),'YMinorTick','on')

set(gca(),'XTick',10.^[-1:1:+3])
set(gca(),'XTickLabel',10.^[-1:1:+3])

grid on
% grid minor
set(gca, 'XMinorGrid', 'off')

% Panel Label
xl = xlim
yl = ylim
text(10.^(log10(xl(2)/xl(1))*0.07+log10(xl(1))),(yl(2)-yl(1))*0.94+yl(1),'b', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

axpos = get(gca(),'Position');
set(gca(),'Position', [axpos(1), axpos(2), axpos(3), axpos(4)+axpos(4)*(0.3)])


set(gcf(), 'paperpositionmode','auto');
print(gcf(), '-depsc2', '-loose', 'Figure A3.eps');

