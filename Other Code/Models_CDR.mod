set A;
set F ordered;
set L{i in A};
set P := {i in A, j in A : i<j && card(L[i] inter L[j]) > 0};
param PI := 3.141592;
param eps := 1e-4;
param radius >= 0; 
param n >= 0;
param d >= 0;
param v0{A};
param cap{A};
param theta0{i in A} := if cap[i] >= PI then cap[i] - 2*PI else cap[i];
param x0{i in A} default -radius*cos((i-1)*2*PI/n + PI);
param y0{i in A} default -radius*sin((i-1)*2*PI/n + PI);
param xr0{(i,j) in P} default x0[i]-x0[j];
param yr0{(i,j) in P} default y0[i]-y0[j];
param vrx0{(i,j) in P} := v0[i]*cos(theta0[i]) - v0[j]*cos(theta0[j]);
param vry0{(i,j) in P} := v0[i]*sin(theta0[i]) - v0[j]*sin(theta0[j]);
param Rxr0{(i,j) in P} := -yr0[i,j];
param Ryr0{(i,j) in P} := xr0[i,j];
param dist0{(i,j) in P} := sqrt(xr0[i,j]**2 + yr0[i,j]**2);
param omega{(i,j) in P} = if abs(xr0[i,j]) > eps then atan(yr0[i,j]/xr0[i,j]) else if yr0[i,j] > 0 then -PI/2 else PI/2; 
param tm0{(i,j) in P} := -(xr0[i,j]*vrx0[i,j] + yr0[i,j]*vry0[i,j]);#(vrx0[i,j]**2 + vry0[i,j]**2);
param sep0{(i,j) in P} := (yr0[i,j]**2 - d**2)*vrx0[i,j]**2 + (xr0[i,j]**2 - d**2)*vry0[i,j]**2 - 2*vrx0[i,j]*vry0[i,j]*xr0[i,j]*yr0[i,j];
param tm{(i,j) in P};
param sep{(i,j) in P};
param nf;
param cf := 1e0;
param l0{A} default 1;
param nc;
param max_q;
param sum_q;
param sum_p;
param max_p;
param time_disjunctive;
param time_shadow;
param gain;

# speed rate bounds
param qmin := 0.9;
param qmax := 1.1;

# disjunctive parameters
param uvrx{(i,j) in P};
param uvry{(i,j) in P};
param lvrx{(i,j) in P};
param lvry{(i,j) in P};

param a{(i,j) in P} := yr0[i,j]**2 - d**2;
param b{(i,j) in P} := xr0[i,j]**2 - d**2;
param c{(i,j) in P} := 2*xr0[i,j]*yr0[i,j];

param Mbin11 {(i,j) in P};
param Mbin12 {(i,j) in P};
param Mbin311 {(i,j) in P : xr0[i,j] >=0 and yr0[i,j] < 0};
param Mbin312 {(i,j) in P : xr0[i,j] >=0 and yr0[i,j] < 0};
param Mbin321 {(i,j) in P : xr0[i,j] < 0 and yr0[i,j] >=0};
param Mbin322 {(i,j) in P : xr0[i,j] < 0 and yr0[i,j] >=0};
param Mbin331 {(i,j) in P : xr0[i,j] >=0 and yr0[i,j] >=0};
param Mbin332 {(i,j) in P : xr0[i,j] >=0 and yr0[i,j] >=0};
param Mbin341 {(i,j) in P : xr0[i,j] < 0 and yr0[i,j] < 0};
param Mbin342 {(i,j) in P : xr0[i,j] < 0 and yr0[i,j] < 0};

# shadow parameters
param l{(i,j) in P} = omega[i,j] + asin(d/dist0[i,j]); # l
param r{(i,j) in P} = omega[i,j] - asin(d/dist0[i,j]); # g
param Mvelocidade{(i,j) in P} := qmax*(v0[i] + v0[j]);
set R := {(i,j) in P: abs(xr0[i,j]) <= d};
set NR := P diff R;
param Mtanl{(i,j) in NR} := qmax*(v0[i] + v0[j])*(abs(tan(l[i,j])) + 1);
param Mtanlr{(i,j) in R} := qmax*(v0[i] + v0[j])*(abs(1/tan(l[i,j])) + 1);
param Mtanr{(i,j) in NR} := qmax*(v0[i] + v0[j])*(abs(tan(r[i,j])) + 1);
param Mtanrr{(i,j) in R} := qmax*(v0[i] + v0[j])*(abs(1/tan(r[i,j])) + 1);
param Mdvgn{(i,j) in NR} := (abs(xr0[i,j]) + abs(yr0[i,j]))*Mvelocidade[i,j];
param Mdvgnr{(i,j) in R} := (abs(Rxr0[i,j]) + abs(Ryr0[i,j]))*Mvelocidade[i,j];

# variables
var sigma{(i,j) in P,k in L[i] inter L[j]} binary;
var rho{(i,j) in P,k in L[i] inter L[j]} binary;
var delta{(i,j) in P,k in L[i] inter L[j]} binary;
var phi{(i,j) in P,k in L[i] inter L[j]} binary;
var tao{(i,j) in P,k in L[i] inter L[j]} binary;
var psi{(i,j) in P,k in L[i] inter L[j]} binary;
var dvrx{(i,j) in P} >= lvrx[i,j] <= uvrx[i,j];
var dvry{(i,j) in P} >= lvry[i,j] <= uvry[i,j];
var vrx{(i,j) in P};
var vry{(i,j) in P};
var q{i in A} >= qmin <= qmax default 1;
var Q{i in A} >=0 <= max(1-qmin,qmax-1);
var f{i in A,k in L[i]} binary;
var p{i in A} >= 0;
var fp{(i,j) in P} binary;
var z{(i,j) in P} binary;

# objective functions
minimize PLin: sum{i in A} (Q[i] + cf*p[i]);
minimize Quad: sum{i in A} ((1-q[i])**2 + cf*((p[i])**2));


####################### Disjunctive

s.t. cdvrx{(i,j) in P} : dvrx[i,j] = q[i]*v0[i]*cos(theta0[i]) - q[j]*v0[j]*cos(theta0[j]);
s.t. cdvry{(i,j) in P} : dvry[i,j] = q[i]*v0[i]*sin(theta0[i]) - q[j]*v0[j]*sin(theta0[j]);

s.t. bin11{(i,j) in P} : xr0[i,j]*dvry[i,j] - yr0[i,j]*dvrx[i,j] <= (1 - z[i,j]+ 1 - fp[i,j])*Mbin11[i,j];
s.t. bin12{(i,j) in P} : xr0[i,j]*dvry[i,j] - yr0[i,j]*dvrx[i,j] >=  -(1 + z[i,j]- fp[i,j])*Mbin12[i,j];
s.t. bin311{(i,j) in P : xr0[i,j] >=0 and yr0[i,j] < 0} : 2*a[i,j]*dvrx[i,j] - dvry[i,j]*(c[i,j] - sqrt(c[i,j]**2 - 4*a[i,j]*b[i,j])) <=  (2 - z[i,j]- fp[i,j])*Mbin311[i,j];
s.t. bin312{(i,j) in P : xr0[i,j] >=0 and yr0[i,j] < 0} : -2*b[i,j]*dvry[i,j] + dvrx[i,j]*(c[i,j] - sqrt(c[i,j]**2 - 4*a[i,j]*b[i,j])) <= (1 + z[i,j] - fp[i,j])*Mbin312[i,j];
s.t. bin321{(i,j) in P : xr0[i,j] < 0 and yr0[i,j] >=0} : -2*a[i,j]*dvrx[i,j] + dvry[i,j]*(c[i,j] - sqrt(c[i,j]**2 - 4*a[i,j]*b[i,j])) <= (2 - z[i,j]- fp[i,j])*Mbin321[i,j];
s.t. bin322{(i,j) in P : xr0[i,j] < 0 and yr0[i,j] >=0} : 2*b[i,j]*dvry[i,j] - dvrx[i,j]*(c[i,j] - sqrt(c[i,j]**2 - 4*a[i,j]*b[i,j])) <= (1 + z[i,j] - fp[i,j])*Mbin322[i,j];
s.t. bin331{(i,j) in P : xr0[i,j] >=0 and yr0[i,j] >=0} : 2*b[i,j]*dvry[i,j] - dvrx[i,j]*(c[i,j] + sqrt(c[i,j]**2 - 4*a[i,j]*b[i,j])) <=  (2 - z[i,j]- fp[i,j])*Mbin331[i,j];
s.t. bin332{(i,j) in P : xr0[i,j] >=0 and yr0[i,j] >=0} : 2*a[i,j]*dvrx[i,j] - dvry[i,j]*(c[i,j] + sqrt(c[i,j]**2 - 4*a[i,j]*b[i,j])) <= (1 + z[i,j] - fp[i,j])*Mbin332[i,j];
s.t. bin341{(i,j) in P : xr0[i,j] < 0 and yr0[i,j] < 0} : -2*b[i,j]*dvry[i,j] + dvrx[i,j]*(c[i,j] + sqrt(c[i,j]**2 - 4*a[i,j]*b[i,j])) <= (2 - z[i,j]- fp[i,j])*Mbin341[i,j];
s.t. bin342{(i,j) in P : xr0[i,j] < 0 and yr0[i,j] < 0} : -2*a[i,j]*dvrx[i,j] + dvry[i,j]*(c[i,j] + sqrt(c[i,j]**2 - 4*a[i,j]*b[i,j])) <= (1 + z[i,j] - fp[i,j])*Mbin342[i,j];

s.t. flight3{(i,j) in P, k in L[i] inter L[j]}: f[i,k] + f[j,k] <= 1 + fp[i,j];
s.t. flight4{i in A}: sum{k in L[i]} f[i,k] = 1;
s.t. flight6{i in A}: p[i] >= sum{k in L[i]}k*f[i,k] - l0[i];
s.t. flight7{i in A}: p[i] >= l0[i] - sum{k in L[i]}k*f[i,k];

# linear objective function constraint
s.t. cq1{i in A}: Q[i] >= 1- q[i];
s.t. cq2{i in A}: Q[i] >= q[i] - 1;


####################### Shadow

s.t. cvrx{(i,j) in NR}: vrx[i,j] = ((v0[i]*q[i])*(cos(theta0[i])) - ((v0[j]*q[j])*(cos(theta0[j]))));
s.t. cvry{(i,j) in NR}: vry[i,j] = ((v0[i]*q[i])*(sin(theta0[i])) - ((v0[j]*q[j])*(sin(theta0[j]))));
								
s.t. sigma1{(i,j) in NR,k in L[i] inter L[j]}: -vrx[i,j] <= Mvelocidade[i,j]*(1 - sigma[i,j,k]);						
s.t. sigma2{(i,j) in NR,k in L[i] inter L[j]}: vrx[i,j]*tan(l[i,j]) - vry[i,j] <= (1 - sigma[i,j,k])*Mtanl[i,j];
s.t. delta1{(i,j) in NR,k in L[i] inter L[j]}: -vrx[i,j] <= Mvelocidade[i,j]*(1 - delta[i,j,k]);
s.t. delta2{(i,j) in NR,k in L[i] inter L[j]}: -vrx[i,j]*tan(r[i,j]) + vry[i,j] <= Mtanr[i,j]*(1 - delta[i,j,k]);
s.t. phi1{(i,j) in NR,k in L[i] inter L[j]}: vrx[i,j] <= Mvelocidade[i,j]*(1 - phi[i,j,k]);						
s.t. phi2{(i,j) in NR,k in L[i] inter L[j]}: -vrx[i,j]*tan(l[i,j]) + vry[i,j] <= Mtanl[i,j]*(1 - phi[i,j,k]);
s.t. rho1{(i,j) in NR,k in L[i] inter L[j]}: vrx[i,j] <= Mvelocidade[i,j]*(1 - rho[i,j,k]);
s.t. rho2{(i,j) in NR,k in L[i] inter L[j]}: vrx[i,j]*tan(r[i,j]) - vry[i,j] <= Mtanr[i,j]*(1 - rho[i,j,k]);
s.t. dvg{(i,j) in NR,k in L[i] inter L[j]} : -(xr0[i,j]*vrx[i,j] + yr0[i,j]*vry[i,j]) <= (1 - tao[i,j,k])*Mdvgn[i,j];

s.t. rcvrx{(i,j) in R}: vrx[i,j] = ((v0[i]*q[i])*(cos(theta0[i] + PI/2)) - ((v0[j]*q[j])*(cos(theta0[j] + PI/2))));
s.t. rcvry{(i,j) in R}: vry[i,j] = ((v0[i]*q[i])*(sin(theta0[i] + PI/2)) - ((v0[j]*q[j])*(sin(theta0[j] + PI/2))));						

s.t. rsigma1{(i,j) in R,k in L[i] inter L[j]}: -vrx[i,j] <= Mvelocidade[i,j]*(1 - sigma[i,j,k]);						
s.t. rsigma2{(i,j) in R,k in L[i] inter L[j]}: vrx[i,j]*tan(l[i,j] + PI/2) - vry[i,j] <= (1 - sigma[i,j,k])*Mtanlr[i,j];
s.t. rdelta1{(i,j) in R,k in L[i] inter L[j]}: -vrx[i,j] <= Mvelocidade[i,j]*(1 - delta[i,j,k]);						
s.t. rdelta2{(i,j) in R,k in L[i] inter L[j]}: -vrx[i,j]*tan(r[i,j] + PI/2) + vry[i,j] <= Mtanrr[i,j]*(1 - delta[i,j,k]);
s.t. rphi1{(i,j) in R,k in L[i] inter L[j]}: vrx[i,j] <= Mvelocidade[i,j]*(1 - phi[i,j,k]);
s.t. rphi2{(i,j) in R,k in L[i] inter L[j]}: -vrx[i,j]*tan(l[i,j] + PI/2) + vry[i,j] <= Mtanlr[i,j]*(1 - phi[i,j,k]);
s.t. rrho1{(i,j) in R,k in L[i] inter L[j]}: vrx[i,j] <= Mvelocidade[i,j]*(1 - rho[i,j,k]);
s.t. rrho2{(i,j) in R,k in L[i] inter L[j]}: vrx[i,j]*tan(r[i,j] + PI/2) - vry[i,j] <= Mtanrr[i,j]*(1 - rho[i,j,k]);
s.t. rdvg{(i,j) in R,k in L[i] inter L[j]} : -(Rxr0[i,j]*vrx[i,j] + Ryr0[i,j]*vry[i,j]) <= (1 - tao[i,j,k])*Mdvgnr[i,j];

s.t. level1{i in A}: sum{k in L[i]} f[i,k] = 1;
s.t. level2{(i,j) in P,k in L[i] inter L[j]}: f[i,k] + f[j,k] + psi[i,j,k] <= 2;
s.t. level2b{(i,j) in P,k in L[i] inter L[j]}: f[i,k] + f[j,k] >=  psi[i,j,k] - 1;
s.t. level3{i in A}: p[i] >= sum{k in L[i]}k*f[i,k] - l0[i];
s.t. level4{i in A}: p[i] >= l0[i] - sum{k in L[i]}k*f[i,k];

s.t. r1{(i,j) in P,k in L[i] inter L[j]}: sigma[i,j,k] + delta[i,j,k] + rho[i,j,k] + phi[i,j,k] + psi[i,j,k] + tao[i,j,k] >= 1;
s.t. r2{(i,j) in P,k in L[i] inter L[j]}: sigma[i,j,k] + delta[i,j,k] + rho[i,j,k] + phi[i,j,k] + psi[i,j,k] - tao[i,j,k] <= 1;
