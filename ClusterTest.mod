set A;
set P := {i in A, j in A : i < j};
set C := {(i,j) in P};


param wd;
param n >= 0; #Number of Aircraft
param d >= 0; #Separation Distance (5nm)
param radius >= 0; #Radius of Circle
param level default 0; #Flight Level - Not used in 2D
param PI := 3.141592;
param eps := 1e-3;
param M := 1e4;
param v0{A}; #Initial Velocity of Aircraft
param cap{A}; #Initial Heading of Aircraft
param theta0{i in A} := if cap[i] >= PI then cap[i] - 2*PI else cap[i]; 
param x0{i in A} default -radius*cos((i-1)*2*PI/n + PI); #Initial x position
param y0{i in A} default -radius*sin((i-1)*2*PI/n + PI); #Initial y position
param xr0{i in A, j in A : i<j} := x0[i]-x0[j];  #Initial Relative x position (i,j)
param yr0{i in A, j in A : i<j} := y0[i]-y0[j]; #Initial Relative y position (i,j)

# control bounds
param qmin := 0.94; #Speed Adjustment - Lower Bound
param qmax := 1.03; #Speed Adjustment - Upper Bound
param hmin := -PI/6; #Heading Adjustment - Lower Bound
param hmax := +PI/6; #Heading Adjustment - Upper Bound

param ucc := 1;
param lcc := min(cos(hmin),cos(hmax));
param uss := sin(hmax);
param lss := sin(hmin);

# bigM values
param Mbin1 {(i,j) in P};
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

#bounds and coefficients
param uvrx{(i,j) in P}; #Relative Velocity x - Upper Bound
param uvry{(i,j) in P}; #Relative Velocity y - Upper Bound
param lvrx{(i,j) in P}; #Relative Velocity x - Lower Bound
param lvry{(i,j) in P}; #Relative Velocity y - Lower Bound
param a{(i,j) in P}; 
param b{(i,j) in P};
param c{(i,j) in P};
param gammal{(i,j) in P}; #Coefficient of relative velocity components for aircraft pair
param gammau{(i,j) in P}; #Coefficient of relative velocity components for aircraft pair
param phil{(i,j) in P}; #Coefficient of relative velocity components for aircraft pair
param phiu{(i,j) in P}; #Coefficient of relative velocity components for aircraft pair

# variables
var q{i in A} >= qmin <= qmax; #Speed Adjustment within Bounds
var theta{i in A} >= hmin <= hmax; #Heading Adjustment within Bounds
var vrx{(i,j) in P} >= lvrx[i,j] <= uvrx[i,j]; #Relative Velocity in x within Bounds
var vry{(i,j) in P} >= lvry[i,j] <= uvry[i,j]; #Relative Velocity in y within Bounds
var z{(i,j) in P} binary; #Disjunctive Binary Variable for Aircraft Pairs

#objective function
min 

#constraints
s.t. cvrx{(i,j) in P} : vrx[i,j] = q[i]*v0[i]*cos(theta[i] + theta0[i]) - q[j]*v0[j]*cos(theta[j]+theta0[j]);
s.t. cvry{(i,j) in P} : vry[i,j] = q[i]*v0[i]*sin(theta[i] + theta0[i]) - q[j]*v0[j]*sin(theta[j]+theta0[j]);
s.t. bin11{(i,j) in P} : xr0[i,j]*vry[i,j] - yr0[i,j]*vrx[i,j] <= (1 - z[i,j])*Mbin11[i,j];
s.t. bin12{(i,j) in P} : xr0[i,j]*vry[i,j] - yr0[i,j]*vrx[i,j] >=  -(z[i,j])*Mbin12[i,j];
s.t. bin311{(i,j) in P : xr0[i,j] >=0 and yr0[i,j] < 0} : gammal[i,j]*vrx[i,j] - phil[i,j]*vry[i,j] <= (1 - z[i,j])*Mbin311[i,j];
s.t. bin312{(i,j) in P : xr0[i,j] >=0 and yr0[i,j] < 0} : gammau[i,j]*vry[i,j] + phiu[i,j]*vrx[i,j] <= z[i,j]*Mbin312[i,j];
s.t. bin321{(i,j) in P : xr0[i,j] < 0 and yr0[i,j] >=0} : gammal[i,j]*vrx[i,j] + phil[i,j]*vry[i,j] <= (1 - z[i,j])*Mbin321[i,j];
s.t. bin322{(i,j) in P : xr0[i,j] < 0 and yr0[i,j] >=0} : gammau[i,j]*vry[i,j] - phiu[i,j]*vrx[i,j] <= z[i,j]*Mbin322[i,j];
s.t. bin331{(i,j) in P : xr0[i,j] >=0 and yr0[i,j] >=0} : gammal[i,j]*vry[i,j] - phil[i,j]*vrx[i,j] <= (1 - z[i,j])*Mbin331[i,j];
s.t. bin332{(i,j) in P : xr0[i,j] >=0 and yr0[i,j] >=0} : gammau[i,j]*vrx[i,j] - phiu[i,j]*vry[i,j] <= z[i,j]*Mbin332[i,j];
s.t. bin341{(i,j) in P : xr0[i,j] < 0 and yr0[i,j] < 0} : gammal[i,j]*vry[i,j] + phil[i,j]*vrx[i,j] <= (1 - z[i,j])*Mbin341[i,j];
s.t. bin342{(i,j) in P : xr0[i,j] < 0 and yr0[i,j] < 0} : gammau[i,j]*vrx[i,j] + phiu[i,j]*vry[i,j] <= z[i,j]*Mbin342[i,j];



