function param = setparameters_uniaxial_local()
%In this function, set the required parameters for the input to the
%registration code

%degree of B-splines
pU = 2;
pV = 2;
pW = 2;

% maximum number of levels
maxlevel = 3;

% order of gauss quadrature
orderGauss = 4;

%initial grid: number of elements in x,y,z direction for level 1
m_var = 12;
n_var = 12;
o_var = 12;

%regularization parameters 
mu = 1.0;
lambda = 1.0;
alpha = 10.0;
beta = 1.0;

%gamma in g(x)term
smallNumber = 10^-12;

%number of elements in each direction at each level
nelemU = zeros(maxlevel,1);
nelemV = zeros(maxlevel,1);
nelemW = zeros(maxlevel,1);

%number of splines in each direction at each level
nobU = zeros(maxlevel,1);
nobV = zeros(maxlevel,1);
nobW = zeros(maxlevel,1);

%number of knot vectors in each direction at each level
kU = zeros(maxlevel,1); 
kV = zeros(maxlevel,1);
kW = zeros(maxlevel,1);

rho = zeros(maxlevel,1); %refinement parameter
timestep = zeros(maxlevel,1); %timestep for each refinement level

%number of elements in each direction
for level = 1:maxlevel
nelemU(level,1) = m_var*2^(level-1);
nelemV(level,1) = n_var*2^(level-1);
nelemW(level,1) = o_var*2^(level-1);

kU(level,1) = nelemU(level,1)+2*pU+1;
kV(level,1) = nelemV(level,1)+2*pU+1;
kW(level,1) = nelemW(level,1)+2*pU+1;

nobU(level,1) = kU(level,1) - pU - 1;
nobV(level,1) = kV(level,1) - pV - 1;
nobW(level,1) = kW(level,1) - pW - 1;
end

rho(1,1) = 1; %level 2 refinement
rho(2,1) = 3; %level 3 refinement
rho(3,1) = 6;

%lamda0--0.01
%lamda100--1e-6
%lambda2--0.005
%lambda2.5--0.003
%lambda5--0.0008
%lambda3--0.008
%lambda1--0.008

timestep(1,1) = 0.001;
timestep(2,1) = 0.001;
timestep(3,1) = 0.001;

%make a struct variable param, with all the parameters
param = struct('pU',pU,'pV',pV,'pW',pW,'maxlevel',maxlevel,'nelemU',nelemU,'nelemV',nelemV,'nelemW',nelemW,...
    'orderGauss',orderGauss,'kU',kU,'kV',kV,'kW',kW,'nobU',nobU,...
    'nobV',nobV,'nobW',nobW,'rho',rho,'timestep',timestep,'smallNumber',smallNumber,'mu',mu,'lambda',lambda,'alpha',alpha,'beta',beta);
end
