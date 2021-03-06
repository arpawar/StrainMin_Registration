%script from automatically generating compiled *.mex files from the matlab
%*.m files using the Matlab Coder toolbox (codegen command). 

clear all
close all
clc

%adding paths
addpaths();


%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
%cfg.ReportPotentialDifferences = false;
cfg.GlobalDataSyncMethod = 'NoSync';

cd ../iteration_loop_func
%% Define argument types for entry-point 'computenewPoints'.
ARGS = cell(1,1);
ARGS{1} = cell(7,1);
ARGS{1}{1} = struct;
ARGS{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{1} = coder.typeof(ARGS{1}{1},[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{3} = struct;
ARGS{1}{3}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{3} = coder.typeof(ARGS{1}{3},[Inf  1],[1 0]);
ARGS{1}{4} = struct;
ARGS{1}{4}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{4} = coder.typeof(ARGS{1}{4},[Inf  1],[1 0]);
ARGS{1}{5} = struct;
ARGS{1}{5}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{5} = coder.typeof(ARGS{1}{5},[Inf  1],[1 0]);
ARGS{1}{6} = struct;
ARGS{1}{6}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{6} = coder.typeof(ARGS{1}{6},[Inf  1],[1 0]);
ARGS{1}{7} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg computenewPoints -args ARGS{1}


%% Define argument types for entry-point 'compute_Integ_Domain_hyper_elastic'.
ARGS = cell(1,1);
ARGS{1} = cell(31,1);
ARGS{1}{1} = struct;
ARGS{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{1} = coder.typeof(ARGS{1}{1},[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{3} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{5} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{6} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{7} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{8} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{9} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{10} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{11} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{12} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{13} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{14} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{15} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{16} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{17} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{18} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{19} = struct;
ARGS{1}{19}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{19} = coder.typeof(ARGS{1}{19},[Inf  1],[1 0]);
ARGS{1}{20} = struct;
ARGS{1}{20}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{20} = coder.typeof(ARGS{1}{20},[Inf  1],[1 0]);
ARGS{1}{21} = struct;
ARGS{1}{21}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{21} = coder.typeof(ARGS{1}{21},[Inf  1],[1 0]);
ARGS{1}{22} = struct;
ARGS{1}{22}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{22} = coder.typeof(ARGS{1}{22},[Inf  1],[1 0]);
ARGS{1}{23} = coder.typeof(0);
ARGS{1}{24} = coder.typeof(0);
ARGS{1}{25} = coder.typeof(0);
ARGS{1}{26} = coder.typeof(0);
ARGS{1}{27} = coder.typeof(0,[4 1]);
ARGS{1}{28} = coder.typeof(0,[4 1]);
ARGS{1}{29} = coder.typeof(0,[4 1]);
ARGS{1}{30} = coder.typeof(single(0),[Inf  3],[1 0]);
ARGS{1}{31} = coder.typeof(0,[Inf  4  4],[1 0 0]);

%% Invoke MATLAB Coder.
codegen -config cfg compute_Integ_Domain_hyper_elastic -args ARGS{1}


%% Define argument types for entry-point 'compute_Integ_Domain_hyper_elastic_mod'.
ARGS = cell(1,1);
ARGS{1} = cell(32,1);
ARGS{1}{1} = struct;
ARGS{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{1} = coder.typeof(ARGS{1}{1},[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{3} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{5} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{6} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{7} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{8} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{9} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{10} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{11} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{12} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{13} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{14} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{15} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{16} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{17} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{18} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{19} = struct;
ARGS{1}{19}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{19} = coder.typeof(ARGS{1}{19},[Inf  1],[1 0]);
ARGS{1}{20} = struct;
ARGS{1}{20}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{20} = coder.typeof(ARGS{1}{20},[Inf  1],[1 0]);
ARGS{1}{21} = struct;
ARGS{1}{21}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{21} = coder.typeof(ARGS{1}{21},[Inf  1],[1 0]);
ARGS{1}{22} = struct;
ARGS{1}{22}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{22} = coder.typeof(ARGS{1}{22},[Inf  1],[1 0]);
ARGS{1}{23} = coder.typeof(0);
ARGS{1}{24} = coder.typeof(0);
ARGS{1}{25} = coder.typeof(0);
ARGS{1}{26} = coder.typeof(0);
ARGS{1}{27} = coder.typeof(0);
ARGS{1}{28} = coder.typeof(0,[4 1]);
ARGS{1}{29} = coder.typeof(0,[4 1]);
ARGS{1}{30} = coder.typeof(0,[4 1]);
ARGS{1}{31} = coder.typeof(single(0),[Inf  3],[1 0]);
ARGS{1}{32} = coder.typeof(0,[Inf  4  4],[1 0 0]);

%% Invoke MATLAB Coder.
codegen -config cfg compute_Integ_Domain_hyper_elastic_mod -args ARGS{1}

%% Define argument types for entry-point 'compute_Integ_Domain_hyper_elastic_theta'.
ARGS = cell(1,1);
ARGS{1} = cell(33,1);
ARGS{1}{1} = struct;
ARGS{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{1} = coder.typeof(ARGS{1}{1},[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{3} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{5} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{6} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{7} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{8} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{9} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{10} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{11} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{12} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{13} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{14} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{15} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{16} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{17} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{18} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{19} = struct;
ARGS{1}{19}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{19} = coder.typeof(ARGS{1}{19},[Inf  1],[1 0]);
ARGS{1}{20} = struct;
ARGS{1}{20}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{20} = coder.typeof(ARGS{1}{20},[Inf  1],[1 0]);
ARGS{1}{21} = struct;
ARGS{1}{21}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{21} = coder.typeof(ARGS{1}{21},[Inf  1],[1 0]);
ARGS{1}{22} = struct;
ARGS{1}{22}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{22} = coder.typeof(ARGS{1}{22},[Inf  1],[1 0]);
ARGS{1}{23} = coder.typeof(0);
ARGS{1}{24} = coder.typeof(0);
ARGS{1}{25} = coder.typeof(0);
ARGS{1}{26} = coder.typeof(0);
ARGS{1}{27} = coder.typeof(0);
ARGS{1}{28} = coder.typeof(0,[4 1]);
ARGS{1}{29} = coder.typeof(0,[4 1]);
ARGS{1}{30} = coder.typeof(0,[4 1]);
ARGS{1}{31} = coder.typeof(single(0),[Inf  3],[1 0]);
ARGS{1}{32} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{33} = coder.typeof(0,[Inf  4  4 4],[1 0 0]);
%% Invoke MATLAB Coder.
codegen -config cfg compute_Integ_Domain_hyper_elastic_theta -args ARGS{1}

%% Define argument types for entry-point 'compute_Integ_Domain_hyper_elastic_skin_theta'.
ARGS = cell(1,1);
ARGS{1} = cell(36,1);
ARGS{1}{1} = struct;
ARGS{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{1} = coder.typeof(ARGS{1}{1},[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{3} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{5} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{6} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{7} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{8} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{9} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{10} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{11} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{12} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{13} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{14} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{15} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{16} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{17} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{18} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{19} = struct;
ARGS{1}{19}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{19} = coder.typeof(ARGS{1}{19},[Inf  1],[1 0]);
ARGS{1}{20} = struct;
ARGS{1}{20}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{20} = coder.typeof(ARGS{1}{20},[Inf  1],[1 0]);
ARGS{1}{21} = struct;
ARGS{1}{21}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{21} = coder.typeof(ARGS{1}{21},[Inf  1],[1 0]);
ARGS{1}{22} = struct;
ARGS{1}{22}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{22} = coder.typeof(ARGS{1}{22},[Inf  1],[1 0]);
ARGS{1}{23} = coder.typeof(0);
ARGS{1}{24} = coder.typeof(0);
ARGS{1}{25} = coder.typeof(0);
ARGS{1}{26} = coder.typeof(0);
ARGS{1}{27} = coder.typeof(0);
ARGS{1}{28} = coder.typeof(0,[4 1]);
ARGS{1}{29} = coder.typeof(0,[4 1]);
ARGS{1}{30} = coder.typeof(0,[4 1]);
ARGS{1}{31} = coder.typeof(single(0),[Inf  3],[1 0]);
ARGS{1}{32} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{33} = coder.typeof(0,[Inf  4  4 4],[1 0 0]);
ARGS{1}{34} = coder.typeof(0,[Inf  4  4 4],[1 0 0]);
ARGS{1}{35} = coder.typeof(0,[Inf  4  4 4],[1 0 0]);
ARGS{1}{36} = coder.typeof(0,[Inf  4  4 4],[1 0 0]);
%% Invoke MATLAB Coder.
codegen -config cfg compute_Integ_Domain_hyper_elastic_skin_theta -args ARGS{1}

%% Define argument types for entry-point 'compute_Integ_Domain_hyper_elastic_growth_projection'.
ARGS = cell(1,1);
ARGS{1} = cell(39,1);
ARGS{1}{1} = struct;
ARGS{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{1} = coder.typeof(ARGS{1}{1},[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{3} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{5} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{6} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{7} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{8} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{9} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{10} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{11} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{12} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{13} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{14} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{15} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{16} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{17} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{18} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{19} = struct;
ARGS{1}{19}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{19} = coder.typeof(ARGS{1}{19},[Inf  1],[1 0]);
ARGS{1}{20} = struct;
ARGS{1}{20}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{20} = coder.typeof(ARGS{1}{20},[Inf  1],[1 0]);
ARGS{1}{21} = struct;
ARGS{1}{21}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{21} = coder.typeof(ARGS{1}{21},[Inf  1],[1 0]);
ARGS{1}{22} = struct;
ARGS{1}{22}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{22} = coder.typeof(ARGS{1}{22},[Inf  1],[1 0]);
ARGS{1}{23} = coder.typeof(0);
ARGS{1}{24} = coder.typeof(0);
ARGS{1}{25} = coder.typeof(0);
ARGS{1}{26} = coder.typeof(0);
ARGS{1}{27} = coder.typeof(0);
ARGS{1}{28} = coder.typeof(0,[4 1]);
ARGS{1}{29} = coder.typeof(0,[4 1]);
ARGS{1}{30} = coder.typeof(0,[4 1]);
ARGS{1}{31} = coder.typeof(single(0),[Inf  3],[1 0]);
ARGS{1}{32} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{33} = coder.typeof(0);
ARGS{1}{34} = coder.typeof(0);
ARGS{1}{35} = coder.typeof(0,[Inf  4  4 4],[1 0 0]);
ARGS{1}{36} = coder.typeof(0);
ARGS{1}{37} = coder.typeof(0);
ARGS{1}{38} = coder.typeof(0,[Inf  4  4 4],[1 0 0]);
ARGS{1}{39} = coder.typeof(0,[Inf  4  4 4],[1 0 0]);
ARGS{1}{40} = coder.typeof(0,[Inf  4  4 4],[1 0 0]);
%% Invoke MATLAB Coder.
codegen -config cfg compute_Integ_Domain_skin_growth_projection -args ARGS{1}

%% Define argument types for entry-point 'compute_Integ_Domain_hyper_elastic_growth'.
ARGS = cell(1,1);
ARGS{1} = cell(34,1);
ARGS{1}{1} = struct;
ARGS{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{1} = coder.typeof(ARGS{1}{1},[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{3} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{5} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{6} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{7} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{8} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{9} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{10} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{11} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{12} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{13} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{14} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{15} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{16} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{17} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{18} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{19} = struct;
ARGS{1}{19}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{19} = coder.typeof(ARGS{1}{19},[Inf  1],[1 0]);
ARGS{1}{20} = struct;
ARGS{1}{20}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{20} = coder.typeof(ARGS{1}{20},[Inf  1],[1 0]);
ARGS{1}{21} = struct;
ARGS{1}{21}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{21} = coder.typeof(ARGS{1}{21},[Inf  1],[1 0]);
ARGS{1}{22} = struct;
ARGS{1}{22}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{22} = coder.typeof(ARGS{1}{22},[Inf  1],[1 0]);
ARGS{1}{23} = coder.typeof(0);
ARGS{1}{24} = coder.typeof(0);
ARGS{1}{25} = coder.typeof(0);
ARGS{1}{26} = coder.typeof(0);
ARGS{1}{27} = coder.typeof(0,[4 1]);
ARGS{1}{28} = coder.typeof(0,[4 1]);
ARGS{1}{29} = coder.typeof(0,[4 1]);
ARGS{1}{30} = coder.typeof(single(0),[Inf  3],[1 0]);
ARGS{1}{31} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{32} = coder.typeof(0);
ARGS{1}{33} = coder.typeof(0);
ARGS{1}{34} = coder.typeof(0,[Inf  4  4 4],[1 0 0]);
 
%% Invoke MATLAB Coder.
codegen -config cfg compute_Integ_Domain_hyper_elastic_growth -args ARGS{1}

%% Define argument types for entry-point 'tripleIterLoop'.
ARGS = cell(1,1);
ARGS{1} = cell(4,1);
ARGS{1}{1} = coder.typeof(0,[1 3]);
ARGS{1}{2} = struct;
ARGS{1}{2}.active_cell = coder.typeof(0);
ARGS{1}{2}.phi = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2}.phiu = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2}.phiv = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2}.phiw = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2},[Inf  1],[1 0]);
ARGS{1}{3} = struct;
ARGS{1}{3}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(ARGS{1}{3},[Inf  1],[1 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{5} = coder.typeof(0,[Inf Inf Inf]);
ARGS{1}{6} = coder.typeof(0);
ARGS{1}{7} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg tripleIterLoop -args ARGS{1}
%
%% Define argument types for entry-point 'extract_normal'.
ARGS = cell(1,1);
ARGS{1} = cell(4,1);
ARGS{1}{1} = coder.typeof(0,[Inf Inf Inf]);
ARGS{1}{2} = coder.typeof(0,[Inf Inf Inf]);
ARGS{1}{3} = coder.typeof(0,[Inf Inf Inf]);
ARGS{1}{4} = coder.typeof(0,[Inf Inf Inf]);
%% Invoke MATLAB Coder.
codegen -config cfg extract_normal -args ARGS{1}

%% Define argument types for entry-point 'tripleIterLoop_skin'.
ARGS = cell(1,1);
ARGS{1} = cell(10,1);
ARGS{1}{1} = coder.typeof(0,[1 3]);
ARGS{1}{2} = struct;
ARGS{1}{2}.active_cell = coder.typeof(0);
ARGS{1}{2}.phi = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2}.phiu = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2}.phiv = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2}.phiw = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2},[Inf  1],[1 0]);
ARGS{1}{3} = struct;
ARGS{1}{3}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(ARGS{1}{3},[Inf  1],[1 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{5} = coder.typeof(0,[Inf Inf Inf]);
ARGS{1}{6} = coder.typeof(0);
ARGS{1}{7} = coder.typeof(0);
ARGS{1}{8} = coder.typeof(0,[Inf Inf Inf]);
ARGS{1}{9} = coder.typeof(0,[Inf Inf Inf]);
ARGS{1}{10} = coder.typeof(0,[Inf Inf Inf]);

%% Invoke MATLAB Coder.
codegen -config cfg tripleIterLoop_skin -args ARGS{1}

% 
%% Define argument types for entry-point 'tripleIterLoop_skin_growth'.
ARGS = cell(1,1);
ARGS{1} = cell(11,1);
ARGS{1}{1} = coder.typeof(0,[1 3]);
ARGS{1}{2} = struct;
ARGS{1}{2}.active_cell = coder.typeof(0);
ARGS{1}{2}.phi = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2}.phiu = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2}.phiv = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2}.phiw = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2},[Inf  1],[1 0]);
ARGS{1}{3} = struct;
ARGS{1}{3}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(ARGS{1}{3},[Inf  1],[1 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{5} = coder.typeof(0,[Inf Inf Inf]);
ARGS{1}{6} = coder.typeof(0);
ARGS{1}{7} = coder.typeof(0);
ARGS{1}{8} = coder.typeof(0,[Inf Inf Inf]);
ARGS{1}{9} = coder.typeof(0,[Inf Inf Inf]);
ARGS{1}{10} = coder.typeof(0,[Inf Inf Inf]);
ARGS{1}{11} = coder.typeof(0,[Inf  1],[1 0]);
%% Invoke MATLAB Coder.
codegen -config cfg tripleIterLoop_skin_growth -args ARGS{1}

%% Define argument types for entry-point 'tripleIterLoop_growth2'.
ARGS = cell(1,1);
ARGS{1} = cell(8,1);
ARGS{1}{1} = coder.typeof(0,[1 3]);
ARGS{1}{2} = struct;
ARGS{1}{2}.active_cell = coder.typeof(0);
ARGS{1}{2}.phi = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2}.phiu = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2}.phiv = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2}.phiw = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2},[Inf  1],[1 0]);
ARGS{1}{3} = struct;
ARGS{1}{3}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(ARGS{1}{3},[Inf  1],[1 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{5} = coder.typeof(0,[Inf Inf Inf]);
ARGS{1}{6} = coder.typeof(0);
ARGS{1}{7} = coder.typeof(0);
ARGS{1}{8} = coder.typeof(0,[Inf  1],[1 0]);

%% Invoke MATLAB Coder.
codegen -config cfg tripleIterLoop_growth2 -args ARGS{1}

%% Define argument types for entry-point 'storePixelPhi'.
cd ../store_phi_functions

ARGS = cell(1,1);
ARGS{1} = cell(9,1);
ARGS{1}{1} = coder.typeof(int64(0));
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0,[Inf  3],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{4} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{5} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{6} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARGS{1}{7} = struct;
ARGS{1}{7}.knot_ind = coder.typeof(0,[Inf  3  2],[1 0 0]);
ARGS{1}{7}.flag_active = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.IEN = coder.typeof(0,[Inf  27],[1 0]);
ARGS{1}{7}.chdElem = coder.typeof(0,[Inf  8],[1 0]);
ARGS{1}{7}.cell_centre = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{7}.node = coder.typeof(0);
ARGS{1}{7}.parElem = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.actE = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nodes = coder.typeof(0,[Inf  8],[1 0]);
ARGS{1}{7} = coder.typeof(ARGS{1}{7},[Inf  1],[1 0]);
ARGS{1}{8} = struct;
ARGS{1}{8}.mat = coder.typeof(single(0),[Inf  27],[1 0]);
ARGS{1}{8} = coder.typeof(ARGS{1}{8},[Inf  1],[1 0]);
ARGS{1}{9} = struct;
ARGS{1}{9}.pU = coder.typeof(0);
ARGS{1}{9}.pV = coder.typeof(0);
ARGS{1}{9}.pW = coder.typeof(0);
ARGS{1}{9}.maxlevel = coder.typeof(0);
ARGS{1}{9}.nelemU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nelemV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nelemW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.orderGauss = coder.typeof(0);
ARGS{1}{9}.kU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.kV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.kW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nobU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nobV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nobW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.rho = coder.typeof(0,[3 1]);
ARGS{1}{9}.timestep = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.smallNumber = coder.typeof(0);
ARGS{1}{9}.mu = coder.typeof(0);
ARGS{1}{9}.lambda = coder.typeof(0);
ARGS{1}{9}.alpha = coder.typeof(0);
ARGS{1}{9}.beta = coder.typeof(0);
ARGS{1}{9} = coder.typeof(ARGS{1}{9});

%% Invoke MATLAB Coder.
codegen -config cfg storePixelPhi -args ARGS{1}

%% Define argument types for entry-point 'GaussPhi'.
ARGS = cell(1,1);
ARGS{1} = cell(7,1);
ARGS{1}{1} = coder.typeof(0,[Inf  2],[1 0]);
ARGS{1}{2} = struct;
ARGS{1}{2}.knot_ind = coder.typeof(0,[Inf  3  2],[1 0 0]);
ARGS{1}{2}.flag_active = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.IEN = coder.typeof(0,[Inf  27],[1 0]);
ARGS{1}{2}.chdElem = coder.typeof(0,[Inf  8],[1 0]);
ARGS{1}{2}.cell_centre = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{2}.node = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.parElem = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.actE = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.nodes = coder.typeof(0,[Inf  8],[1 0]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2},[Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{3} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{4} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{5} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARGS{1}{6} = struct;
ARGS{1}{6}.mat = coder.typeof(single(0),[Inf  27],[1 0]);
ARGS{1}{6} = coder.typeof(ARGS{1}{6},[Inf  1],[1 0]);
ARGS{1}{7} = struct;
ARGS{1}{7}.pU = coder.typeof(0);
ARGS{1}{7}.pV = coder.typeof(0);
ARGS{1}{7}.pW = coder.typeof(0);
ARGS{1}{7}.maxlevel = coder.typeof(0);
ARGS{1}{7}.nelemU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nelemV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nelemW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.orderGauss = coder.typeof(0);
ARGS{1}{7}.kU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.kV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.kW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nobU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nobV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nobW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.rho = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.timestep = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.smallNumber = coder.typeof(0);
ARGS{1}{7}.mu = coder.typeof(0);
ARGS{1}{7}.lambda = coder.typeof(0);
ARGS{1}{7}.alpha = coder.typeof(0);
ARGS{1}{7}.beta = coder.typeof(0);
ARGS{1}{7} = coder.typeof(ARGS{1}{7});

%% Invoke MATLAB Coder.
codegen -config cfg GaussPhi -args ARGS{1}

cd ../bspline-util

%% Define argument types for entry-point 'FindSpan'.
ARGS = cell(1,1);
ARGS{1} = cell(4,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0,[1 Inf],[0 1]);

%% Invoke MATLAB Coder.
codegen -config cfg FindSpan -args ARGS{1}

cd ../thb_refinement
%% Define argument types for entry-point 'computeCoeffMat'.
ARGS = cell(1,1);
ARGS{1} = cell(7,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0,[27  1]);
ARGS{1}{5} = coder.typeof(0,[27  1]);
ARGS{1}{6} = coder.typeof(0,[27 64  2]);
ARGS{1}{7} = coder.typeof(0,[1 27]);

%% Invoke MATLAB Coder.
codegen -config cfg computeCoeffMat -args ARGS{1}
cd ../run
