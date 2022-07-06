
% where is your eigen?
cg_eigen_dir = '/usr/local/include/eigen3/'; 


% output
outdir = 'solvers/';


%% Euclidean 3-point
% which solver to compile?
fname = 'solvers/solver_mini_estrigid_v2.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)

mex getdata_mini_estrigid.cpp 

%% Similarity 4-point
% which solver to compile?
fname = 'solvers/solver_mini_estsim.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)

