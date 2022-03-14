function [ Path_Used ] = Create_Path_Used()
%Create_Path_Used£ºSave the function and file path used in the function

Path_Used=struct();
Path_Used.map_path='..\data\Recon3D_301\map_ENSG_to_Entrez.xlsx';
Path_Used.cobra_refine_path='MATLAB_PATH\toolbox\cobratoolbox\src\reconstruction\refinement';
Path_Used.raven_core_path='MATLAB_PATH\toolbox\RAVEN\core';
Path_Used.raven_solver_path='MATLAB_PATH\toolbox\RAVEN\solver';
Path_Used.healthy_path='..\data\Healthy';
Path_Used.acknow_path='..\data\Recon3D_301\Acknow.mat';
end