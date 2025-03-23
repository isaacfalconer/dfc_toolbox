function FC_project = get_nroi(FC_project)

disp('Getting number of ROIs...')
cd(FC_project.CONN_results_dir)
FC_project.nroi = numel(load('resultsROI_Subject001_Condition001.mat').names);
disp(['Number of ROIs = ',num2str(FC_project.nroi)])

end