function FC_project = get_ncon_nsub(FC_project)

disp('Getting number of conditions...')

cd(FC_project.CONN_preprocessing_dir)

% Get number of conditions
i = 1;
j = 1;
while i==1
    if numel(dir(sprintf('*_Condition%03d.mat',j))) < 1
        ncon = j-1; % number of conditions in CONN output folder
        i = 0;
    end
    j = j+1;
end

disp(['Number of conditions = ',num2str(ncon)])
disp('Getting number of subjects...')

% Get number of subjects
i = 1;
j = 1;
while i==1
    if numel(dir(sprintf('ROI_Subject%03d_*',j))) < 1
        nsub = j-1;
        i = 0;
    end
    j = j+1;
end

disp(['Number of subjects = ',num2str(nsub)])

FC_project.ncon = ncon;
FC_project.nsub = nsub;

end