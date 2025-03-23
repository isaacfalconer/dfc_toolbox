if ~exist('dFC_project','var') || ~isfield(dFC_project,'dFC')
    msgbox('You must load or create a dFC project before running connectivity states analyis')
else
    maxk = inputdlg('What maximum value of k would you like to use? ', 'Choose max k value to try...', [1 45]);
    [condition,tf] = listdlg('SelectionMode','single','ListString',dFC_project.cond_names);
    if tf
        [dFC_project] = dFC_kmeans_app(dFC_project, condition, str2double(maxk{1}));
    end
end
save_path = which('findme');
old_path = pwd;
cd(save_path(1:length(save_path)-9))

if whos('dFC_project').bytes<2000000000
    save([dFC_project.project_name,'.mat'],'dFC_project')
else
    save([dFC_project.project_name,'.mat'], 'dFC_project','-v7.3')
end

cd(old_path)