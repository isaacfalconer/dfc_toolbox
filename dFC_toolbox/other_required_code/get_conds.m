function FC_project = get_conds(FC_project)

condition_choices = split(num2str([1:FC_project.ncon]),'  ');
condition_prompt = 'Which conditions would you like to include?';
tf = 0;
while tf==0
    if ~isfield(FC_project,'include')
        [conds,tf] = listdlg('PromptString',condition_prompt,'ListString',condition_choices);
        if tf==0
            FC_project.conds = [1:FC_project.ncon];
        else
            FC_project.conds = conds;
        end
    end
end

end