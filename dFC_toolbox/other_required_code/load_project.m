[project_file,project_path] = uigetfile(pwd,'Select dFC project file',which('findme.m'));

app.DynamicFunctionalConnectivityToolboxUIFigure.Pointer = 'watch';

x = app.center_x;
y = app.center_y;
fig = uifigure('Position',[x-85 y-35 170 70]);
uitextarea(fig,'Value','Loading Project...','Position',[10 10 150 50]);
drawnow
figure(app.DynamicFunctionalConnectivityToolboxUIFigure)
figure(fig)
pause(0.01)

if ~(project_file(1) == 0)
    try
        dFC_project = load([project_path,project_file]).dFC_project_save;
        % dFC_project.dFC = dFC_project.dFC_pre_normalize;
    catch
        dFC_project = load([project_path,project_file]).dFC_project;
    end
end

close(fig)

app.DynamicFunctionalConnectivityToolboxUIFigure.Pointer = 'arrow';