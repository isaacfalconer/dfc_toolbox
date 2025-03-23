
cla(app.UIAxes4)

threshold = app.ThresholdSpinner.Value;

node_labels = 0;

% figure('Position',[100,150,300,250])

% mygraph(mygraph<0) = 0;
A = binarize(mygraph,threshold);
mygraph = mygraph - min(mygraph,[],'all');
if ~exist('communities','var') || isempty(communities)
    app.UIAxes4.Title.String = 'Community Structure';
    rng_state = rng();
    rng(123)
    % [communities, ~] = community_louvain(mygraph);
    if state~=1
        if app.dFC_project.which_hemi~=1
            mygraph2 = app.dFC_project.state_dFC_mean{1}(ismember(app.dFC_project.include_hemi_sorted(hemi),rois),ismember(app.dFC_project.include_hemi_sorted(hemi),rois));
        else
            mygraph2 = app.dFC_project.state_dFC_mean{1}(ismember(app.dFC_project.include_hemi_sorted,rois),ismember(app.dFC_project.include_hemi_sorted,rois));
        end
        % mygraph2(mygraph2<0) = 0;
        % mygraph2 = mygraph2 - min(mygraph2,[],'all');
        % [communities2,~] = community_louvain(mygraph2);
        % communities = harmonize_comm(communities2,communities);
    end
    communities = comm;
    rng(rng_state)
else
    app.UIAxes4.Title.String = 'Networks';
end
ngroup = numel(unique(communities));
[~,index] = sort(communities);
mygraph = mygraph(index,index);
A = A(index,index);

numN = sqrt(numel(A));
% numN_group = numN/ngroup;
step = 2*pi()/ngroup; % Angle between groups
space = 2*pi()/(((sqrt(ngroup)*7)));
group_theta = 0:step:(2*pi()-step);
if ngroup==2
    group_theta = group_theta + pi()/4;
end

group_size = zeros(ngroup,1);
for i=1:ngroup
    group_size(i) = numel(communities(communities==i));
end

node_coor_group = cell(ngroup,1);
node_coor = zeros(numN,2);
label_coor_group = cell(ngroup,1);
label_coor = zeros(numN,2);
last = 0;

for i=1:ngroup
    node_coor_group{i} = zeros(group_size(i),2);
    numN_group = group_size(i);
    step_group = 2*pi()/numN_group; % Angle between nodes in group
    node_coor_i = 0:step_group:(2*pi()-step_group);
    group_x = cos(group_theta(i)); % x coordinate of group i
    group_y = sin(group_theta(i)); % y coordinate of group i
    node_x = (group_size(i)/max(group_size))*space*cos(node_coor_i - (((ngroup/2+1)-i)/ngroup)*2*pi()); % Relative x coordinates of group i nodes
    node_y = (group_size(i)/max(group_size))*space*sin(node_coor_i - (((ngroup/2+1)-i)/ngroup)*2*pi()); % Relative y coordinates of group i nodes
    node_coor_group{i}(:,1) = group_x + node_x; %(1:numN_group));
    node_coor_group{i}(:,2) = group_y + node_y; %(1:numN_group));
    label_coor_group{i}(:,1) = group_x + node_x + 0.15*cos(node_coor_i - (((ngroup/2+1)-i)/ngroup)*2*pi()); %(1:numN_group));
    label_coor_group{i}(:,2) = group_y + node_y + 0.15*sin(node_coor_i - (((ngroup/2+1)-i)/ngroup)*2*pi()); %(1:numN_group));
%     node_coor([1:numN_group]+last,:) = node_coor_group{i};
%     last = last + numN_group;
%     for j=1:numel(node_coor_group{i}(:,1))
%         if and(node_coor_group{i}(j,1) == 0, node_coor_group{i}(j,2) == 0)
%             disp(['Node ',num2str(j),' of Cluster ',num2str(i),' is at the origin'])
%         end
%     end
end

node_coor = cat(1,node_coor_group{:});
label_coor = cat(1,label_coor_group{:});

labels = cell(1,numN);
for i=1:numN
    labels{i} = num2str(i);
end

hold(app.UIAxes4,'on')

% colors = get_comm_colors();
if color==0
    colors = repmat([0 0 0],256,1);
else
    colors = hsv;
end
% colors = colors([1:34,64:256],:);

comm_colors = zeros(numel(communities),3);
for i=1:max(communities)
    % comm_colors(communities==i,:) = repmat(colors(1+ceil(size(colors,1)/ngroup)*(i-1),:),sum(communities==i),1);
    comm_colors(communities==i,:) = repmat(colors(i*floor(256/max(communities)),:),[sum(communities==i) 1]);
    % comm_colors(communities==(max(communities)+1-i),:) = repmat(colors(1+ceil(size(colors,1)/ngroup)*(i-1),:),sum(communities==(max(communities)+1-i)),1); % Reverse color order
end

comm_colors = color;

max_alpha = 1.5;
% max_alpha = max(mygraph,[],'all');


stop= 0;
for i = 1:(length(node_coor))
    for j=1:length(node_coor)
        if and(node_coor(i,1)==0, node_coor(i,2)==0) && stop==0
            disp(['i = ', num2str(i)])
            disp(['j = ', num2str(j)])
            stop = 1;
        end
        if i~=j && A(i,j)==1
            if mygraph(i,j) > max_alpha
                mygraph_ij = max_alpha;
            else
                mygraph_ij = mygraph(i,j);
            end
            patchline([node_coor(i,1) node_coor(j,1)],[node_coor(i,2) node_coor(j,2)],...
                app.UIAxes4,...
                'linewidth',1.5,...
                'edgealpha',(mygraph_ij-threshold)/(max_alpha-threshold),...
                'edgecolor',mean([comm_colors(index(j),:);comm_colors(index(i),:)],1))
        end
    end
end

comm_colors = comm_colors(index,:);

scatter(app.UIAxes4,node_coor(:,1),node_coor(:,2),5/(numN/(numN+90)),comm_colors,'filled')

if node_labels
    text(label_coor(:,1),label_coor(:,2),cellfun(@num2str,num2cell([1:numN]),'UniformOutput',false),'HorizontalAlignment','center')
end

if max(communities)>7
    disp('Warning: Network contains more than 7 communities, so colors were repeated')
end
% scatter(node_coor(:,1),node_coor(:,2),15,'MarkerFaceColor','k','MarkerEdgeColor','k')
% text(node_coor(:,1),node_coor(:,2),labels,'VerticalAlignment','top','HorizontalAlignment','left')

hold(app.UIAxes, 'off')

plotcenter = mean(node_coor,1);

if exist('names','var')
    names = names(index);
    % xlim(app.UIAxes4,[-(1+1.5*space),(2+1.5*space) + fix(numN/20)])
    % ylim(app.UIAxes4,[-(1+1.5*space),(1+1.5*space)])
    xlim(app.UIAxes4,[plotcenter(1)-1 plotcenter(1)+1])
    ylim(app.UIAxes4,[plotcenter(2)-1 plotcenter(2)+1])
    y_offset = 0;
    x_offset = 0;
    for i=1:numel(names)
        if y_offset > 2.5 && communities(index(i))~= communities(index(i-1))
            x_offset = 1;
            y_offset = 0;
        end
        if i>1 && communities(index(i))~= communities(index(i-1)) && y_offset ~= 0
            y_offset = y_offset + 0.1;
        end
        text(1.75 + x_offset ,1.5 - y_offset,names{i},'Color',comm_colors(i,:),'FontSize',6)
        y_offset = y_offset + 0.1;
    end
else
    % xlim(app.UIAxes4,[min(label_coor(:,1)),1.2*max(label_coor(:,1))])
    % ylim(app.UIAxes4,[1.1*min(label_coor(:,2)),1.1*max(label_coor(:,2))])
    xlim(app.UIAxes4,[min(node_coor(:,1))-0.3 max(node_coor(:,1))+0.3])
    ylim(app.UIAxes4,[min(node_coor(:,2))-0.3 max(node_coor(:,2))+0.3])
%     xlim([-(1+space),(1+1*space)])
%     ylim([-(1+1*space),(1+1*space)])
end


drawnow
pause(0.01)

if strcmp(app.plot_color,'default') && exist ('comm', 'var') && max(comm)>7
    warndlg(['The graph displayed has more communities ' ...
        'than colors available in the default color scheme. ' ...
        'Colors have been recycled. Choose another color' ...
        ' scheme to assign unique colors to each community.'])
end