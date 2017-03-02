function varargout = View_Edit_SPADE_tree_annotation(varargin)
% VIEW_EDIT_SPADE_TREE_ANNOTATION M-file for
% View_Edit_SPADE_tree_annotation.fig
%      VIEW_EDIT_SPADE_TREE_ANNOTATION, by itself, creates a new VIEW_EDIT_SPADE_TREE_ANNOTATION or raises the existing
%      singleton*.
%
%      H = VIEW_EDIT_SPADE_TREE_ANNOTATION returns the handle to a new VIEW_EDIT_SPADE_TREE_ANNOTATION or the handle to
%      the existing singleton*.
%
%      VIEW_EDIT_SPADE_TREE_ANNOTATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_EDIT_SPADE_TREE_ANNOTATION.M with the given input arguments.
%
%      VIEW_EDIT_SPADE_TREE_ANNOTATION('Property','Value',...) creates a new VIEW_EDIT_SPADE_TREE_ANNOTATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before View_Edit_SPADE_tree_annotation_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to View_Edit_SPADE_tree_annotation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help View_Edit_SPADE_tree_annotation

% Last Modified by GUIDE v2.5 19-Jun-2015 02:38:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @View_Edit_SPADE_tree_annotation_OpeningFcn, ...
                   'gui_OutputFcn',  @View_Edit_SPADE_tree_annotation_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before View_Edit_SPADE_tree_annotation is made visible.
function View_Edit_SPADE_tree_annotation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to View_Edit_SPADE_tree_annotation (see VARARGIN)

% Choose default command line output for View_Edit_SPADE_tree_annotation
handles.output = hObject;
if length(varargin)~=0
    handles.cluster_mst_upsample_filename =  varargin{1}{1};
else
    handles.cluster_mst_upsample_filename =  'SPADE_cluster_mst_upsample_result.mat';
end
load(handles.cluster_mst_upsample_filename);
% initialize the values available in this window
handles.all_fcs_filenames = all_fcs_filenames;
handles.all_assign = all_assign;
handles.file_annot = file_annot;
handles.clustering_data = data;
handles.clustering_idx  = idx;
handles.clustering_local_density = local_density;
handles.clustering_marker_names = marker_names;
handles.clustering_used_markers = used_markers;
handles.clustering_marker_node_average = marker_node_average;
handles.mst_tree = mst_tree;
handles.node_positions = node_positions;
handles.node_size = node_size;
handles.tree_annotations = tree_annotations;
handles.tree_bubble_contour = tree_bubble_contour;
if exist('user_defined_color_thresholds')
    handles.user_defined_color_thresholds = user_defined_color_thresholds;
else
    handles.user_defined_color_thresholds = zeros(0,4);
end

% get all file names from the marker_node_average cell array, preserve the order of the files
[C,I] = unique(marker_node_average(:,1)); handles.all_existing_files = marker_node_average(sort(I),1);
% get all marker names from the marker_node_average cell array, the code is different from above because I of unique give the last index not the first, so we need the following to preserve ordering of the markers 
handles.all_existing_markers = []; tmp = marker_node_average(:,2);
while length(tmp)~=0 
    handles.all_existing_markers = [handles.all_existing_markers; tmp(1)];
    tmp(ismember(tmp,tmp(1))==1)=[];
end

% internal control variables
handles.edge_handle = [];
handles.edge_begin_end = [];
handles.node_handle = [];
handles.GETRECT_H1 = [];
handles.GETRECT_H2 = [];
handles.mouse_selected_nodes = [];
handles.mouse_down_position=[];
handles.is_mouse_down = 0; % 0 means mouse not down, while 1 means mouse down
handles.is_mouse_down_and_moved=0;
handles.is_ctrl_down = 0;  
handles.is_shift_down = 0; 
handles.mouse_mode = []; % selection, move 
handles.show_annotation = 0;  % 0 = no show; 1 = show all; 2 = show selected
handles.color_definition = 0; % 0 = expr; 1 = ratio; 2 = cell freq
handles.color_scheme = 0;     % 0 = JET; 1 = half JET; 2 = gray scale
guidata(hObject, handles);
% % draw the tree
axes(handles.Axes_mst); handles = draw_SPADE_tree_when_open(handles);
set(handles.Axes_mst,'XTick',[],'YTick',[]);
guidata(hObject, handles);
% % update the parameters in this window
% % biaxial plots
set(handles.Popup_x,'string',handles.clustering_marker_names,'value',1);
set(handles.Popup_y,'string',handles.clustering_marker_names,'value',2);
update_baxial_plots(handles);
% % listboxes of marker, file, file, and ref
set(handles.listbox_select_a_marker,'string',handles.all_existing_markers,'value',1);
set(handles.listbox_select_a_file,'string',handles.all_existing_files,'value',1);
set(handles.listbox_select_ref_files,'string',handles.all_existing_files(1),'value',1);
% % tree annotations
if ~isempty(handles.tree_annotations)
    tmp = [];
    for i=1:length(handles.tree_annotations)
        tmp{i} = num2str(handles.tree_annotations{i});
    end
    set(handles.listbox_annotations,'string',tmp,'value',1);
end





% --- Outputs from this function are returned to the command line.
function varargout = View_Edit_SPADE_tree_annotation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% my own function to draw the two biaxial plots
function update_baxial_plots(handles)
tmp = get(handles.Popup_x,'string');
marker_x = tmp(get(handles.Popup_x,'value'));
tmp = get(handles.Popup_y,'string');
marker_y = tmp(get(handles.Popup_y,'value'));
x = handles.clustering_data(ismember(handles.clustering_marker_names,marker_x),:);
y = handles.clustering_data(ismember(handles.clustering_marker_names,marker_y),:);
ind = ((~isnan(x)) & (~isnan(y)));
x = x(ind);
y = y(ind);
idx = handles.clustering_idx(ind);
selected_cells = zeros(size(idx));
for i=1:length(handles.mouse_selected_nodes)
    selected_cells(idx==handles.mouse_selected_nodes(i))=1;
end
axes(handles.Axes_node_scatter);
plot(x,y,'g.',x(selected_cells==1),y(selected_cells==1),'b.');
axis([min(x)-(max(x)-min(x))*0.05,max(x)+(max(x)-min(x))*0.05,min(y)-(max(y)-min(y))*0.05,max(y)+(max(y)-min(y))*0.05])
axis_lim = axis;
axes(handles.Axes_node_contour);
hold off;
SPADE_contour2D(x,y);
hold on;
plot(x(selected_cells==1),y(selected_cells==1),'b.');
hold off;
axis(axis_lim);




% my own function to draw the SPADE tree according to the spec's in GUI
function handles = draw_SPADE_tree_when_open(handles)
% clear the plot
hold off; plot(0,0,'visible','off'); hold on;
% draw edges
adj = handles.mst_tree;
coeff = handles.node_positions;
pairs = SPADE_find_matrix_big_element(triu(adj,1),1);
for k=1:size(pairs,1), 
    handle_tmp = line(coeff(1,pairs(k,:)), coeff(2,pairs(k,:)),'color','g'); 
    handles.edge_handle = [handles.edge_handle; handle_tmp];
    handles.edge_begin_end = [handles.edge_begin_end; pairs(k,:)];
end

% % show annotation contour
%???NOTE show the bubbles 
switch handles.show_annotation % 0 = no show; 1 = show all; 2 = show selected
    case 0 
        'no show skip bubble';
    case 1 
        if length(handles.tree_annotations)~=0
            for k=1:length(handles.tree_bubble_contour)
                if ~isempty(handles.tree_bubble_contour{k})
                    outer_points = handles.tree_bubble_contour{k};
                    line([outer_points(:,1); outer_points(1,1)],[outer_points(:,2);outer_points(1,2)],'color',0.3+zeros(1,3));
                end
            end
        end
    case 2 
        if length(handles.tree_annotations)~=0
            k = get(handles.listbox_annotations,'value');
            if k~=0 && ~isempty(handles.tree_bubble_contour{k})
                outer_points = handles.tree_bubble_contour{k};
                line([outer_points(:,1); outer_points(1,1)],[outer_points(:,2);outer_points(1,2)],'color',0.3+zeros(1,3));
            end
        end
    otherwise
        'do nothing';
end

% % show node numbers
if get(handles.checkbox_show_node_index,'value')==1
    if isempty(handles.mouse_selected_nodes)
        for k=1:size(coeff,2), text(coeff(1,k)+2, coeff(2,k), num2str(k), 'FontSize', 7); end
    else
        for k=1:size(coeff,2)
            if ismember(k,handles.mouse_selected_nodes)==1
                text(coeff(1,k)+2, coeff(2,k), num2str(k), 'FontSize', 7);
            end
        end
    end
end


% % about color
node_size = handles.node_size;
node_color = zeros(length(node_size),3);
if get(handles.checkbox_is_use_color,'value')==0
    node_color = zeros(length(node_size),3) + 0.3;
else
    marker_file_selection_invalid = 0;
    % get color map
    switch handles.color_scheme     % 0 = JET; 1 = half JET; 2 = gray scale
        case 0
            cmap_tmp = get_JET_color_map;
        case 1
            cmap_tmp = get_half_JET_color_map;
        case 2
            cmap_tmp = get_gray_color_map;
    end

    switch handles.color_definition % 0 = expr; 1 = ratio; 2 = cell freq
        case 0
            tmp = get(handles.listbox_select_a_marker,'string');
            selected_marker = tmp(get(handles.listbox_select_a_marker,'value'));
            tmp = get(handles.listbox_select_a_file,'string');
            selected_file = tmp(get(handles.listbox_select_a_file,'value'));
            ind = find(ismember(handles.clustering_marker_node_average(:,1),selected_file)==1 & ismember(handles.clustering_marker_node_average(:,2),selected_marker)==1);
            if isempty(ind)
                marker_file_selection_invalid=1;
            else
                color_code_data = handles.clustering_marker_node_average{ind,3};
            end
            % normalize color code data
            if isequal(selected_marker, {'CellFreq'})
                color_code_data = color_code_data/sum(color_code_data)*100;
                prc95 = max(color_code_data);
                prc05 = min(color_code_data); 
            else
                if length(unique(color_code_data))<=3 || isequal(prctile(color_code_data,95), prctile(color_code_data,5))
                    prc95 = max(color_code_data);
                    prc05 = min(color_code_data); 
                else
                    prc95 = prctile(color_code_data,95);
                    prc05 = prctile(color_code_data,5); 
                end
            end
            color_code_data = (color_code_data - prc05)/(prc95-prc05);
            color_code_data(color_code_data<0.01)=0.01;
            color_code_data(color_code_data>0.99)=0.99;
            draw_colorbar([20,-55,20,5],cmap_tmp,{' ',' ',' '});
            text(38,-48,num2str(prc95,'%.2f')); text(18,-48,num2str(prc05,'%.2f'))
        case 1
            tmp = get(handles.listbox_select_a_marker,'string');
            selected_marker = tmp(get(handles.listbox_select_a_marker,'value'));
            tmp = get(handles.listbox_select_a_file,'string');
            selected_file = tmp(get(handles.listbox_select_a_file,'value'));
            ind = find(ismember(handles.clustering_marker_node_average(:,1),selected_file)==1 & ismember(handles.clustering_marker_node_average(:,2),selected_marker)==1);
            tmp = get(handles.listbox_select_ref_files,'string');
            ref_files = tmp(get(handles.listbox_select_ref_files,'value'));
            ind_ref = find(ismember(handles.clustering_marker_node_average(:,1),ref_files)==1 & ismember(handles.clustering_marker_node_average(:,2),selected_marker)==1);
            if isempty(ind) || isempty(ind_ref)
                marker_file_selection_invalid=1;
            else
                selected_file_data = handles.clustering_marker_node_average{ind,3};
                ref_file_data=[];
                for k=1:length(ind_ref)
                    ref_file_data = [ref_file_data;handles.clustering_marker_node_average{ind_ref(k),3}];
                end
                color_code_data = selected_file_data - nanmean(ref_file_data,1);
            end
            % normalize color code data
            if isequal(selected_marker, {'CellFreq'})
                color_code_data = (selected_file_data/sum(selected_file_data) - nanmean(ref_file_data./repmat(sum(ref_file_data,2),1,size(ref_file_data,2)),1))*100;
                prc95 = max(abs(color_code_data));
            else
                prc95 = prctile(abs(color_code_data),95);
            end
            color_code_data = (color_code_data/prc95 + 1)/2;
            color_code_data(color_code_data<0.01)=0.01;
            color_code_data(color_code_data>0.99)=0.99;
            draw_colorbar([20,-55,20,5],cmap_tmp,{' ',' ',' '});
            text(38,-48,num2str(prc95,'%.2f')); text(18,-48,num2str(-prc95,'%.2f'))
        case 2
            tmp = get(handles.listbox_select_a_file,'string');
            selected_file = tmp(get(handles.listbox_select_a_file,'value'));
            if get(handles.listbox_select_a_file,'value')==1 % selected the pooled file
                [dummy, color_code_data] = SPADE_compute_one_marker_group_mean(ones(1,length(handles.clustering_idx)),handles.clustering_idx);    
            else
                ind = find(ismember(handles.file_annot,selected_file)==1);
                % [dummy, color_code_data] = SPADE_compute_one_marker_group_mean(ones(1,length(handles.all_assign{ind})),handles.all_assign{ind});    
                [group_avg, counts, group_idx_values] = SPADE_compute_one_marker_group_mean(ones(1,length(handles.all_assign{ind})),handles.all_assign{ind}); % the following few lines are for the purpose that: some file may not have any cell belong to one particular node, and therefore, the "group_avg" does not have information for every node
                color_code_data = zeros(1,max(handles.clustering_idx));
                color_code_data(group_idx_values)=counts;
            end
            color_code_data = color_code_data/sum(color_code_data)*100; % normalize from counts to percentage
            % normalize color code data
            prc95 = max(color_code_data);
            prc05 = min(color_code_data); 
            color_code_data = (color_code_data - prc05)/(prc95-prc05);
            color_code_data(color_code_data<0.01)=0.01;
            color_code_data(color_code_data>0.99)=0.99;
            draw_colorbar([20,-55,20,5],cmap_tmp,{' ',' ',' '});
            text(38,-48,num2str(prc95,'%.2f')); text(18,-48,num2str(prc05,'%.2f'))
    end

    if marker_file_selection_invalid==1 % this selection is invalid, the data that color codes according to the selection does not exist
        node_color = zeros(length(node_size),3); 
        title('invalid marker file selection');
    elseif length(unique(color_code_data))==1
        node_color = zeros(length(node_size),3); 
        title('no variation');
    else
        % compute the color according to color_code_data, using the settings of color scheme
        % get color map
        switch handles.color_scheme     % 0 = JET; 1 = half JET; 2 = gray scale
            case 0
                cmap_tmp = get_JET_color_map;
            case 1
                cmap_tmp = get_half_JET_color_map;
            case 2
                cmap_tmp = get_gray_color_map;
        end
        for k=1:size(coeff,2), 
            node_color(k,:) = interp1(((1:size(cmap_tmp,1))'-1)/(size(cmap_tmp,1)-1),cmap_tmp,color_code_data(k));  
            if sum(isnan(node_color(k,:)))~=0
                node_color(k,:) = [1,1,1];
            end
        end
    end
end



% % draw nodes
if isempty(handles.mouse_selected_nodes)
    for k=1:length(node_size)
        handle_tmp = plot(coeff(1,k),coeff(2,k),'o','markersize',node_size(k), 'color',node_color(k,:), 'markerfacecolor',node_color(k,:),'markeredgecolor',node_color(k,:));
        set(handle_tmp,'ButtonDownFcn','View_Edit_SPADE_tree_annotation(''edit_tree_NodeButtonDownFcn'',gcbo,[],guidata(gcbo))');
        handles.node_handle = [handles.node_handle; handle_tmp];
    end
else
    for k=1:length(node_size)
        if ismember(k,handles.mouse_selected_nodes)==0
            handle_tmp = plot(coeff(1,k),coeff(2,k),'o','markersize',node_size(k), 'markerfacecolor',node_color(k,:), 'markerfacecolor',node_color(k,:),'markeredgecolor',node_color(k,:));
            if ismember(handle_tmp,get(handles.Axes_mst,'child'))
                set(handle_tmp,'ButtonDownFcn','View_Edit_SPADE_tree_annotation(''edit_tree_NodeButtonDownFcn'',gcbo,[],guidata(gcbo))');
            end
            handles.node_handle = [handles.node_handle; handle_tmp];
        else
            % this node is selected, and needs to be processed differently
            handle_tmp = plot(coeff(1,k),coeff(2,k),'s','markersize',node_size(k), 'markerfacecolor',node_color(k,:), 'markerfacecolor',node_color(k,:),'markeredgecolor',[0,0,0]);
            if ismember(handle_tmp,get(handles.Axes_mst,'child'))
                set(handle_tmp,'ButtonDownFcn','View_Edit_SPADE_tree_annotation(''edit_tree_NodeButtonDownFcn'',gcbo,[],guidata(gcbo))');
            end
            handles.node_handle = [handles.node_handle; handle_tmp];
        end
    end
end

% % control axis limits
axis_lims = reshape([-max(abs(coeff)'); +max(abs(coeff)')], 1, 4)*1.1;
for i=1:4
    if abs(axis_lims(i))<55
        if axis_lims(i)>=0
            axis_lims(i)=55;
        else
            axis_lims(i)=-55;
        end
    end
end
axis(axis_lims);
set(handles.Axes_mst,'ButtonDownFcn','View_Edit_SPADE_tree_annotation(''Axes_mst_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
% guidata(handles.button_show_tree_new_window,handles);




% my own function to draw the SPADE tree according to the spec's in GUI
function draw_SPADE_tree_when_update(handles)
% components to be cleared: anything other than the edges, nodes, and two boxes for mouse selection 
components_to_clear = setdiff(get(handles.Axes_mst,'children'),[handles.node_handle;handles.edge_handle;handles.GETRECT_H1;handles.GETRECT_H2]);
for i=1:length(components_to_clear)
    delete(components_to_clear(i));
end
% update edges
adj = handles.mst_tree;
coeff = handles.node_positions;
pairs = SPADE_find_matrix_big_element(triu(adj,1),1);
for k=1:size(pairs,1), 
    handle_tmp = handles.edge_handle(k);
    if ~isequal(get(handle_tmp,'XData'),coeff(1,pairs(k,:))) ||  ~isequal(get(handle_tmp,'YData'),coeff(2,pairs(k,:)))
        set(handle_tmp,'XData',coeff(1,pairs(k,:)),'YData',coeff(2,pairs(k,:)));
    end
end

% % show annotation contour
%???NOTE show the bubbles 
switch handles.show_annotation % 0 = no show; 1 = show all; 2 = show selected
    case 0 
        'no show skip bubble';
    case 1 
        if length(handles.tree_annotations)~=0
            for k=1:length(handles.tree_bubble_contour)
                if ~isempty(handles.tree_bubble_contour{k})
                    outer_points = handles.tree_bubble_contour{k};
                    line([outer_points(:,1); outer_points(1,1)],[outer_points(:,2);outer_points(1,2)],'color',0.3+zeros(1,3));
                end
            end
        end
    case 2 
        if length(handles.tree_annotations)~=0
            k = get(handles.listbox_annotations,'value');
            if k~=0 && ~isempty(handles.tree_bubble_contour{k})
                outer_points = handles.tree_bubble_contour{k};
                line([outer_points(:,1); outer_points(1,1)],[outer_points(:,2);outer_points(1,2)],'color',0.3+zeros(1,3));
            end
        end
    otherwise
        'do nothing';
end

% % show node numbers
if get(handles.checkbox_show_node_index,'value')==1
    if isempty(handles.mouse_selected_nodes)
        for k=1:size(coeff,2), text(coeff(1,k)+2, coeff(2,k), num2str(k), 'FontSize', 7); end
    else
        for k=1:size(coeff,2)
            if ismember(k,handles.mouse_selected_nodes)==1
                text(coeff(1,k)+2, coeff(2,k), num2str(k), 'FontSize', 7);
            end
        end
    end
end


% % about color
node_size = handles.node_size;
node_color = zeros(length(node_size),3);
if get(handles.checkbox_is_use_color,'value')==0
    node_color = zeros(length(node_size),3) + 0.3;
    set(handles.edit_lower_threshold,'enable','off','string','');
    set(handles.edit_upper_threshold,'enable','off','string','');
    set(handles.button_default_colorbar_threshold,'enable','off');
else
    marker_file_selection_invalid = 0;
    % get color map
    switch handles.color_scheme     % 0 = JET; 1 = half JET; 2 = gray scale
        case 0
            cmap_tmp = get_JET_color_map;
        case 1
            cmap_tmp = get_half_JET_color_map;
        case 2
            cmap_tmp = get_gray_color_map;
    end

    switch handles.color_definition % 0 = expr; 1 = ratio; 2 = cell freq
        case 0
            tmp = get(handles.listbox_select_a_marker,'string');
            selected_marker = tmp(get(handles.listbox_select_a_marker,'value'));
            tmp = get(handles.listbox_select_a_file,'string');
            selected_file = tmp(get(handles.listbox_select_a_file,'value'));
            ind = find(ismember(handles.clustering_marker_node_average(:,1),selected_file)==1 & ismember(handles.clustering_marker_node_average(:,2),selected_marker)==1);
            if isempty(ind)
                marker_file_selection_invalid=1;
            else
                color_code_data = handles.clustering_marker_node_average{ind,3};
            end
            % normalize color code data
            if isequal(selected_marker, {'CellFreq'})
                color_code_data = color_code_data/sum(color_code_data)*100;
                prc95 = max(color_code_data);
                prc05 = min(color_code_data); 
            else
                if length(unique(color_code_data))<=3 || isequal(prctile(color_code_data,95), prctile(color_code_data,5))
                    prc95 = max(color_code_data);
                    prc05 = min(color_code_data); 
                else
                    prc95 = prctile(color_code_data,95);
                    prc05 = prctile(color_code_data,5); 
                end
            end
            
            ind_user_defined_thresholds = find(ismember(handles.user_defined_color_thresholds(:,1),handles.color_definition) & ismember(handles.user_defined_color_thresholds(:,2),get(handles.listbox_select_a_marker,'value')));
            if ~isempty(ind_user_defined_thresholds)
                prc05 = handles.user_defined_color_thresholds(ind_user_defined_thresholds,3);
                prc95 = handles.user_defined_color_thresholds(ind_user_defined_thresholds,4);
                set(handles.button_default_colorbar_threshold,'enable','on');
            else
                set(handles.button_default_colorbar_threshold,'enable','off');
            end
            
            color_code_data = (color_code_data - prc05)/(prc95-prc05);
            color_code_data(color_code_data<0.01)=0.01;
            color_code_data(color_code_data>0.99)=0.99;
            draw_colorbar([20,-55,20,5],cmap_tmp,{' ',' ',' '});
            text(38,-48,num2str(prc95,'%.2f')); text(18,-48,num2str(prc05,'%.2f'))
            set(handles.edit_lower_threshold,'enable','on','string',num2str(prc05,'%.2f'),'value',prc05);
            set(handles.edit_upper_threshold,'enable','on','string',num2str(prc95,'%.2f'),'value',prc95);
        case 1
            tmp = get(handles.listbox_select_a_marker,'string');
            selected_marker = tmp(get(handles.listbox_select_a_marker,'value'));
            tmp = get(handles.listbox_select_a_file,'string');
            selected_file = tmp(get(handles.listbox_select_a_file,'value'));
            ind = find(ismember(handles.clustering_marker_node_average(:,1),selected_file)==1 & ismember(handles.clustering_marker_node_average(:,2),selected_marker)==1);
            tmp = get(handles.listbox_select_ref_files,'string');
            ref_files = tmp(get(handles.listbox_select_ref_files,'value'));
            ind_ref = find(ismember(handles.clustering_marker_node_average(:,1),ref_files)==1 & ismember(handles.clustering_marker_node_average(:,2),selected_marker)==1);
            if isempty(ind) || isempty(ind_ref)
                marker_file_selection_invalid=1;
            else
                selected_file_data = handles.clustering_marker_node_average{ind,3};
                ref_file_data=[];
                for k=1:length(ind_ref)
                    ref_file_data = [ref_file_data;handles.clustering_marker_node_average{ind_ref(k),3}];
                end
                color_code_data = selected_file_data - nanmean(ref_file_data,1);
            end
            % normalize color code data
            if isequal(selected_marker, {'CellFreq'})
                color_code_data = (selected_file_data/sum(selected_file_data) - nanmean(ref_file_data./repmat(sum(ref_file_data,2),1,size(ref_file_data,2)),1))*100;
                prc95 = max(abs(color_code_data));
            else
                prc95 = prctile(abs(color_code_data),95);
            end
            
            ind_user_defined_thresholds = find(ismember(handles.user_defined_color_thresholds(:,1),handles.color_definition) & ismember(handles.user_defined_color_thresholds(:,2),get(handles.listbox_select_a_marker,'value')));
            if ~isempty(ind_user_defined_thresholds) && prc95~=0
                prc95 = handles.user_defined_color_thresholds(ind_user_defined_thresholds,4);
                set(handles.button_default_colorbar_threshold,'enable','on');
            else
                set(handles.button_default_colorbar_threshold,'enable','off');
            end
            
            color_code_data = (color_code_data/prc95 + 1)/2;
            color_code_data(color_code_data<0.01)=0.01;
            color_code_data(color_code_data>0.99)=0.99;
            draw_colorbar([20,-55,20,5],cmap_tmp,{' ',' ',' '});
            text(38,-48,num2str(prc95,'%.2f')); text(18,-48,num2str(-prc95,'%.2f'))
            set(handles.edit_lower_threshold,'enable','on','string',num2str(-prc95,'%.2f'),'value',-prc95);
            set(handles.edit_upper_threshold,'enable','on','string',num2str(prc95,'%.2f'),'value',prc95);
        case 2
            tmp = get(handles.listbox_select_a_file,'string');
            selected_file = tmp(get(handles.listbox_select_a_file,'value'));
            if get(handles.listbox_select_a_file,'value')==1 % selected the pooled file
                [dummy, color_code_data] = SPADE_compute_one_marker_group_mean(ones(1,length(handles.clustering_idx)),handles.clustering_idx);    
            else
                ind = find(ismember(handles.file_annot,selected_file)==1);
                % [dummy, color_code_data] = SPADE_compute_one_marker_group_mean(ones(1,length(handles.all_assign{ind})),handles.all_assign{ind});    
                [group_avg, counts, group_idx_values] = SPADE_compute_one_marker_group_mean(ones(1,length(handles.all_assign{ind})),handles.all_assign{ind}); % the following few lines are for the purpose that: some file may not have any cell belong to one particular node, and therefore, the "group_avg" does not have information for every node
                color_code_data = zeros(1,max(handles.clustering_idx));
                color_code_data(group_idx_values)=counts;
            end
            color_code_data = color_code_data/sum(color_code_data)*100; % normalize from counts to percentage
            % normalize color code data
            prc95 = max(color_code_data);
            prc05 = min(color_code_data); 

            ind_user_defined_thresholds = find(ismember(handles.user_defined_color_thresholds(:,1),handles.color_definition) & ismember(handles.user_defined_color_thresholds(:,2),get(handles.listbox_select_a_marker,'value')));
            if ~isempty(ind_user_defined_thresholds)
                prc05 = handles.user_defined_color_thresholds(ind_user_defined_thresholds,3);
                prc95 = handles.user_defined_color_thresholds(ind_user_defined_thresholds,4);
                set(handles.button_default_colorbar_threshold,'enable','on');
            else
                set(handles.button_default_colorbar_threshold,'enable','off');
            end

            
            color_code_data = (color_code_data - prc05)/(prc95-prc05);
            color_code_data(color_code_data<0.01)=0.01;
            color_code_data(color_code_data>0.99)=0.99;
            draw_colorbar([20,-55,20,5],cmap_tmp,{' ',' ',' '});
            text(38,-48,num2str(prc95,'%.2f')); text(18,-48,num2str(prc05,'%.2f'))
            set(handles.edit_lower_threshold,'enable','on','string',num2str(prc05,'%.2f'),'value',prc05);
            set(handles.edit_upper_threshold,'enable','on','string',num2str(prc95,'%.2f'),'value',prc95);
    end

    if marker_file_selection_invalid==1 % this selection is invalid, the data that color codes according to the selection does not exist
        node_color = zeros(length(node_size),3); 
        title('invalid marker file selection');
    elseif length(unique(color_code_data))==1
        node_color = zeros(length(node_size),3); 
        title('no variation');
    else
        % compute the color according to color_code_data, using the settings of color scheme
        % get color map
        switch handles.color_scheme     % 0 = JET; 1 = half JET; 2 = gray scale
            case 0
                cmap_tmp = get_JET_color_map;
            case 1
                cmap_tmp = get_half_JET_color_map;
            case 2
                cmap_tmp = get_gray_color_map;
        end
        for k=1:size(coeff,2), 
            node_color(k,:) = interp1(((1:size(cmap_tmp,1))'-1)/(size(cmap_tmp,1)-1),cmap_tmp,color_code_data(k));  
            if sum(isnan(node_color(k,:)))~=0
                node_color(k,:) = [1,1,1];
            end
        end
    end
end



% % draw nodes
if isempty(handles.mouse_selected_nodes)
    for k=1:length(node_size)
        handle_tmp = handles.node_handle(k);
        if ~isequal(get(handle_tmp,'XData'),coeff(1,k)) || ~isequal(get(handle_tmp,'YData'),coeff(2,k)) || ~isequal(get(handle_tmp,'Marker'),'o') || ~isequal(get(handle_tmp,'MarkerSize'),node_size(k)) || ~isequal(get(handle_tmp,'color'),node_color(k,:)) || ~isequal(get(handle_tmp,'markerfaceColor'),node_color(k,:)) || ~isequal(get(handle_tmp,'markeredgecolor'),node_color(k,:)) 
            set(handle_tmp,'XData',coeff(1,k),'YData',coeff(2,k),'Marker','o','markersize',node_size(k), 'color', node_color(k,:), 'markerfacecolor',node_color(k,:),'markeredgecolor',node_color(k,:));
        end
    end
else
    for k=1:length(node_size)
        if ismember(k,handles.mouse_selected_nodes)==0
            handle_tmp = handles.node_handle(k);
            if ~isequal(get(handle_tmp,'XData'),coeff(1,k)) || ~isequal(get(handle_tmp,'YData'),coeff(2,k)) || ~isequal(get(handle_tmp,'Marker'),'o') || ~isequal(get(handle_tmp,'MarkerSize'),node_size(k)) || ~isequal(get(handle_tmp,'color'),node_color(k,:)) || ~isequal(get(handle_tmp,'markerfaceColor'),node_color(k,:)) || ~isequal(get(handle_tmp,'markeredgecolor'),node_color(k,:)) 
                set(handle_tmp,'XData',coeff(1,k),'YData',coeff(2,k),'Marker','o','markersize',node_size(k), 'color', node_color(k,:), 'markerfacecolor',node_color(k,:),'markeredgecolor',node_color(k,:));
            end
        else
            % this node is selected, and needs to be processed differently
            handle_tmp = handles.node_handle(k);
            if ~isequal(get(handle_tmp,'XData'),coeff(1,k)) || ~isequal(get(handle_tmp,'YData'),coeff(2,k)) || ~isequal(get(handle_tmp,'Marker'),'s') || ~isequal(get(handle_tmp,'MarkerSize'),node_size(k)) || ~isequal(get(handle_tmp,'color'),node_color(k,:)) || ~isequal(get(handle_tmp,'markerfaceColor'),node_color(k,:)) || ~isequal(get(handle_tmp,'markeredgecolor'),[0 0 0]) 
                set(handle_tmp,'XData',coeff(1,k),'YData',coeff(2,k),'Marker','s','markersize',node_size(k), 'color', node_color(k,:), 'markerfacecolor',node_color(k,:),'markeredgecolor',[0 0 0]);
            end
        end
    end
end

% % control axis limits
axis_lims = reshape([-max(abs(coeff)'); +max(abs(coeff)')], 1, 4)*1.1;
for i=1:4
    if abs(axis_lims(i))<55
        if axis_lims(i)>=0
            axis_lims(i)=55;
        else
            axis_lims(i)=-55;
        end
    end
end
axis(axis_lims);
set(handles.Axes_mst,'ButtonDownFcn','View_Edit_SPADE_tree_annotation(''Axes_mst_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
guidata(handles.button_show_tree_new_window,handles);





function SPADE_cmap = get_JET_color_map

SPADE_cmap = [   0         0    0.5625
                 0         0    0.6250
                 0         0    0.6875
                 0         0    0.7500
                 0         0    0.8125
                 0         0    0.8750
                 0         0    0.9375
                 0         0    1.0000
                 0    0.0625    1.0000
                 0    0.1250    1.0000
                 0    0.1875    1.0000
                 0    0.2500    1.0000
                 0    0.3125    1.0000
                 0    0.3750    1.0000
                 0    0.4375    1.0000
                 0    0.5000    1.0000
                 0    0.5625    1.0000
                 0    0.6250    1.0000
                 0    0.6875    1.0000
                 0    0.7500    1.0000
                 0    0.8125    1.0000
                 0    0.8750    1.0000
                 0    0.9375    1.0000
                 0    1.0000    1.0000
            0.0625    1.0000    0.9375
            0.1250    1.0000    0.8750
            0.1875    1.0000    0.8125
            0.2500    1.0000    0.7500
            0.3125    1.0000    0.6875
            0.3750    1.0000    0.6250
            0.4375    1.0000    0.5625
            0.5000    1.0000    0.5000
            0.5625    1.0000    0.4375
            0.6250    1.0000    0.3750
            0.6875    1.0000    0.3125
            0.7500    1.0000    0.2500
            0.8125    1.0000    0.1875
            0.8750    1.0000    0.1250
            0.9375    1.0000    0.0625
            1.0000    1.0000         0
            1.0000    0.9375         0
            1.0000    0.8750         0
            1.0000    0.8125         0
            1.0000    0.7500         0
            1.0000    0.6875         0
            1.0000    0.6250         0
            1.0000    0.5625         0
            1.0000    0.5000         0
            1.0000    0.4375         0
            1.0000    0.3750         0
            1.0000    0.3125         0
            1.0000    0.2500         0
            1.0000    0.1875         0
            1.0000    0.1250         0
            1.0000    0.0625         0
            1.0000         0         0
            0.9375         0         0
            0.8750         0         0
            0.8125         0         0
            0.7500         0         0
            0.6875         0         0
            0.6250         0         0
            0.5625         0         0
            0.5000         0         0];
        
    
function SPADE_cmap = get_half_JET_color_map

SPADE_cmap = [  0.5625    1.0000    0.4375
                0.6250    1.0000    0.3750
                0.6875    1.0000    0.3125
                0.7500    1.0000    0.2500
                0.8125    1.0000    0.1875
                0.8750    1.0000    0.1250
                0.9375    1.0000    0.0625
                1.0000    1.0000         0
                1.0000    0.9375         0
                1.0000    0.8750         0
                1.0000    0.8125         0
                1.0000    0.7500         0
                1.0000    0.6875         0
                1.0000    0.6250         0
                1.0000    0.5625         0
                1.0000    0.5000         0
                1.0000    0.4375         0
                1.0000    0.3750         0
                1.0000    0.3125         0
                1.0000    0.2500         0
                1.0000    0.1875         0
                1.0000    0.1250         0
                1.0000    0.0625         0
                1.0000         0         0
                0.9375         0         0
                0.8750         0         0
                0.8125         0         0
                0.7500         0         0
                0.6875         0         0
                0.6250         0         0
                0.5625         0         0
                0.5000         0         0];
        
    
function SPADE_cmap = get_gray_color_map
baseMap = [220   220    220;
           [165   165    165]*1.4;
           [110   110    110]*1.4;
           [55    55     55]*1.4;
           1     1      1]/255*(10/11);
idx1 = linspace(0,1,size(baseMap,1));
idx2 = linspace(0,1,100);
SPADE_cmap = interp1(idx1,baseMap,idx2);
   



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% my own function to draw the SPADE tree according to the spec's in GUI, but with no node mouse click events 
function draw_SPADE_tree_no_mouse_click_events(handles)
% clear the plot
hold off; plot(0,0,'visible','off'); hold on;
% draw edges
adj = handles.mst_tree;
coeff = handles.node_positions;
pairs = SPADE_find_matrix_big_element(triu(adj,1),1);
for k=1:size(pairs,1), 
    handle_tmp = line(coeff(1,pairs(k,:)), coeff(2,pairs(k,:)),'color','g'); 
end

% % show annotation contour
%???NOTE show the bubbles 
switch handles.show_annotation % 0 = no show; 1 = show all; 2 = show selected
    case 0 
        'no show skip bubble';
    case 1 
        if length(handles.tree_annotations)~=0
            for k=1:length(handles.tree_bubble_contour)
                if ~isempty(handles.tree_bubble_contour{k})
                    outer_points = handles.tree_bubble_contour{k};
                    line([outer_points(:,1); outer_points(1,1)],[outer_points(:,2);outer_points(1,2)],'color',0.3+zeros(1,3));
                end
            end
        end
    case 2 
        if length(handles.tree_annotations)~=0
            k = get(handles.listbox_annotations,'value');
            if k~=0 && ~isempty(handles.tree_bubble_contour{k})
                outer_points = handles.tree_bubble_contour{k};
                line([outer_points(:,1); outer_points(1,1)],[outer_points(:,2);outer_points(1,2)],'color',0.3+zeros(1,3));
            end
        end
    otherwise
        'do nothing';
end

% % show node numbers
if get(handles.checkbox_show_node_index,'value')==1
    if isempty(handles.mouse_selected_nodes)
        for k=1:size(coeff,2), text(coeff(1,k)+2, coeff(2,k), num2str(k), 'FontSize', 7); end
    else
        for k=1:size(coeff,2)
            if ismember(k,handles.mouse_selected_nodes)==1
                text(coeff(1,k)+2, coeff(2,k), num2str(k), 'FontSize', 7);
            end
        end
    end
end


% % about color
node_size = handles.node_size;
node_color = zeros(length(node_size),3);
if get(handles.checkbox_is_use_color,'value')==0
    node_color = zeros(length(node_size),3) + 0.3;
else
    marker_file_selection_invalid = 0;
    % get color map
    switch handles.color_scheme     % 0 = JET; 1 = half JET; 2 = gray scale
        case 0
            cmap_tmp = get_JET_color_map;
        case 1
            cmap_tmp = get_half_JET_color_map;
        case 2
            cmap_tmp = get_gray_color_map;
    end

    switch handles.color_definition % 0 = expr; 1 = ratio; 2 = cell freq
        case 0
            tmp = get(handles.listbox_select_a_marker,'string');
            selected_marker = tmp(get(handles.listbox_select_a_marker,'value'));
            tmp = get(handles.listbox_select_a_file,'string');
            selected_file = tmp(get(handles.listbox_select_a_file,'value'));
            ind = find(ismember(handles.clustering_marker_node_average(:,1),selected_file)==1 & ismember(handles.clustering_marker_node_average(:,2),selected_marker)==1);
            if isempty(ind)
                marker_file_selection_invalid=1;
            else
                color_code_data = handles.clustering_marker_node_average{ind,3};
            end
            % normalize color code data
            if isequal(selected_marker, {'CellFreq'})
                color_code_data = color_code_data/sum(color_code_data)*100;
                prc95 = max(color_code_data);
                prc05 = min(color_code_data); 
            else
                if length(unique(color_code_data))<=3
                    prc95 = max(color_code_data);
                    prc05 = min(color_code_data); 
                else
                    prc95 = prctile(color_code_data,95);
                    prc05 = prctile(color_code_data,5); 
                end
            end
            
            ind_user_defined_thresholds = find(ismember(handles.user_defined_color_thresholds(:,1),handles.color_definition) & ismember(handles.user_defined_color_thresholds(:,2),get(handles.listbox_select_a_marker,'value')));
            if ~isempty(ind_user_defined_thresholds)
                prc05 = handles.user_defined_color_thresholds(ind_user_defined_thresholds,3);
                prc95 = handles.user_defined_color_thresholds(ind_user_defined_thresholds,4);
            end
            
            
            color_code_data = (color_code_data - prc05)/(prc95-prc05);
            color_code_data(color_code_data<0.01)=0.01;
            color_code_data(color_code_data>0.99)=0.99;
            draw_colorbar([20,-55,20,5],cmap_tmp,{' ',' ',' '});
            text(38,-48,num2str(prc95,'%.2f')); text(18,-48,num2str(prc05,'%.2f'))
        case 1
            tmp = get(handles.listbox_select_a_marker,'string');
            selected_marker = tmp(get(handles.listbox_select_a_marker,'value'));
            tmp = get(handles.listbox_select_a_file,'string');
            selected_file = tmp(get(handles.listbox_select_a_file,'value'));
            ind = find(ismember(handles.clustering_marker_node_average(:,1),selected_file)==1 & ismember(handles.clustering_marker_node_average(:,2),selected_marker)==1);
            tmp = get(handles.listbox_select_ref_files,'string');
            ref_files = tmp(get(handles.listbox_select_ref_files,'value'));
            ind_ref = find(ismember(handles.clustering_marker_node_average(:,1),ref_files)==1 & ismember(handles.clustering_marker_node_average(:,2),selected_marker)==1);
            if isempty(ind) || isempty(ind_ref)
                marker_file_selection_invalid=1;
            else
                selected_file_data = handles.clustering_marker_node_average{ind,3};
                ref_file_data=[];
                for k=1:length(ind_ref)
                    ref_file_data = [ref_file_data;handles.clustering_marker_node_average{ind_ref(k),3}];
                end
                color_code_data = selected_file_data - nanmean(ref_file_data,1);
            end
            % normalize color code data
            if isequal(selected_marker, {'CellFreq'})
                color_code_data = (selected_file_data/sum(selected_file_data) - nanmean(ref_file_data./repmat(sum(ref_file_data,2),1,size(ref_file_data,2)),1))*100;
                prc95 = max(abs(color_code_data));
            else
                prc95 = prctile(abs(color_code_data),95);
            end
            
            ind_user_defined_thresholds = find(ismember(handles.user_defined_color_thresholds(:,1),handles.color_definition) & ismember(handles.user_defined_color_thresholds(:,2),get(handles.listbox_select_a_marker,'value')));
            if ~isempty(ind_user_defined_thresholds) && prc95~=0
                prc95 = handles.user_defined_color_thresholds(ind_user_defined_thresholds,4);
            end
            
            color_code_data = (color_code_data/prc95 + 1)/2;
            color_code_data(color_code_data<0.01)=0.01;
            color_code_data(color_code_data>0.99)=0.99;
            draw_colorbar([20,-55,20,5],cmap_tmp,{' ',' ',' '});
            text(38,-48,num2str(prc95,'%.2f')); text(18,-48,num2str(-prc95,'%.2f'))
        case 2
            tmp = get(handles.listbox_select_a_file,'string');
            selected_file = tmp(get(handles.listbox_select_a_file,'value'));
            if get(handles.listbox_select_a_file,'value')==1 % selected the pooled file
                [dummy, color_code_data] = SPADE_compute_one_marker_group_mean(ones(1,length(handles.clustering_idx)),handles.clustering_idx);    
            else
                ind = find(ismember(handles.file_annot,selected_file)==1);
                % [dummy, color_code_data] = SPADE_compute_one_marker_group_mean(ones(1,length(handles.all_assign{ind})),handles.all_assign{ind});    
                [group_avg, counts, group_idx_values] = SPADE_compute_one_marker_group_mean(ones(1,length(handles.all_assign{ind})),handles.all_assign{ind}); % the following few lines are for the purpose that: some file may not have any cell belong to one particular node, and therefore, the "group_avg" does not have information for every node
                color_code_data = zeros(1,max(handles.clustering_idx));
                color_code_data(group_idx_values)=counts;
            end
            color_code_data = color_code_data/sum(color_code_data)*100; % normalize from counts to percentage
            % normalize color code data
            prc95 = max(color_code_data);
            prc05 = min(color_code_data); 

            ind_user_defined_thresholds = find(ismember(handles.user_defined_color_thresholds(:,1),handles.color_definition) & ismember(handles.user_defined_color_thresholds(:,2),get(handles.listbox_select_a_marker,'value')));
            if ~isempty(ind_user_defined_thresholds)
                prc05 = handles.user_defined_color_thresholds(ind_user_defined_thresholds,3);
                prc95 = handles.user_defined_color_thresholds(ind_user_defined_thresholds,4);
            end
            
            
            color_code_data = (color_code_data - prc05)/(prc95-prc05);
            color_code_data(color_code_data<0.01)=0.01;
            color_code_data(color_code_data>0.99)=0.99;
            draw_colorbar([20,-55,20,5],cmap_tmp,{' ',' ',' '});
            text(38,-48,num2str(prc95,'%.2f')); text(18,-48,num2str(prc05,'%.2f'))
    end

    if marker_file_selection_invalid==1 % this selection is invalid, the data that color codes according to the selection does not exist
        node_color = zeros(length(node_size),3); 
        title('invalid marker file selection');
    elseif length(unique(color_code_data))==1
        node_color = zeros(length(node_size),3); 
        title('no variation');
    else
        % compute the color according to color_code_data, using the settings of color scheme
        % get color map
        switch handles.color_scheme     % 0 = JET; 1 = half JET; 2 = gray scale
            case 0
                cmap_tmp = get_JET_color_map;
            case 1
                cmap_tmp = get_half_JET_color_map;
            case 2
                cmap_tmp = get_gray_color_map;
        end
        for k=1:size(coeff,2), 
            node_color(k,:) = interp1(((1:size(cmap_tmp,1))'-1)/(size(cmap_tmp,1)-1),cmap_tmp,color_code_data(k));  
            if sum(isnan(node_color(k,:)))~=0
                node_color(k,:) = [1,1,1];
            end
        end
    end
end



% % draw nodes
if isempty(handles.mouse_selected_nodes)
    for k=1:length(node_size)
        handle_tmp = plot(coeff(1,k),coeff(2,k),'o','markersize',node_size(k), 'color',node_color(k,:), 'markerfacecolor',node_color(k,:),'markeredgecolor',node_color(k,:));
    end
else
    for k=1:length(node_size)
        if ismember(k,handles.mouse_selected_nodes)==0
            handle_tmp = plot(coeff(1,k),coeff(2,k),'o','markersize',node_size(k), 'markerfacecolor',node_color(k,:), 'markerfacecolor',node_color(k,:),'markeredgecolor',node_color(k,:));
        else
            % this node is selected, and needs to be processed differently
            handle_tmp = plot(coeff(1,k),coeff(2,k),'s','markersize',node_size(k), 'markerfacecolor',node_color(k,:), 'markerfacecolor',node_color(k,:),'markeredgecolor',[0,0,0]);
        end
    end
end

% % control axis limits
axis_lims = reshape([-max(abs(coeff)'); +max(abs(coeff)')], 1, 4)*1.1;
for i=1:4
    if abs(axis_lims(i))<55
        if axis_lims(i)>=0
            axis_lims(i)=55;
        else
            axis_lims(i)=-55;
        end
    end
end
axis(axis_lims);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% --- Executes on selection change in Popup_y.
function Popup_y_Callback(hObject, eventdata, handles)
update_baxial_plots(handles);


% --- Executes during object creation, after setting all properties.
function Popup_y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Popup_x.
function Popup_x_Callback(hObject, eventdata, handles)
update_baxial_plots(handles);

% --- Executes during object creation, after setting all properties.
function Popup_x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_update_plots.
function button_update_plots_Callback(hObject, eventdata, handles)
update_baxial_plots(handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on button press in checkbox_show_node_index.
function checkbox_show_node_index_Callback(hObject, eventdata, handles)
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);


% --- Executes on button press in button_show_tree_new_window.
function button_show_tree_new_window_Callback(hObject, eventdata, handles)
figure;draw_SPADE_tree_no_mouse_click_events(handles);
axis off;
s = get(handles.listbox_select_a_marker,'string');title(s{get(handles.listbox_select_a_marker,'value')},'interpret','none') ;


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in listbox_annotations.
function listbox_annotations_Callback(hObject, eventdata, handles)
tmp = get(handles.listbox_annotations,'string');
if length(tmp)==0
    return
end
handles.mouse_selected_nodes = str2num(tmp{get(handles.listbox_annotations,'value')});
set(handles.edit_mouse_selected_nodes,'string',num2str(handles.mouse_selected_nodes));
guidata(hObject, handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);




% --- Executes during object creation, after setting all properties.
function listbox_annotations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_annotations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

 
% % --- Executes on button press in radiobutton_annotation_no_show.
% function radiobutton_annotation_no_show_Callback(hObject, eventdata, handles)
% if handles.show_annotation == 0 % this click does not change anything
%     set(handles.radiobutton_annotation_no_show,'value',1);
%     return
% end
% handles.show_annotation = 0; % 0 = no show; 1 = show all; 2 = show selected
% guidata(hObject, handles);
% axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
% 
% 
% % --- Executes on button press in radiobutton_annotation_show_all.
% function radiobutton_annotation_show_all_Callback(hObject, eventdata, handles)
% if handles.show_annotation == 1 % this click does not change anything
%     set(handles.radiobutton_annotation_show_all,'value',1);
%     return
% end
% handles.show_annotation = 1; % 0 = no show; 1 = show all; 2 = show selected
% %???NOTE : compute bubbles if needed
% guidata(hObject, handles);
% axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
% 
% 
% % --- Executes on button press in radiobutton_annotation_show_selected.
% function radiobutton_annotation_show_selected_Callback(hObject, eventdata, handles)
% if handles.show_annotation == 2 % this click does not change anything
%     set(handles.radiobutton_annotation_show_selected,'value',1);
%     return
% end
% handles.show_annotation = 2; % 0 = no show; 1 = show all; 2 = show selected
% %???NOTE : compute bubbles if needed
% guidata(hObject, handles);
% axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);



% --- Executes on button press in checkbox_is_use_color.
function checkbox_is_use_color_Callback(hObject, eventdata, handles)
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);




% --- Executes on selection change in listbox_select_a_marker.
function listbox_select_a_marker_Callback(hObject, eventdata, handles)
if get(handles.checkbox_is_use_color,'value')==1 && handles.color_definition~=2
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end



% --- Executes during object creation, after setting all properties.
function listbox_select_a_marker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_select_a_marker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_select_a_file.
function listbox_select_a_file_Callback(hObject, eventdata, handles)
if get(handles.checkbox_is_use_color,'value')==1 
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end



% --- Executes during object creation, after setting all properties.
function listbox_select_a_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_select_a_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in listbox_select_ref_files.
function listbox_select_ref_files_Callback(hObject, eventdata, handles)
if get(handles.checkbox_is_use_color,'value')==1 
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end


% --- Executes during object creation, after setting all properties.
function listbox_select_ref_files_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_select_ref_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_add_ref_file.
function button_add_ref_file_Callback(hObject, eventdata, handles)
all_existing_files = get(handles.listbox_select_a_file, 'string');
file_to_be_added = all_existing_files(get(handles.listbox_select_a_file, 'value'));
if ismember(file_to_be_added, get(handles.listbox_select_ref_files,'string'))
    return
else
    new_ref_files = [get(handles.listbox_select_ref_files,'string');file_to_be_added];
    [C,IA,IB] = intersect(all_existing_files,new_ref_files);
    set(handles.listbox_select_ref_files,'string',all_existing_files(sort(IA)),'value',1);
    set(handles.button_remove_from_ref,'enable','on');
    if get(handles.checkbox_is_use_color,'value')==1 && handles.color_definition==1
        axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
    end
end


% --- Executes on button press in button_remove_from_ref.
function button_remove_from_ref_Callback(hObject, eventdata, handles)
ref_files = get(handles.listbox_select_ref_files,'string');
selected_ind = get(handles.listbox_select_ref_files,'value');
ref_files(selected_ind)=[];
set(handles.listbox_select_ref_files,'string',ref_files);
if selected_ind>length(ref_files)
    set(handles.listbox_select_ref_files,'value',selected_ind-1);
end
if length(ref_files)==1
    set(handles.button_remove_from_ref,'enable','off');
end
if get(handles.checkbox_is_use_color,'value')==1 && handles.color_definition==1
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end



% --- Executes on button press in radiobutton_color_definition_expr.
function radiobutton_color_definition_expr_Callback(hObject, eventdata, handles)
if handles.color_definition == 0 % this click does not change anything
    set(handles.radiobutton_color_definition_expr,'value',1);
    return
end
handles.color_definition = 0; % 0 = expr; 1 = ratio; 2 = cell freq
guidata(hObject, handles);
if get(handles.checkbox_is_use_color,'value')==1
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end


% --- Executes on button press in radiobutton_color_definition_ratio.
function radiobutton_color_definition_ratio_Callback(hObject, eventdata, handles)
if handles.color_definition == 1 % this click does not change anything
    set(handles.radiobutton_color_definition_ratio,'value',1);
    return
end
handles.color_definition = 1; % 0 = expr; 1 = ratio; 2 = cell freq
guidata(hObject, handles);
if get(handles.checkbox_is_use_color,'value')==1
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end


% --- Executes on button press in radiobutton_color_definition_cell_freq.
function radiobutton_color_definition_cell_freq_Callback(hObject, eventdata, handles)
if handles.color_definition == 2 % this click does not change anything
    set(handles.radiobutton_color_definition_cell_freq,'value',1);
    return
end
handles.color_definition = 2; % 0 = expr; 1 = ratio; 2 = cell freq
guidata(hObject, handles);
if get(handles.checkbox_is_use_color,'value')==1
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end



% --- Executes on button press in radiobutton_color_scheme_JET.
function radiobutton_color_scheme_JET_Callback(hObject, eventdata, handles)
if handles.color_scheme == 0 % this click does not change anything
    set(handles.radiobutton_color_scheme_JET,'value',1);
    return
end
handles.color_scheme = 0;     % 0 = JET; 1 = half JET; 2 = gray scale
guidata(hObject, handles);
if get(handles.checkbox_is_use_color,'value')==1
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end



% --- Executes on button press in radiobutton_color_scheme_half_JET.
function radiobutton_color_scheme_half_JET_Callback(hObject, eventdata, handles)
if handles.color_scheme == 1 % this click does not change anything
    set(handles.radiobutton_color_scheme_half_JET,'value',1);
    return
end
handles.color_scheme = 1;     % 0 = JET; 1 = half JET; 2 = gray scale
guidata(hObject, handles);
if get(handles.checkbox_is_use_color,'value')==1
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end



% --- Executes on button press in radiobutton_color_scheme_gray.
function radiobutton_color_scheme_gray_Callback(hObject, eventdata, handles)
if handles.color_scheme == 2 % this click does not change anything
    set(handles.radiobutton_color_scheme_gray,'value',1);
    return
end
handles.color_scheme = 2;     % 0 = JET; 1 = half JET; 2 = gray scale
guidata(hObject, handles);
if get(handles.checkbox_is_use_color,'value')==1
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end



function edit_mouse_selected_nodes_Callback(hObject, eventdata, handles)
user_input_index = str2num(get(handles.edit_mouse_selected_nodes,'string'));
if isempty(user_input_index) % user put invalid node index, we will ignore that and reset
    set(handles.edit_mouse_selected_nodes,'string',num2str(handles.mouse_selected_nodes));
    return
end
user_input_index = round(unique(user_input_index)); % in case there are non-integers and duplicates
user_input_index(user_input_index<1)=[];
user_input_index(user_input_index>length(handles.node_size))=[];
if isequal(sort(user_input_index),sort(handles.mouse_selected_nodes))  % if this event is invoked but the user did not actually update anything here
    set(handles.edit_mouse_selected_nodes,'string',num2str(handles.mouse_selected_nodes));
    return
else
    handles.mouse_selected_nodes = user_input_index;
    set(handles.edit_mouse_selected_nodes,'string',num2str(handles.mouse_selected_nodes));
    guidata(hObject, handles);
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end




% button down invoked by clicking a node
function edit_tree_NodeButtonDownFcn(hObject, eventdata, handles)
[dummy,node_num] = min(sum(abs(handles.node_positions-repmat([get(hObject,'XData');get(hObject,'YData')],1,size(handles.node_positions,2)))));
handles.clicked_node = node_num;
node_num
selectiontype = get(handles.figure1,'selectiontype');
if isequal(selectiontype,'alt') || isequal(selectiontype,'extend') % shift or Ctrl is down
    handles.mouse_down_position = handles.node_positions(:,node_num)';
    handles.mouse_mode='selection';
    handles.GETRECT_H1 = line('Parent', handles.Axes_mst, ...
                              'XData', [0 0 0 0 0]+handles.mouse_down_position(1), ...
                              'YData', [0 0 0 0 0]+handles.mouse_down_position(2), ...
                              'Visible', 'off', ...
                              'Clipping', 'off', ...
                              'Color', 'k', ...
                              'LineStyle', '-');
    handles.GETRECT_H2 = line('Parent', handles.Axes_mst, ...
                              'XData', [0 0 0 0 0]+handles.mouse_down_position(1), ...
                              'YData', [0 0 0 0 0]+handles.mouse_down_position(2), ...
                              'Visible', 'off', ...
                              'Clipping', 'off', ...
                              'Color', 'w', ...
                              'LineStyle', ':');
    set(handles.figure1,'WindowButtonMotionFcn','View_Edit_SPADE_tree_annotation(''edit_tree_WindowButtonMotionFcn'',gcbo,[],guidata(gcbo))');
    set(handles.figure1,'WindowButtonUpFcn','View_Edit_SPADE_tree_annotation(''edit_tree_WindowButtonUpFcn'',gcbo,[],guidata(gcbo))');
    handles.is_mouse_down=1;
    handles.is_mouse_down_and_moved=0;
    guidata(handles.button_show_tree_new_window, handles);    
    return
end

if ~ismember(node_num,handles.mouse_selected_nodes) % ctrl and alt not used, clicked a node that was not selected
    handles.mouse_selected_nodes = node_num;
    tmp = get(handles.Axes_mst,'currentPoint'); handles.mouse_down_position = tmp(1,1:2);
    handles.mouse_mode='move';
    guidata(handles.button_show_tree_new_window, handles);
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
    set(handles.edit_mouse_selected_nodes,'string',num2str(handles.mouse_selected_nodes));    handles.mouse_mode='move';
    set(handles.figure1,'WindowButtonMotionFcn','View_Edit_SPADE_tree_annotation(''edit_tree_WindowButtonMotionFcn'',gcbo,[],guidata(gcbo))');
    set(handles.figure1,'WindowButtonUpFcn','View_Edit_SPADE_tree_annotation(''edit_tree_WindowButtonUpFcn'',gcbo,[],guidata(gcbo))');
    handles.is_mouse_down=1;
    handles.is_mouse_down_and_moved=0;
    guidata(handles.button_show_tree_new_window, handles);    
    return
end

if ismember(node_num,handles.mouse_selected_nodes) % ctrl and alt not used, clicked a node that was selected
    tmp = get(handles.Axes_mst,'currentPoint'); handles.mouse_down_position = tmp(1,1:2);
    handles.mouse_mode='move';
    set(handles.figure1,'WindowButtonMotionFcn','View_Edit_SPADE_tree_annotation(''edit_tree_WindowButtonMotionFcn'',gcbo,[],guidata(gcbo))');
    set(handles.figure1,'WindowButtonUpFcn','View_Edit_SPADE_tree_annotation(''edit_tree_WindowButtonUpFcn'',gcbo,[],guidata(gcbo))');
    handles.is_mouse_down=1;
    handles.is_mouse_down_and_moved=0;
    guidata(handles.button_show_tree_new_window, handles);    
    return
end






% button down invoked by clicking empty space in this axis
function Axes_mst_ButtonDownFcn(hObject, eventdata, handles)
% get the position of the point when clicked
tmp = get(handles.Axes_mst,'currentPoint');
handles.mouse_down_position = tmp(1,1:2);
% get selection type, indicating the keyboard status
selectiontype = get(handles.figure1,'selectiontype');
if isequal(selectiontype,'open') % double clicked just now
    return
end
if isequal(selectiontype,'normal') 
    handles.mouse_selected_nodes = [];
end
if isequal(selectiontype,'alt') || isequal(selectiontype,'extend') % shift or Ctrl is down
    1; % do nothing
end
handles.mouse_mode='selection';
handles.GETRECT_H1 = line('Parent', handles.Axes_mst, ...
                          'XData', [0 0 0 0 0], ...
                          'YData', [0 0 0 0 0], ...
                          'Visible', 'off', ...
                          'Clipping', 'off', ...
                          'Color', 'k', ...
                          'LineStyle', '-');
handles.GETRECT_H2 = line('Parent', handles.Axes_mst, ...
                          'XData', [0 0 0 0 0], ...
                          'YData', [0 0 0 0 0], ...
                          'Visible', 'off', ...
                          'Clipping', 'off', ...
                          'Color', 'w', ...
                          'LineStyle', ':');
set(handles.figure1,'WindowButtonMotionFcn','View_Edit_SPADE_tree_annotation(''edit_tree_WindowButtonMotionFcn'',gcbo,[],guidata(gcbo))');
set(handles.figure1,'WindowButtonUpFcn','View_Edit_SPADE_tree_annotation(''edit_tree_WindowButtonUpFcn'',gcbo,[],guidata(gcbo))');
handles.is_mouse_down=1;
handles.is_mouse_down_and_moved=0;
guidata(hObject, handles);




function edit_tree_WindowButtonMotionFcn(hObject, eventdata, handles)
if isequal(handles.mouse_mode,'selection')
    tmp = get(handles.Axes_mst,'currentPoint');
    current_position = tmp(1,1:2);
    starting_position = handles.mouse_down_position;
    x_min = min(starting_position(1),current_position(1));
    x_max = max(starting_position(1),current_position(1));
    y_min = min(starting_position(2),current_position(2));
    y_max = max(starting_position(2),current_position(2));
    xdata = [x_min,x_max,x_max,x_min,x_min];
    ydata = [y_min,y_min,y_max,y_max,y_min];
    set(handles.GETRECT_H1, 'XData', xdata, 'YData', ydata, 'visible', 'on');
    set(handles.GETRECT_H2, 'XData', xdata, 'YData', ydata, 'visible', 'on');
end
if isequal(handles.mouse_mode,'move')
    set(handles.figure1,'WindowButtonMotionFcn',[]);
    tmp = get(handles.Axes_mst,'currentPoint');
    current_position = tmp(1,1:2);
    starting_position = handles.mouse_down_position;
    handles.node_positions(:,handles.mouse_selected_nodes) = handles.node_positions(:,handles.mouse_selected_nodes) + repmat(current_position(:)-starting_position(:),1,length(handles.mouse_selected_nodes));
    handles.is_mouse_down_and_moved=1;
    handles.mouse_down_position = current_position;
    for k=1:length(handles.tree_bubble_contour)% clear all the bubble contours
        handles.tree_bubble_contour{k}=[];
    end
    set(handles.radiobutton_tree_annotation_no_show,'value',1); handles.show_annotation=0;
    guidata(hObject, handles);    
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
    set(handles.figure1,'WindowButtonMotionFcn','View_Edit_SPADE_tree_annotation(''edit_tree_WindowButtonMotionFcn'',gcbo,[],guidata(gcbo))');
end



function edit_tree_WindowButtonUpFcn(hObject, eventdata, handles)
set(handles.figure1,'WindowButtonUpFcn',[]);
set(handles.figure1,'WindowButtonMotionFcn',[]);
% get selection type, indicating the keyboard status
selectiontype = get(handles.figure1,'selectiontype');
% update selections if in selection mode
if isequal(handles.mouse_mode,'selection')
    set(handles.GETRECT_H1, 'visible', 'off');
    set(handles.GETRECT_H2, 'visible', 'off');
    xdata = get(handles.GETRECT_H1,'XData');
    ydata = get(handles.GETRECT_H1,'yData');
    x_min = xdata(1); x_max = xdata(2);
    y_min = ydata(1); y_max = ydata(3);
    node_in_box = find(handles.node_positions(1,:)>=x_min & handles.node_positions(1,:)<=x_max & handles.node_positions(2,:)>=y_min & handles.node_positions(2,:)<=y_max);
    if isequal(selectiontype,'normal') 
        handles.mouse_selected_nodes = node_in_box;
        guidata(hObject, handles);
        axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
        set(handles.edit_mouse_selected_nodes,'string',num2str(handles.mouse_selected_nodes));
    elseif isequal(selectiontype,'extend')
        handles.mouse_selected_nodes = union(handles.mouse_selected_nodes,node_in_box);
        guidata(hObject, handles);
        axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
        set(handles.edit_mouse_selected_nodes,'string',num2str(handles.mouse_selected_nodes));
    elseif isequal(selectiontype,'alt')
        handles.mouse_selected_nodes = union(setdiff(handles.mouse_selected_nodes,node_in_box),setdiff(node_in_box,handles.mouse_selected_nodes));
        guidata(hObject, handles);
        axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
        set(handles.edit_mouse_selected_nodes,'string',num2str(handles.mouse_selected_nodes));
    end
end
if isequal(handles.mouse_mode,'move')
    if handles.is_mouse_down_and_moved==0 % this means the user clicked one node that was selected, but his intention was to select this node and de-select the others
        handles.mouse_selected_nodes = handles.clicked_node;
        guidata(hObject, handles);
        axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
        set(handles.edit_mouse_selected_nodes,'string',num2str(handles.mouse_selected_nodes));
    end
end
handles.is_mouse_down=0;
handles.is_mouse_down_and_moved=0;
handles.mouse_mode=[];
guidata(hObject, handles);





% --- Executes on button press in button_add_selected_nodes_to_annotation.
function button_add_selected_nodes_to_annotation_Callback(hObject, eventdata, handles)
if isempty(handles.mouse_selected_nodes)
    return
end
handles.tree_annotations{end+1} = handles.mouse_selected_nodes;
handles.tree_bubble_contour{end+1}=cell(0);
% clear all contours because one new bubble is added
for k=1:length(handles.tree_bubble_contour) 
    handles.tree_bubble_contour{k}=[];
end
% update the listbox
tmp = [];
for i=1:length(handles.tree_annotations)
    tmp{i} = num2str(handles.tree_annotations{i});
end
set(handles.listbox_annotations,'string',tmp,'value',length(handles.tree_annotations));
% update the viz option, no_show
set(handles.radiobutton_tree_annotation_no_show,'value',1); handles.show_annotation=0;
guidata(hObject,handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);  



% --- Executes on button press in button_remove_annotation.
function button_remove_annotation_Callback(hObject, eventdata, handles)
if isempty(handles.tree_annotations)
    return
end
ind_to_remove = get(handles.listbox_annotations,'value');
handles.tree_annotations(ind_to_remove) = [];
handles.tree_bubble_contour(ind_to_remove) = [];
tmp = [];
for i=1:length(handles.tree_annotations)
    tmp{i} = num2str(handles.tree_annotations{i});
end
set(handles.listbox_annotations,'string',tmp);
if ind_to_remove>length(handles.tree_annotations)
    set(handles.listbox_annotations,'value',max(1,ind_to_remove-1));
end
guidata(hObject,handles);
listbox_annotations_Callback(handles.listbox_annotations, [], handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);



% --- Executes on button press in radiobutton_tree_annotation_no_show.
function radiobutton_tree_annotation_no_show_Callback(hObject, eventdata, handles)
if handles.show_annotation == 0 % this click does not change anything
    set(handles.radiobutton_tree_annotation_no_show,'value',1);
    return
end
handles.show_annotation = 0; % 0 = no show; 1 = show all; 2 = show selected
guidata(hObject, handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);




% --- Executes on button press in radiobutton_tree_annotation_show_all.
function radiobutton_tree_annotation_show_all_Callback(hObject, eventdata, handles)
if handles.show_annotation == 1 % this click does not change anything
    set(handles.radiobutton_tree_annotation_show_all,'value',1);
    return
end
handles.show_annotation = 1; % 0 = no show; 1 = show all; 2 = show selected
for k=1:length(handles.tree_bubble_contour)
    if ~isempty(handles.tree_bubble_contour{k})
        continue;
    else
        bubble_contours = compute_annotation_contour_here(handles);
        handles.tree_bubble_contour = bubble_contours;
        break; 
    end
end
guidata(hObject, handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);




% --- Executes on button press in radiobutton_tree_annotation_show_selected.
function radiobutton_tree_annotation_show_selected_Callback(hObject, eventdata, handles)
if handles.show_annotation == 2 % this click does not change anything
    set(handles.radiobutton_tree_annotation_show_selected,'value',1);
    return
end
handles.show_annotation = 2; % 0 = no show; 1 = show all; 2 = show selected
for k=1:length(handles.tree_bubble_contour)
    if ~isempty(handles.tree_bubble_contour{k})
        continue;
    else
        bubble_contours = compute_annotation_contour_here(handles);
        handles.tree_bubble_contour = bubble_contours;
        break; 
    end
end
guidata(hObject, handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);





function slider2_Callback(hObject, eventdata, handles)

function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider3_Callback(hObject, eventdata, handles)

function slider3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider4_Callback(hObject, eventdata, handles)

function slider4_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on slider movement.
function slider_scale_Callback(hObject, eventdata, handles)
if get(hObject,'Value')>0
    % enlarge
    if length(handles.mouse_selected_nodes)>2
        center_positions = mean(handles.node_positions(:,handles.mouse_selected_nodes),2);
        for i=1:length(handles.mouse_selected_nodes)
            handles.node_positions(:,handles.mouse_selected_nodes(i)) = center_positions + (handles.node_positions(:,handles.mouse_selected_nodes(i))-center_positions)*1.05;
        end
    end
else
    % shrink
    if length(handles.mouse_selected_nodes)>2
        center_positions = mean(handles.node_positions(:,handles.mouse_selected_nodes),2);
        for i=1:length(handles.mouse_selected_nodes)
            handles.node_positions(:,handles.mouse_selected_nodes(i)) = center_positions + (handles.node_positions(:,handles.mouse_selected_nodes(i))-center_positions)/1.05;
        end        
    end
end
set(hObject,'Value',0);
guidata(hObject,handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);



% --- Executes on slider movement.
function slider_rotate_Callback(hObject, eventdata, handles)
if get(hObject,'Value')>0
    % clock-wise
    if length(handles.mouse_selected_nodes)>2
        r = [cos(pi/18), sin(pi/18);-sin(pi/18), cos(pi/18)];
        center_positions = mean(handles.node_positions(:,handles.mouse_selected_nodes),2);
        for i=1:length(handles.mouse_selected_nodes)
            handles.node_positions(:,handles.mouse_selected_nodes(i)) = center_positions + r*(handles.node_positions(:,handles.mouse_selected_nodes(i))-center_positions);
        end
    end
else
    % counter clock-wise
    if length(handles.mouse_selected_nodes)>2
        r = [cos(pi/18), -sin(pi/18);+sin(pi/18), cos(pi/18)];
        center_positions = mean(handles.node_positions(:,handles.mouse_selected_nodes),2);
        for i=1:length(handles.mouse_selected_nodes)
            handles.node_positions(:,handles.mouse_selected_nodes(i)) = center_positions + r*(handles.node_positions(:,handles.mouse_selected_nodes(i))-center_positions);
        end        
    end
end
set(hObject,'Value',0);
guidata(hObject,handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);



% --- Executes on slider movement.
function slider_change_node_size_Callback(hObject, eventdata, handles)
if get(hObject,'Value')>0
    handles.node_size = handles.node_size*1.1;
else
    handles.node_size = handles.node_size*0.9;
end
set(hObject,'Value',0);
guidata(hObject,handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

node_positions = handles.node_positions;
node_size = handles.node_size;
tree_annotations = handles.tree_annotations;
tree_bubble_contour = handles.tree_bubble_contour;
user_defined_color_thresholds = handles.user_defined_color_thresholds;
save(handles.cluster_mst_upsample_filename,'node_positions','node_size','tree_annotations','tree_bubble_contour','user_defined_color_thresholds','-append');
delete(hObject);


% --- Executes on button press in button_export_selected.
function button_export_selected_Callback(hObject, eventdata, handles)
% hObject    handle to button_export_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_update_plots_Callback(handles.button_update_plots, [], handles);
if isempty(handles.mouse_selected_nodes)
    return
end
[PATHSTR,NAME,EXT]  = fileparts(handles.cluster_mst_upsample_filename);
if exist(fullfile(PATHSTR,'exported_results'))~=7
    mkdir(fullfile(PATHSTR,'exported_results'));
end
for i=1:length(handles.all_existing_files)
    if isequal(handles.all_existing_files{i},'POOLED')
        idx_tmp = handles.clustering_idx;
        data_tmp = handles.clustering_data;
        marker_names_tmp = handles.clustering_marker_names;
    else
        idx_tmp = handles.all_assign{find(ismember(handles.file_annot,handles.all_existing_files(i))==1)};
        [data_tmp, marker_names_tmp] = readfcs(fullfile(PATHSTR,handles.all_fcs_filenames{find(ismember(handles.file_annot,handles.all_existing_files(i))==1)}));
    end
    ind = zeros(1,length(idx_tmp));
    for k=1:length(handles.mouse_selected_nodes)
        ind(idx_tmp==handles.mouse_selected_nodes(k))=1;
    end
    data_tmp = data_tmp(:,ind==1);
    % fcs write data marker_names filename
    if exist(fullfile(PATHSTR,'exported_results','selected_cells_in_fcs'))~=7
        mkdir(fullfile(PATHSTR,'exported_results','selected_cells_in_fcs'))
    end
    writefcs(fullfile(PATHSTR,'exported_results','selected_cells_in_fcs',[handles.all_existing_files{i},'_selected_nodes.fcs']), data_tmp, marker_names_tmp,marker_names_tmp);
end



function edit_annotation_bubble_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_annotation_bubble_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_value = str2num(get(handles.edit_annotation_bubble_size,'string')); 
if isempty(new_value)
    set(handles.edit_annotation_bubble_size,'string',num2str(get(handles.edit_annotation_bubble_size,'value')))
    return
end
set(handles.edit_annotation_bubble_size,'value',new_value);
for k=1:length(handles.tree_bubble_contour)
    handles.tree_bubble_contour{k}=[];
    guidata(hObject, handles);
end
if handles.show_annotation~=0
    bubble_contours = compute_annotation_contour_here(handles);
    handles.tree_bubble_contour = bubble_contours;
    guidata(hObject, handles);
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end


% --- Executes during object creation, after setting all properties.
function edit_annotation_bubble_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_annotation_bubble_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%%%%%%%%%%%%%%%%%%% below is the function to calculate the annotationa bubble
function bubble_contours = compute_annotation_contour_here(handles) 

mst_tree = handles.mst_tree;
node_positions = handles.node_positions;
node_size = handles.node_size;
bubbles = handles.tree_annotations;
r = ones(1,length(bubbles))*str2num(get(handles.edit_annotation_bubble_size,'string')); 
resolution = 100;

axes(handles.Axes_mst); canvas_axis = axis; 
dist = zeros(length(canvas_axis(1):(canvas_axis(2)-canvas_axis(1))/resolution:canvas_axis(2))*length(canvas_axis(3):(canvas_axis(4)-canvas_axis(3))/resolution:canvas_axis(4)),size(node_positions,2));
grid_position = zeros(size(dist,1),2);
delta_x = (canvas_axis(2)-canvas_axis(1))/resolution;
delta_y = (canvas_axis(4)-canvas_axis(3))/resolution;
counter = 1;
for x=canvas_axis(1):delta_x:canvas_axis(2)
    for y=canvas_axis(3):delta_y:canvas_axis(4)
        dist(counter,:) = sqrt((node_positions(1,:)-x).^2 + (node_positions(2,:)-y).^2);
        grid_position(counter,:) = [x,y];
        counter = counter + 1;
    end
end


for i=1:length(bubbles)
    indicators = (sum(dist(:,bubbles{i})<r(i),2)~=0) & (min(dist(:,bubbles{i})')'+r(i)/10<min(dist(:,setdiff(1:end,bubbles{i}))')');
    % it is possible that there are holes in the region labeled by indicators, we want to fill those holes so that the bubbles look nice.   
    tmp = reshape(double(indicators),sqrt(length(indicators)),sqrt(length(indicators)));
    tmp(1,tmp(1,:)==0)=-1;
    tmp(end,tmp(end,:)==0)=-1;
    tmp(tmp(:,1)==0,1)=-1;
    tmp(tmp(:,end)==0,end)=-1;
    old_tmp = tmp;
    while 1
        tmp([-ones(1,size(tmp,2));tmp(1:end-1,:)]==-1 & tmp==0)=-1;
        tmp([tmp(2:end,:);-ones(1,size(tmp,2))]==-1 & tmp==0)=-1;
        tmp([-ones(size(tmp,1),1),tmp(:,1:end-1)]==-1 & tmp==0)=-1;
        tmp([tmp(:,2:end),-ones(size(tmp,1),1)]==-1 & tmp==0)=-1;
        if isequal(old_tmp,tmp)
            break;
        end
        old_tmp=tmp;
        % figure(1);subplot(2,2,4);imagesc(flipud(tmp))
    end
    indicators(tmp(:)==0)=1;
    
    outer_points = get_outer_points(grid_position, indicators, delta_x, delta_y);
%     line([outer_points(:,1); outer_points(1,1)],[outer_points(:,2);outer_points(1,2)])
%     for i=1:size(outer_points,1), plot(outer_points(i,1),outer_points(i,2),'g.'); pause(0.05); end
%     for i=1:size(outer_points,1)-1, line(outer_points(i:i+1,1),outer_points(i:i+1,2)); pause(0.05); end
    bubble_contours{i} = outer_points;
end



function outer_points = get_outer_points(grid_position, indicators, delta_x, delta_y)

patch_squareform = reshape(indicators,sqrt(length(indicators)), sqrt(length(indicators)));
contour_squareform = zeros(size(patch_squareform));

for i=1:size(patch_squareform,1)
    ind_tmp = (patch_squareform(i,:)~=[-1,patch_squareform(i,1:end-1)] | patch_squareform(i,:)~=[patch_squareform(i,2:end),-1]) & patch_squareform(i,:)==1;
    contour_squareform(i,ind_tmp)=1;
end
for j=1:size(patch_squareform,2)
    ind_tmp = (patch_squareform(:,j)~=[-1;patch_squareform(1:end-1,j)] | patch_squareform(:,j)~=[patch_squareform(2:end,j);-1]) & patch_squareform(:,j)==1;
    contour_squareform(ind_tmp,j)=1;
end

contour_squareform = contour_squareform(:);
outer_points = grid_position(contour_squareform==1,:);

% order them
ordered_points = outer_points(1,:); 
flag_outer_points = zeros(size(outer_points,1),1);flag_outer_points(1)=1;
while sum(flag_outer_points)~=length(flag_outer_points)
    dist = sum((repmat(ordered_points(end,:),size(outer_points,1),1) - outer_points).^2,2);
    dist(flag_outer_points==1) = max(dist) + 1;
    [dummy,I] = min(dist);
%     if dummy > sqrt(delta_x^2+delta_y^2)*2;
%         flag_outer_points(I)=1; continue;
%     end
    ordered_points = [ordered_points; outer_points(I,:)]; flag_outer_points(I)=1;
end
outer_points = ordered_points;


% --- Executes on button press in pushbutton_Tree_figure_to_ps.
function pushbutton_Tree_figure_to_ps_Callback(hObject, eventdata, handles)
[PATHSTR,NAME,EXT]  = fileparts(handles.cluster_mst_upsample_filename);
if exist(fullfile(PATHSTR,'exported_results'))~=7
    mkdir(fullfile(PATHSTR,'exported_results'));
end
% write the plot without color
set(handles.checkbox_is_use_color,'value',0)
set(handles.radiobutton_color_definition_expr,'value',1); handles.color_definition=0; % 0 = expr; 1 = ratio; 2 = cell freq
% set(handles.radiobutton_color_scheme_JET,'value',1)
figure(1); close(1);
h = figure(1); draw_SPADE_tree_no_mouse_click_events(handles); axis off; 
title('SPADE tree'); drawnow; 
print(h,fullfile(PATHSTR,'exported_results','SPADE_tree_colored_by_markers.ps'),'-dpsc2');
% write colored versions
set(handles.checkbox_is_use_color,'value',1);
markers_to_show = get(handles.listbox_select_a_marker,'string');
for k=1:length(markers_to_show)-2
    set(handles.listbox_select_a_marker,'value',k)
    set(handles.listbox_select_a_file,'value',1)
    h = figure(1); draw_SPADE_tree_no_mouse_click_events(handles); axis off; 
    title(markers_to_show{k},'interpret','none'); drawnow; 
    print(h,fullfile(PATHSTR,'exported_results','SPADE_tree_colored_by_markers.ps'),'-dpsc2','-append');
end
close(1);
% ps2pdf('psfile', fullfile(PATHSTR,'exported_results','SPADE_tree_colored_by_markers.ps'), 'pdffile', fullfile(PATHSTR,'exported_results','SPADE_tree_colored_by_markers.pdf'), 'gspapersize', 'a4');
% refresh the output panel 
set(handles.checkbox_is_use_color,'value',0)
set(handles.radiobutton_color_definition_expr,'value',1)
set(handles.radiobutton_color_scheme_JET,'value',1)
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);


% --- Executes on button press in pushbutton_Tree_to_GML_TXT.
function pushbutton_Tree_to_GML_TXT_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Tree_to_GML_TXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[PATHSTR,NAME,EXT]  = fileparts(handles.cluster_mst_upsample_filename);
if exist(fullfile(PATHSTR,'exported_results'))~=7
    mkdir(fullfile(PATHSTR,'exported_results'));
end
node_names = cell(size(handles.mst_tree,1),1); for k=1:length(node_names), node_names{k} = num2str(k); end
% write to GML
write_undirected_adj_position_to_gml(fullfile(PATHSTR,'exported_results','SPADE_tree.gml'), handles.mst_tree, handles.node_positions, node_names);
% write node positions to txt
write_to_txt(fullfile(PATHSTR,'exported_results','SPADE_tree_node_position.txt'), [{'node'},{'x position'},{'y position'}], node_names, handles.node_positions', char(9))
% write node adj to txt
write_to_txt(fullfile(PATHSTR,'exported_results','SPADE_tree_adjacency.txt'), [{'node'},node_names'], node_names, full(handles.mst_tree), char(9))




% --- Executes on button press in pushbutton_node_annot_expr_cell_freq_to_TXT.
function pushbutton_node_annot_expr_cell_freq_to_TXT_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_node_annot_expr_cell_freq_to_TXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[PATHSTR,NAME,EXT]  = fileparts(handles.cluster_mst_upsample_filename);
if exist(fullfile(PATHSTR,'exported_results'))~=7
    mkdir(fullfile(PATHSTR,'exported_results'));
end
% define node names
node_names = cell(size(handles.mst_tree,1),1); for k=1:length(node_names), node_names{k} = num2str(k); end
% standardize tree_annotation_definitions
tree_annotation_definitions = get(handles.listbox_annotations,'string');
for k=1:length(tree_annotation_definitions)
    elements_in_annotation = str2num(tree_annotation_definitions{k});
    tree_annotation_definitions{k} = [];
    for m = 1:length(elements_in_annotation)
        tree_annotation_definitions{k} = [tree_annotation_definitions{k}, num2str(elements_in_annotation(m)), ' '];
    end
    tree_annotation_definitions{k} = tree_annotation_definitions{k}(1:end-1);
end
% define bubble names
tree_annotation_names = cell(length(tree_annotation_definitions),1); for k=1:length(tree_annotation_names), tree_annotation_names{k} = ['annotation-',num2str(k)]; end
% write bubble definition
write_to_txt(fullfile(PATHSTR,'exported_results','SPADE_tree_annotation_definition.txt'), [], [tree_annotation_names, tree_annotation_definitions], [], char(9))

% write node/bubble marker expr into txt according to POOLED and overlapping markers
tmp = handles.clustering_marker_node_average(find(ismember(handles.clustering_marker_node_average(:,1),'POOLED') & ~ismember(handles.clustering_marker_node_average(:,2),[{'FileInd'},{'CellFreq'}])),:);
d = []; 
for k=1:size(tmp,1), 
    d = [d; tmp{k,3}]; 
end; 
for k=1:length(tree_annotation_definitions)
    elements_in_annotation = str2num(tree_annotation_definitions{k});
    d = [d, sum(d(:,elements_in_annotation),2)/length(elements_in_annotation)];
end
d = d';
write_to_txt(fullfile(PATHSTR,'exported_results','SPADE_tree_node_annot_marker_expr.txt'), [{'node/annot'}, tmp(:,2)'] , [node_names;tree_annotation_names], d, char(9));

% file specific marker expr for each node and annotation
if exist(fullfile(PATHSTR,'exported_results','file_specific_marker_expr'))~=7
    mkdir(fullfile(PATHSTR,'exported_results','file_specific_marker_expr'));
end
all_possible_marker_names = setdiff(unique(handles.clustering_marker_node_average(:,2)),{'CellFreq','FileInd'}');
for j = 1:length(all_possible_marker_names)
    tmp = handles.clustering_marker_node_average(find(~ismember(handles.clustering_marker_node_average(:,1),'POOLED') & ismember(handles.clustering_marker_node_average(:,2),all_possible_marker_names(j))),:);
    d = []; 
    for k=1:size(tmp,1), 
        d = [d; tmp{k,3}]; 
    end 
    for k=1:length(tree_annotation_definitions)
        elements_in_annotation = str2num(tree_annotation_definitions{k});
        d = [d, nanmean(d(:,elements_in_annotation),2)];
    end
    d = d';
    write_to_txt(fullfile(PATHSTR,'exported_results','file_specific_marker_expr',[all_possible_marker_names{j},'_expr_file_specific.txt']), [{'node/annot'}, tmp(:,1)'] , [node_names;tree_annotation_names], d, char(9));
end


% write node/bubble cell freq according to individual files 
tmp = handles.clustering_marker_node_average(find(~ismember(handles.clustering_marker_node_average(:,1),'POOLED') & ismember(handles.clustering_marker_node_average(:,2),[{'CellFreq'}])),:);
d = []; 
for k=1:size(tmp,1), 
    d = [d; tmp{k,3};tmp{k,3}/sum(tmp{k,3})*100]; 
end; 
for k=1:length(tree_annotation_definitions)
    elements_in_annotation = str2num(tree_annotation_definitions{k});
    d = [d, sum(d(:,elements_in_annotation),2)];
end
d = d';
write_to_txt(fullfile(PATHSTR,'exported_results','SPADE_tree_node_annot_cell_freq.txt'), [{'node/annot'}, reshape([tmp(:,1),strcat(tmp(:,1),' - %')]',1,size(tmp,1)*2)] , [node_names;tree_annotation_names], d, char(9));





% --- Executes on button press in pushbutton_Clustering_result_to_FCS.
function pushbutton_Clustering_result_to_FCS_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Clustering_result_to_FCS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf('\nSaving clustering results to fcs ...\n')
[PATHSTR,NAME,EXT]  = fileparts(handles.cluster_mst_upsample_filename);
if exist(fullfile(PATHSTR,'exported_results'))~=7
    mkdir(fullfile(PATHSTR,'exported_results'));
end
for i=1:length(handles.all_existing_files)
    if isequal(handles.all_existing_files{i},'POOLED')
        continue;
    else
        idx_tmp = handles.all_assign{find(ismember(handles.file_annot,handles.all_existing_files(i))==1)};
        [data_tmp, marker_names_tmp, channel_names_tmp] = readfcs(fullfile(PATHSTR,handles.all_fcs_filenames{find(ismember(handles.file_annot,handles.all_existing_files(i))==1)}));
    end
    data_tmp = [data_tmp(:,1:length(idx_tmp));idx_tmp];
    marker_names_tmp{end+1} = 'cluster-assign';
    channel_names_tmp{end+1} = 'cluster-assign';
    if ~isempty(handles.tree_annotations)
        annotation_tmp = zeros(1,size(data_tmp,2));
        for k=1:length(handles.tree_annotations)
            annotation_tmp = annotation_tmp + double(ismember(idx_tmp,handles.tree_annotations{k}))*k;
        end
        data_tmp = [data_tmp;annotation_tmp];
        marker_names_tmp{end+1} = 'annotation';
        channel_names_tmp{end+1} = 'annotation';
    end
    % fcs write data marker_names filename
    if exist(fullfile(PATHSTR,'exported_results','clustering_result_added_in_fcs'))~=7
        mkdir(fullfile(PATHSTR,'exported_results','clustering_result_added_in_fcs'))
    end
    writefcs(fullfile(PATHSTR,'exported_results','clustering_result_added_in_fcs',[handles.all_fcs_filenames{find(ismember(handles.file_annot,handles.all_existing_files(i))==1)}]), data_tmp, marker_names_tmp,channel_names_tmp);
end
fprintf('\nDone!\n')




% --- Executes on button press in pushbutton_Tree_CellFreq_to_ps.
function pushbutton_Tree_CellFreq_to_ps_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Tree_CellFreq_to_ps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[PATHSTR,NAME,EXT]  = fileparts(handles.cluster_mst_upsample_filename);
if exist(fullfile(PATHSTR,'exported_results'))~=7
    mkdir(fullfile(PATHSTR,'exported_results'));
end
% write the plot without color
set(handles.checkbox_is_use_color,'value',0)
set(handles.radiobutton_color_definition_expr,'value',1); ; handles.color_definition=0; % 0 = expr; 1 = ratio; 2 = cell freq
% set(handles.radiobutton_color_scheme_JET,'value',1)
figure(1); close(1);
h = figure(1); draw_SPADE_tree_no_mouse_click_events(handles); axis off; 
title('SPADE tree'); drawnow; 
print(h,fullfile(PATHSTR,'exported_results','SPADE_tree_colored_by_CellFreq.ps'),'-dpsc2');
% write colored versions
set(handles.checkbox_is_use_color,'value',1);
markers_to_show = get(handles.listbox_select_a_marker,'string');
files_to_show = get(handles.listbox_select_a_file,'string');
for k=2:length(files_to_show)
    set(handles.listbox_select_a_marker,'value',length(markers_to_show));
    set(handles.listbox_select_a_file,'value',k)
    h = figure(1); draw_SPADE_tree_no_mouse_click_events(handles); axis off; 
    title([files_to_show{k},' - Cell Freq'],'interpret','none'); drawnow; 
    print(h,fullfile(PATHSTR,'exported_results','SPADE_tree_colored_by_CellFreq.ps'),'-dpsc2','-append');
end
close(1);
% ps2pdf('psfile', fullfile(PATHSTR,'exported_results','SPADE_tree_colored_by_CellFreq.ps'), 'pdffile', fullfile(PATHSTR,'exported_results','SPADE_tree_colored_by_CellFreq.pdf'), 'gspapersize', 'a4');
% refresh the output panel 
set(handles.checkbox_is_use_color,'value',0)
set(handles.radiobutton_color_definition_expr,'value',1)
set(handles.radiobutton_color_scheme_JET,'value',1)
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);


% --- Executes on button press in pushbutton_EMD.
function pushbutton_EMD_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_EMD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[PATHSTR,NAME,EXT]  = fileparts(handles.cluster_mst_upsample_filename);
if exist(fullfile(PATHSTR,'exported_results'))~=7
    mkdir(fullfile(PATHSTR,'exported_results'));
end
number_of_files = length(handles.all_fcs_filenames);
all_fcs_filenames = handles.all_fcs_filenames;

if length(all_fcs_filenames)==1
    fprintf('Pairwise EMD not calculated because you only have one sample!!!\n\n');
    return
end

if exist(fullfile(PATHSTR,'exported_results','earthmoverdist.mat'))~=2
    ind = find(~ismember(handles.clustering_marker_node_average(:,1),'POOLED') & ismember(handles.clustering_marker_node_average(:,2),'CellFreq'));
    file_annot = handles.clustering_marker_node_average(ind,1);
    CellFreq = [];
    for i=1:length(ind)
        CellFreq = [CellFreq; handles.clustering_marker_node_average{ind(i),3}/sum(handles.clustering_marker_node_average{ind(i),3})];
    end
    [is_connected, shortest_hop] = check_graph_connectivity(handles.mst_tree);
    EMD = zeros(number_of_files,number_of_files);
    fprintf('Computing pairwise EMDs for %d samples: %5d,%5d',number_of_files,0,0);
    for i=1:number_of_files-1
        for j=i+1:number_of_files
            fprintf('\b\b\b\b\b\b\b\b\b\b\b%5d,%5d',i,j);
            EMD(i,j) = EarthMoverDist2(CellFreq(i,:),CellFreq(j,:),shortest_hop); 
            drawnow;
        end
    end
    EMD = EMD + EMD';
    fprintf('\nDone!!\n\n');    
    save(fullfile(PATHSTR,'exported_results','earthmoverdist.mat'),'EMD','all_fcs_filenames','file_annot');
    write_to_txt(fullfile(PATHSTR,'exported_results','earthmoverdist.txt'), [{'Samples'},file_annot'], file_annot, EMD, char(9));
else
    load(fullfile(PATHSTR,'exported_results','earthmoverdist.mat'),'EMD','all_fcs_filenames','file_annot');
end

smallest_value_in_upper_triu = min(min(tril(ones(size(EMD)),0)*max(max(EMD)) + EMD));
EMD_display = EMD + eye(size(EMD))*smallest_value_in_upper_triu;
h = figure(1);
[perm_r,perm_c] = HCC_heatmap(EMD_display);
imagesc(EMD_display(perm_r,perm_r));
axis off
for i=1:length(perm_r)
    text(i,length(perm_r)+1,file_annot{perm_r(i)},'fontsize',12,'Rotation',-75,'interpret','none')
    text(-length(perm_r)/7,i,file_annot{perm_r(i)},'fontsize',12,'interpret','none')
end
print(h,fullfile(PATHSTR,'exported_results','earthmoverdist_figure.ps'),'-dpsc2');



% --- Executes during object creation, after setting all properties.
function Axes_mst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Axes_mst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate Axes_mst


% --- Executes on button press in pushbutton_reset_layout.
function pushbutton_reset_layout_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset_layout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.node_positions = arch_layout(handles.mst_tree);
handles.tree_bubble_contour = repmat({[]},1,length(handles.tree_annotations));
guidata(hObject,handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);


% --- Executes on button press in pushbutton_TPE_layout.
function pushbutton_TPE_layout_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_TPE_layout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[C,IA,IB] = intersect(handles.clustering_marker_names, handles.clustering_used_markers);
centers = grpstats(handles.clustering_data(IA,:)',handles.clustering_idx)';
weights = squareform(pdist(centers')); 
dist = graphallshortestpaths(weights.*handles.mst_tree);

node_positions = TPE_2D_contour_only(dist, handles.node_positions);
node_positions = node_positions - repmat(median(node_positions,2),1,size(node_positions,2));
node_positions(1,:) = node_positions(1,:)/max(abs(node_positions(:)))*50;
node_positions(2,:) = node_positions(2,:)/max(abs(node_positions(:)))*50;
handles.node_positions = node_positions;
guidata(hObject,handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);


% --- Executes on button press in pushbutton_Spring_layout.
function pushbutton_Spring_layout_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Spring_layout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

A = handles.mst_tree;

% [idx,centers,centers_ind]=deterministic_kmeans_for_SPADE(handles.clustering_data(ismember(handles.clustering_marker_names, handles.clustering_used_markers),:)',max(handles.clustering_idx),handles.clustering_local_density);
% centers = centers';
% [A] = SPADE_center_density_to_tree(idx,centers,centers_ind, handles.clustering_data(ismember(handles.clustering_marker_names, handles.clustering_used_markers),:), handles.clustering_local_density);
% [center_dist_scaled, center_dist, scaling_factor] = SPADE_center_pairwise_dist(idx,centers,centers_ind, handles.clustering_data(ismember(handles.clustering_marker_names, handles.clustering_used_markers),:), handles.clustering_local_density);
% [A, adj2, cost_value] = mst_from_dist_matrix(center_dist_scaled);
% k=3;
% for i=1:size(A,1)
%     [~,ind] = sort(center_dist_scaled(i,:),'ascend');
%     for j=2:k+1
%         if scaling_factor(i,ind(j))<=1.2
%             A(i,ind(j))=1;
%             A(ind(j),i)=1;
%         end
%     end
% end



[X,spring,distance]=kamada_kawai_spring_layout_mex(...
    A, 1e-30, 20000, 1, ...   % adj, tolerance, max iteration, spring_constant
    [], 1, 0, 'matrix');     % progressive_opt, options.edge_length, edge_weights, edge_weight_opt;

% [X,spring,distance]=kamada_kawai_spring_layout_mex(...
%     squareform(pdist(centers')).*A, 1e-30, 2000, 1, ...   % adj, tolerance, max iteration, spring_constant
%     [], 1, 0, 'matrix');     % progressive_opt, options.edge_length, edge_weights, edge_weight_opt;


% use spring embedding to determine node positoin
node_positions = X';
% normalize node positions
node_positions = node_positions - repmat((max(node_positions,[],2)+min(node_positions,[],2))/2,1,size(node_positions,2));
node_positions = node_positions/max(abs(node_positions(:)))*50;
% rotate so that the highest density point is in the west
node_local_density = accumarray(handles.clustering_idx', handles.clustering_local_density')';
% weight_center = sum(node_positions.*repmat(node_local_density,2,1),2)/sum(node_local_density);
weight_center = find_highest_density_position(node_positions, node_local_density);
tmp_score = zeros(1,360) + Inf;
for i=1:360
    tmp_angle = i/180*pi;
    tmp = [cos(tmp_angle), sin(tmp_angle); -sin(tmp_angle), cos(tmp_angle)]*weight_center(:);
    if tmp(1)<0 
        tmp_score(i)=abs(tmp(2));
    end
end
[~,i] = min(tmp_score);
tmp_angle = i/180*pi;
node_positions = [cos(tmp_angle), sin(tmp_angle); -sin(tmp_angle), cos(tmp_angle)]*node_positions;
% flipup to make sure that more density is in the north half
if sum(node_local_density(node_positions(1,:)>0 & node_positions(2,:)>0)) < sum(node_local_density(node_positions(1,:)>0 & node_positions(2,:)<0))
    node_positions(2,:)=-node_positions(2,:);
end


handles.node_positions = node_positions;
handles.tree_bubble_contour = repmat({[]},1,length(handles.tree_annotations));
guidata(hObject,handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);

% viSNE_color_plots(node_positions, ...
%                   cell2mat(handles.clustering_marker_node_average(ismember(handles.clustering_marker_node_average(:,1),'POOLED') & ismember(handles.clustering_marker_node_average(:,2),handles.clustering_used_markers) , 3)), ...
%                   handles.clustering_marker_node_average(ismember(handles.clustering_marker_node_average(:,1),'POOLED') & ismember(handles.clustering_marker_node_average(:,2),handles.clustering_used_markers), 2), ...
%                   node_local_density);
% h = figure;
% % write colored versions
% set(handles.checkbox_is_use_color,'value',1);
% markers_to_show = intersect(handles.clustering_used_markers, get(handles.listbox_select_a_marker,'string'));
% for k=1:length(markers_to_show)
%     set(handles.listbox_select_a_marker,'value',find(ismember(get(handles.listbox_select_a_marker,'string'),markers_to_show(k))));
%     set(handles.listbox_select_a_file,'value',1)
%     subplot(3,3,k); draw_SPADE_tree_no_mouse_click_events(handles); axis off; 
%     title(markers_to_show{k},'interpret','none'); drawnow; 
% end
% set(handles.checkbox_is_use_color,'value',0)
% set(handles.radiobutton_color_definition_expr,'value',1)
% set(handles.radiobutton_color_scheme_JET,'value',1)
% axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);



% --- Executes on button press in button_auto_suggest_annotation.
function button_auto_suggest_annotation_Callback(hObject, eventdata, handles)
% hObject    handle to button_auto_suggest_annotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.mouse_selected_nodes)
    considered_nodes = 1:size(handles.mst_tree,1);
else
    considered_nodes = handles.mouse_selected_nodes;
end
considered_mst_tree = handles.mst_tree;
for i=1:length(handles.tree_annotations)
    considered_mst_tree(handles.tree_annotations{i},setdiff(1:size(handles.mst_tree,1),handles.tree_annotations{i}))=0;
    considered_mst_tree(setdiff(1:size(handles.mst_tree,1),handles.tree_annotations{i}),handles.tree_annotations{i})=0;
end
for i=setdiff(1:size(handles.mst_tree,1),considered_nodes)
    considered_mst_tree(i,:)=0;
    considered_mst_tree(:,i)=0;
end

data = handles.clustering_data(ismember(handles.clustering_marker_names, handles.clustering_used_markers),:); % pooled downsampled data used for clustering
marker_names = handles.clustering_marker_names(ismember(handles.clustering_marker_names, handles.clustering_used_markers),:); % pooled downsampled data used for clustering
idx = handles.clustering_idx;
remaining_edges = find_matrix_big_element(triu(considered_mst_tree,1),1);
if isempty(remaining_edges)
    return
end
remaining_adj = considered_mst_tree;
score = NaN + zeros(1,size(remaining_edges,1));

[~, ~, ~,components] = extract_connected_component(remaining_adj);
components = components';
parent_var = [];
overall_var = var(data,[],2);
best_score = -Inf;
for i=1:size(components,1)
    parent_var = var(data(:,components(i,idx)==1),[],2)*sum(components(i,idx)==1);
    for j=1:size(remaining_edges,1)
        if sum(components(i,remaining_edges(j,:)))~=2
            continue;
        end
        tmp_adj = remaining_adj;
        tmp_adj(remaining_edges(j,2),remaining_edges(j,1))=0;
        tmp_adj(remaining_edges(j,1),remaining_edges(j,2))=0;
        [~, ~, ~,tmp_components] = extract_connected_component(tmp_adj); 
        tmp_components = tmp_components';
        child_component_1 = tmp_components( find(tmp_components(:,remaining_edges(j,1))), : );
        child_component_2 = tmp_components( find(tmp_components(:,remaining_edges(j,2))), : );
        child_var_1 = var(data(:,child_component_1(idx)==1),[],2)*sum(child_component_1(idx)==1);
        child_var_2 = var(data(:,child_component_2(idx)==1),[],2)*sum(child_component_2(idx)==1);
        score(j) = sum((parent_var - child_var_1 - child_var_2)./overall_var);
        if score(j)>best_score
            best_score = score(j);
            edge_to_remove = remaining_edges(j,:);
            best_child_component_1 = child_component_1;
            best_child_component_2 = child_component_2;
            best_parent_component = child_component_1 + child_component_2;
            best_score_vector = (parent_var - child_var_1 - child_var_2)./overall_var;
        end
    end
end

figure(100)
subplot(5,2,[1 3]);
node_data = handles.clustering_marker_node_average{find(ismember(handles.clustering_marker_node_average(:,1),'POOLED') & ismember(handles.clustering_marker_node_average(:,2),marker_names{best_score_vector==max(best_score_vector)})),3};
tree_annotations=[];
tree_annotations{1} = find(best_child_component_1);
tree_annotations{2} = find(best_child_component_2);
draw_SPADE_tree_annotation(handles.mst_tree, handles.node_positions, handles.node_size, node_data, node_data, 1, 1, best_parent_component, 'jet', tree_annotations);
title(marker_names{best_score_vector==max(best_score_vector)},'interpret','none')
subplot(5,2,[2 4]);
node_data2 = node_data;
node_data2(best_parent_component==0)=NaN;
draw_SPADE_tree_annotation(handles.mst_tree, handles.node_positions, handles.node_size, node_data2, node_data, 1, 1, best_parent_component, 'jet', tree_annotations);
title(marker_names{best_score_vector==max(best_score_vector)},'interpret','none')

data_for_boxplot = [data(:,best_child_component_1(idx)==1), data(:,best_child_component_2(idx)==1), data(:,best_parent_component(idx)==1), data];
group_label_for_boxplot = [ones(1,sum(best_child_component_1(idx))), 2*ones(1, sum(best_child_component_2(idx))), 3*ones(1, sum(best_parent_component(idx))), 4*ones(1,size(data,2))];
for i=1:size(data,1) % length(marker_names)
    h100 = figure(100);
    h0=subplot(5,size(data,1),[size(data,1)*2+i, size(data,1)*3+i]);
    h = boxplot(data_for_boxplot(i,:), group_label_for_boxplot, 'orientation','horizontal','grouporder',[regexp(num2str(4:-1:1),'  ','split')]);
    delete(h([1:4],2:end)); delete(h(7,:));
    axis tight; ylim([0,5]);
    title(marker_names{i},'interpret','none')
    if i==1
        text(min(xlim)-(max(xlim)-min(xlim))/2,1,'ALL')
        text(min(xlim)-(max(xlim)-min(xlim))/2,2,'P')
        text(min(xlim)-(max(xlim)-min(xlim))/2,3,'C2')
        text(min(xlim)-(max(xlim)-min(xlim))/2,4,'C1')
    end
end

handles_accept_button = uicontrol('style', 'pushbutton', 'units', 'normalized', 'position', [0.12 0.1 0.26 0.1], 'Parent',h100, 'String','ACCEPT partition','Callback',{@accept_suggested_annotation, handles, h100,best_child_component_1,best_child_component_2, considered_nodes});
handles_accept_button = uicontrol('style', 'pushbutton', 'units', 'normalized', 'position', [0.39 0.1 0.26 0.1], 'Parent',h100, 'String','ACCEPT & remove parent','Callback',{@accept_suggested_annotation_remove_parent, handles, h100,best_child_component_1,best_child_component_2, considered_nodes});
handles_reject_button = uicontrol('style', 'pushbutton', 'units', 'normalized', 'position', [0.66 0.1 0.26,0.1], 'Parent',h100, 'String','REJECT partition','Callback',{@reject_suggested_annotation, handles, h100, considered_nodes});
1;


function accept_suggested_annotation(hObject, eventdata, handles , h100,best_child_component_1,best_child_component_2, considered_nodes)
% add child component 1
handles.mouse_selected_nodes = find(best_child_component_1);
button_add_selected_nodes_to_annotation_Callback(handles.button_add_selected_nodes_to_annotation, [], handles)
handles = guidata(handles.button_add_selected_nodes_to_annotation);
% add child component 2
handles.mouse_selected_nodes = find(best_child_component_2);
button_add_selected_nodes_to_annotation_Callback(handles.button_add_selected_nodes_to_annotation, [], handles)
handles = guidata(handles.button_add_selected_nodes_to_annotation);
% make parent_annotation selected
parent_annotation = find(best_child_component_1 + best_child_component_2);
handles.mouse_selected_nodes = [];
% if isequal(1:size(handles.mst_tree,1), parent_annotation)
%     handles.mouse_selected_nodes = [];
% else
%     handles.mouse_selected_nodes = parent_annotation;
% end
% set(handles.edit_mouse_selected_nodes,'string',num2str(handles.mouse_selected_nodes));
set(handles.radiobutton_tree_annotation_show_all,'value',1)
radiobutton_tree_annotation_show_all_Callback(handles.radiobutton_tree_annotation_show_all, [], handles);
% % handles = guidata(handles.radiobutton_tree_annotation_show_all);
% % axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
% close the new figure
close(h100);


function accept_suggested_annotation_remove_parent(hObject, eventdata, handles , h100,best_child_component_1,best_child_component_2, considered_nodes)
% add child component 1
handles.mouse_selected_nodes = find(best_child_component_1);
button_add_selected_nodes_to_annotation_Callback(handles.button_add_selected_nodes_to_annotation, [], handles)
handles = guidata(handles.button_add_selected_nodes_to_annotation);
% add child component 2
handles.mouse_selected_nodes = find(best_child_component_2);
button_add_selected_nodes_to_annotation_Callback(handles.button_add_selected_nodes_to_annotation, [], handles)
handles = guidata(handles.button_add_selected_nodes_to_annotation);
% remove the parent
parent_annotation = find(best_child_component_1 + best_child_component_2);
for i=1:length(handles.tree_annotations)
    if isequal(handles.tree_annotations{i},parent_annotation)
        set(handles.listbox_annotations,'value',i);
        button_remove_annotation_Callback(handles.button_remove_annotation,[],handles);
        handles = guidata(handles.button_remove_annotation);
        break;
    end
end
handles.mouse_selected_nodes = [];
% % make parent selected
% if isequal(1:size(handles.mst_tree,1), parent_annotation)
%     handles.mouse_selected_nodes = [];
% else
%     handles.mouse_selected_nodes = parent_annotation;
% end
% set(handles.edit_mouse_selected_nodes,'string',num2str(handles.mouse_selected_nodes));
set(handles.radiobutton_tree_annotation_show_all,'value',1)
radiobutton_tree_annotation_show_all_Callback(handles.radiobutton_tree_annotation_show_all, [], handles);
% % handles = guidata(handles.radiobutton_tree_annotation_show_all);
% % axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
% close the new figure
close(h100);


function reject_suggested_annotation(hObject, eventdata, handles, h100, considered_nodes)
% make considered_nodes selected
if isequal(1:size(handles.mst_tree,1), considered_nodes)
    handles.mouse_selected_nodes = [];
else
    handles.mouse_selected_nodes = considered_nodes;
end
set(handles.edit_mouse_selected_nodes,'string',num2str(handles.mouse_selected_nodes));
set(handles.radiobutton_tree_annotation_show_all,'value',1)
radiobutton_tree_annotation_show_all_Callback(handles.radiobutton_tree_annotation_show_all, [], handles);
% % handles = guidata(handles.radiobutton_tree_annotation_show_all);
% % axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
% close the new figure
close(h100);




% --- Executes on button press in button_view_annotation.
function button_view_annotation_Callback(hObject, eventdata, handles)
% hObject    handle to button_view_annotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.tree_annotations)
    return
end

set(handles.radiobutton_tree_annotation_show_all,'value',1)
radiobutton_tree_annotation_show_all_Callback(handles.radiobutton_tree_annotation_show_all, [], handles);


data = handles.clustering_data(ismember(handles.clustering_marker_names, handles.clustering_used_markers),:); % pooled downsampled data used for clustering
idx = handles.clustering_idx;
marker_names = handles.clustering_marker_names(ismember(handles.clustering_marker_names, handles.clustering_used_markers),:); % pooled downsampled data used for clustering
tree_annotations = handles.tree_annotations;
bubble_contours = handles.tree_bubble_contour; 
mst_tree = handles.mst_tree;
node_positions = handles.node_positions;
node_size = handles.node_size;


h200=figure(200);

all_data = [];
group_id = [];
for i=1:length(tree_annotations)
    ind_tmp = ismember(idx,tree_annotations{i});
    all_data = [all_data, data(:,ind_tmp)];
    group_id = [group_id, i*ones(1,sum(ind_tmp))];
end
all_data = [all_data,data];
group_id = [group_id, length(tree_annotations)+ones(1,size(data,2))];

for i=1:size(data,1)
    handle_boxsubplots(i,1)=subplot(3,length(marker_names),length(marker_names)+i);
    if max(group_id)>=10
        h = boxplot(all_data(i,:), group_id,'orientation','horizontal','grouporder',[regexp(num2str(max(group_id):-1:10),'  ','split'),regexp(num2str(9:-1:1),'  ','split')]);
    else
        h = boxplot(all_data(i,:), group_id,'orientation','horizontal','grouporder',[regexp(num2str(max(group_id):-1:1),'  ','split')]);
    end
    delete(h([1:4],2:end)); delete(h(7,:));
    axis tight; ylim([0,max(group_id)+1]);
    title(marker_names{i},'interpret','none')
    if i==1
        text(min(xlim)-(max(xlim)-min(xlim)),(max(ylim)*0.2+min(ylim)*0.8),'Annotations','rotation',90)
        text(min(xlim)-(max(xlim)-min(xlim))/2,1,'ALL')
    end
end
first_subplot_position = get(handle_boxsubplots(1),'position');
last_subplot_position = get(handle_boxsubplots(2),'position');
handle_boxsubplots_highlight_axes = axes('position',[0,first_subplot_position(2),1,first_subplot_position(4)],'visible','off');
set(handle_boxsubplots_highlight_axes,'ylim',get(handle_boxsubplots(1),'ylim'))
%line([0 1 1 0 0],(max(group_id)+1-1)+0.4*[-1 -1 1 1 -1],'parent',handle_boxsubplots_highlight_axes,'linewidth',1,'color',[0.5,0.5,0.5]);


heatmap_handle = subplot(3,1,3);
annotation_medians = grpstats(all_data',group_id,'median');
for i=1:size(annotation_medians,2) % for each marker
    annotation_medians(:,i) = (annotation_medians(:,i)-min(annotation_medians(:,i)))/(max(annotation_medians(:,i))-min(annotation_medians(:,i)));
end
imagesc(annotation_medians,[0.1,0.9]);
set(heatmap_handle,'YTick',max(group_id),'YTickLabel','ALL','XTick',[])
ylabel('Annotations')
for i=1:size(annotation_medians,2) % for each marker
    text(i,max(group_id)*1.08, marker_names{i},'rotation',-60,'interpret','none')
end
handle_heatmap_highlight_axes = axes('position',get(heatmap_handle,'position'),'visible','off');
set(handle_heatmap_highlight_axes,'xlim',get(heatmap_handle,'xlim'), 'ylim',get(heatmap_handle,'ylim'));


max_group_id = max(group_id);
handles_annotations_listbox  = uicontrol('style', 'listbox', 'units', 'normalized', 'position', [0.5 0.66 0.41 0.3], 'Parent',h200, 'String',get(handles.listbox_annotations,'string'), ...
                                         'Callback',{@highlight_one_annotation, h200, mst_tree, node_positions, node_size, tree_annotations,max_group_id, heatmap_handle,handle_boxsubplots, handle_boxsubplots_highlight_axes, handle_heatmap_highlight_axes});
highlight_one_annotation(handles_annotations_listbox,[],h200, mst_tree, node_positions, node_size, tree_annotations,max_group_id, heatmap_handle,handle_boxsubplots, handle_boxsubplots_highlight_axes, handle_heatmap_highlight_axes);
set(h200,'ResizeFcn',{@highlight_one_annotation_resize, heatmap_handle,handle_boxsubplots, handle_boxsubplots_highlight_axes, handle_heatmap_highlight_axes});



function highlight_one_annotation(hObject, eventdata, h200, mst_tree, node_positions, node_size, tree_annotations,max_group_id, heatmap_handle,handle_boxsubplots, handle_boxsubplots_highlight_axes, handle_heatmap_highlight_axes)
figure(h200);
subplot(3,2,1);
ind = get(hObject,'value');
node_data = zeros(1,size(mst_tree,1));
node_data(tree_annotations{ind})=1;
draw_SPADE_tree(mst_tree, node_positions, node_size, node_data, [-1,1], 1, 0,zeros(1,size(mst_tree,1)), 'JET');

first_subplot_position = get(handle_boxsubplots(1),'position');
last_subplot_position = get(handle_boxsubplots(2),'position');
set(handle_boxsubplots_highlight_axes,'position',[0,first_subplot_position(2),1,first_subplot_position(4)],'visible','off');
set(handle_boxsubplots_highlight_axes,'ylim',get(handle_boxsubplots(1),'ylim'))
axes(handle_boxsubplots_highlight_axes)
delete(get(handle_boxsubplots_highlight_axes,'children'))
line([0 1 1 0 0],(max_group_id+1-ind)+0.4*[-1 -1 1 1 -1],'parent',handle_boxsubplots_highlight_axes,'linewidth',1,'color',[0.5,0.5,0.5]);

heatmap_position = get(heatmap_handle,'position');
set(handle_heatmap_highlight_axes,'position',[0, heatmap_position(2), 1, heatmap_position(4)],'visible','off');
set(handle_heatmap_highlight_axes,'xlim',[0,1], 'ylim',get(heatmap_handle,'ylim'));
axes(handle_heatmap_highlight_axes)
delete(get(handle_heatmap_highlight_axes,'children'))
line([0 10 10 0 0],(max_group_id+1-ind)+0.5*[-1 -1 1 1 -1],'parent',handle_heatmap_highlight_axes,'linewidth',1,'color',[0.5,0.5,0.5]);



function highlight_one_annotation_resize(hObject, eventdata, heatmap_handle,handle_boxsubplots, handle_boxsubplots_highlight_axes, handle_heatmap_highlight_axes)
first_subplot_position = get(handle_boxsubplots(1),'position');
last_subplot_position = get(handle_boxsubplots(2),'position');
set(handle_boxsubplots_highlight_axes,'position',[0,first_subplot_position(2),1,first_subplot_position(4)],'visible','off');
set(handle_boxsubplots_highlight_axes,'ylim',get(handle_boxsubplots(1),'ylim'))

heatmap_position = get(heatmap_handle,'position');
set(handle_heatmap_highlight_axes,'position',[0, heatmap_position(2), 1, heatmap_position(4)],'visible','off');
set(handle_heatmap_highlight_axes,'xlim',[0,1], 'ylim',get(heatmap_handle,'ylim'));




% --- Executes on button press in pushbutton_RadioExpand_layout.
function pushbutton_RadioExpand_layout_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_RadioExpand_layout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

X = handles.clustering_data(ismember(handles.clustering_marker_names,handles.clustering_used_markers),:)';
index = handles.clustering_idx;
clusts = 1:max(handles.clustering_idx);
dist = 'sqeuclidean';
local_density = handles.clustering_local_density;
[centroids, counts, centroids_ind] = gcentroids(X, index, clusts, dist, local_density);
centers = centroids';
node_positions = radio_layout(handles.mst_tree,centers);

% node_positions = radio_layout(handles.mst_tree,cell2mat(handles.clustering_marker_node_average(ismember(handles.clustering_marker_node_average(:,2),handles.clustering_used_markers) & ismember(handles.clustering_marker_node_average(:,1),'POOLED'),3)));

% normalize node positions
node_positions = node_positions - repmat((max(node_positions,[],2)+min(node_positions,[],2))/2,1,size(node_positions,2));
node_positions = node_positions/max(abs(node_positions(:)))*50;
% rotate so that the highest density point is in the west
node_local_density = accumarray(handles.clustering_idx', handles.clustering_local_density')';
% weight_center = sum(node_positions.*repmat(node_local_density,2,1),2)/sum(node_local_density);
weight_center = find_highest_density_position(node_positions, node_local_density);
tmp_score = zeros(1,360) + Inf;
for i=1:360
    tmp_angle = i/180*pi;
    tmp = [cos(tmp_angle), sin(tmp_angle); -sin(tmp_angle), cos(tmp_angle)]*weight_center(:);
    if tmp(1)<0 
        tmp_score(i)=abs(tmp(2));
    end
end
[~,i] = min(tmp_score);
tmp_angle = i/180*pi;
node_positions = [cos(tmp_angle), sin(tmp_angle); -sin(tmp_angle), cos(tmp_angle)]*node_positions;
% flipup to make sure that more density is in the north half
if sum(node_local_density(node_positions(1,:)>0 & node_positions(2,:)>0)) < sum(node_local_density(node_positions(1,:)>0 & node_positions(2,:)<0))
    node_positions(2,:)=-node_positions(2,:);
end


handles.node_positions = node_positions;
handles.tree_bubble_contour = repmat({[]},1,length(handles.tree_annotations));
guidata(hObject,handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);

% viSNE_color_plots(node_positions, ...
%                   cell2mat(handles.clustering_marker_node_average(ismember(handles.clustering_marker_node_average(:,1),'POOLED') & ismember(handles.clustering_marker_node_average(:,2),handles.clustering_used_markers) , 3)), ...
%                   handles.clustering_marker_node_average(ismember(handles.clustering_marker_node_average(:,1),'POOLED') & ismember(handles.clustering_marker_node_average(:,2),handles.clustering_used_markers), 2), ...
%                   node_local_density);
% h = figure;
% % write colored versions
% set(handles.checkbox_is_use_color,'value',1);
% markers_to_show = intersect(handles.clustering_used_markers, get(handles.listbox_select_a_marker,'string'));
% for k=1:length(markers_to_show)
%     set(handles.listbox_select_a_marker,'value',find(ismember(get(handles.listbox_select_a_marker,'string'),markers_to_show(k))));
%     set(handles.listbox_select_a_file,'value',1)
%     subplot(3,3,k); draw_SPADE_tree_no_mouse_click_events(handles); axis off; 
%     title(markers_to_show{k},'interpret','none'); drawnow; 
% end
% set(handles.checkbox_is_use_color,'value',0)
% set(handles.radiobutton_color_definition_expr,'value',1)
% set(handles.radiobutton_color_scheme_JET,'value',1)
% axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);




%copied from the deterministic kmeans code
function [centroids, counts, centroids_ind] = gcentroids(X, index, clusts, dist, local_density)
%GCENTROIDS Centroids and counts stratified by group.
local_density = local_density(:)';
p = size(X,2);
num = length(clusts);
centroids = NaN(num,p);
counts = zeros(num,1);
centroids_ind = zeros(num,1);

for i = 1:num
    members = (index == clusts(i));
    if any(members)
        counts(i) = sum(members);
        switch dist
            case 'sqeuclidean'
                %centroids(i,:) = sum(X(members,:),1) / counts(i);
                tmp = sum(X(members,:).*repmat(local_density(members)',1,size(X,2)),1) / sum(local_density(members));
                [~,ind] = min(sum(abs(X(members,:)-repmat(tmp,counts(i),1)),2));
                centroids_ind(i) = max(find(members==1,ind));
                centroids(i,:) = X(centroids_ind(i),:);
            case 'cityblock'
                % Separate out sorted coords for points in i'th cluster,
                % and use to compute a fast median, component-wise
                members_local_density = local_density(members)';
                [Xsorted,I] = sort(X(members,:),1);
                [~,median_ind] = min(abs(cumsum(members_local_density(I))-sum(members_local_density)/2),[],1);
                centroids(i,:) = diag(Xsorted(median_ind,:));
        end
    end
end



function edit_lower_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lower_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_value = str2num(get(handles.edit_lower_threshold,'string')); 
if isempty(new_value) || new_value>get(handles.edit_upper_threshold,'value')
    set(handles.edit_lower_threshold,'string',num2str(get(handles.edit_lower_threshold,'value'),'%.2f'));
    return
end
set(handles.edit_lower_threshold,'value',new_value,'string',num2str(new_value,'%.2f'));
if handles.color_definition==1
    set(handles.edit_upper_threshold,'value',abs(new_value),'string',num2str(abs(new_value),'%.2f'));
    set(handles.edit_lower_threshold,'value',-abs(new_value),'string',num2str(-abs(new_value),'%.2f'));
end

ind = find(ismember(handles.user_defined_color_thresholds(:,1),handles.color_definition) & ismember(handles.user_defined_color_thresholds(:,2),get(handles.listbox_select_a_marker,'value')));
if isempty(ind)
    handles.user_defined_color_thresholds = [handles.user_defined_color_thresholds; ...
                                             [handles.color_definition, get(handles.listbox_select_a_marker,'value'), get(handles.edit_lower_threshold,'value'), get(handles.edit_upper_threshold,'value')]];
else
    handles.user_defined_color_thresholds(ind,:)=[handles.color_definition, get(handles.listbox_select_a_marker,'value'), get(handles.edit_lower_threshold,'value'), get(handles.edit_upper_threshold,'value')];
end
guidata(hObject, handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);




% --- Executes during object creation, after setting all properties.
function edit_lower_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lower_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_upper_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_upper_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_value = str2num(get(handles.edit_upper_threshold,'string')); 
if isempty(new_value) || new_value<get(handles.edit_lower_threshold,'value')
    set(handles.edit_upper_threshold,'string',num2str(get(handles.edit_upper_threshold,'value'),'%.2f'));
    return
end
set(handles.edit_upper_threshold,'value',new_value,'string',num2str(new_value,'%.2f'));
if handles.color_definition==1
    set(handles.edit_upper_threshold,'value',abs(new_value),'string',num2str(abs(new_value),'%.2f'));
    set(handles.edit_lower_threshold,'value',-abs(new_value),'string',num2str(-abs(new_value),'%.2f'));
end

ind = find(ismember(handles.user_defined_color_thresholds(:,1),handles.color_definition) & ismember(handles.user_defined_color_thresholds(:,2),get(handles.listbox_select_a_marker,'value')));
if isempty(ind)
    handles.user_defined_color_thresholds = [handles.user_defined_color_thresholds; ...
                                             [handles.color_definition, get(handles.listbox_select_a_marker,'value'), get(handles.edit_lower_threshold,'value'), get(handles.edit_upper_threshold,'value')]];
else
    handles.user_defined_color_thresholds(ind,:)=[handles.color_definition, get(handles.listbox_select_a_marker,'value'), get(handles.edit_lower_threshold,'value'), get(handles.edit_upper_threshold,'value')];
end
guidata(hObject, handles);
axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);



% --- Executes during object creation, after setting all properties.
function edit_upper_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_upper_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_default_colorbar_threshold.
function button_default_colorbar_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to button_default_colorbar_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ind = find(ismember(handles.user_defined_color_thresholds(:,1),handles.color_definition) & ismember(handles.user_defined_color_thresholds(:,2),get(handles.listbox_select_a_marker,'value')));
if isempty(ind)
    return
else
    handles.user_defined_color_thresholds(ind,:)=[];
    guidata(hObject, handles);
    axes(handles.Axes_mst); draw_SPADE_tree_when_update(handles);
end
