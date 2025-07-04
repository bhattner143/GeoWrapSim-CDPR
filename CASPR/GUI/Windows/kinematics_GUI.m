% The GUI window for performing kinematics analysis
%
% Author        : Jonathan EDEN
% Created       : 2016
% Description    :

%--------------------------------------------------------------------------
%% Constructor
%--------------------------------------------------------------------------
function varargout = kinematics_GUI(varargin)
    % KINEMATICS_GUI MATLAB code for kinematics_GUI.fig
    %      KINEMATICS_GUI, by itself, creates a new KINEMATICS_GUI or raises the existing
    %      singleton*.
    %
    %      H = KINEMATICS_GUI returns the handle to a new KINEMATICS_GUI or the handle to
    %      the existing singleton*.
    %
    %      KINEMATICS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in KINEMATICS_GUI.M with the given input arguments.
    %
    %      KINEMATICS_GUI('Property','Value',...) creates a new KINEMATICS_GUI or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before kinematics_GUI_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to kinematics_GUI_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help kinematics_GUI

    % Last Modified by GUIDE v2.5 16-Oct-2018 11:32:22

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @kinematics_GUI_OpeningFcn, ...
                       'gui_OutputFcn',  @kinematics_GUI_OutputFcn, ...
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
end

%--------------------------------------------------------------------------
%% GUI Setup Functions
%--------------------------------------------------------------------------
% --- Executes just before kinematics_GUI is made visible.
function kinematics_GUI_OpeningFcn(hObject, ~, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to kinematics_GUI (see VARARGIN)

    % Choose default command line output for kinematics_GUI
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
    % Load the state
    loadState(handles);
    GUIOperations.CreateTabGroup(handles);

    % UIWAIT makes kinematics_GUI wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = kinematics_GUI_OutputFcn(~, ~, handles)
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, handles) %#ok<DEFNU>
    % hObject    handle to figure1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: delete(hObject) closes the figure
    saveState(handles);
    delete(hObject);
end

%--------------------------------------------------------------------------
%% Menu Functions
%--------------------------------------------------------------------------
% --------------------------------------------------------------------
function FileMenu_Callback(~, ~, ~) %#ok<DEFNU>
    % hObject    handle to FileMenu (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function OpenMenuItem_Callback(~, ~, ~) %#ok<DEFNU>
    % hObject    handle to OpenMenuItem (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    file = uigetfile('*.fig');
    if ~isequal(file, 0)
        open(file);
    end
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(~, ~, handles) %#ok<DEFNU>
    % hObject    handle to PrintMenuItem (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    printdlg(handles.figure1)
end

% --------------------------------------------------------------------
function CloseMenuItem_Callback(~, ~, handles) %#ok<DEFNU>
    % hObject    handle to CloseMenuItem (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                         ['Close ' get(handles.figure1,'Name') '...'],...
                         'Yes','No','Yes');
    if strcmp(selection,'No')
        return;
    end

    delete(handles.figure1)
end

%--------------------------------------------------------------------------
%% Popups
%--------------------------------------------------------------------------
% --- Executes on selection change in popupmenu1.
% --- Executes on selection change in plot_type_popup.
function plot_type_popup_Callback(hObject, ~, ~) 
    % hObject    handle to plot_type_popup (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns plot_type_popup contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from plot_type_popup
    settingsXMLObj = GUIOperations.GetSettings('/GUI/XML/kinematicsXML.xml');
    plotsObj = settingsXMLObj.getElementsByTagName('simulator').item(0).getElementsByTagName('plot_functions').item(0).getElementsByTagName('plot_function');
    contents = get(hObject,'Value');
    plotObj = plotsObj.item(contents-1);
    setappdata(hObject,'num_plots',plotObj.getElementsByTagName('figure_quantity').item(0).getFirstChild.getData);
end

% --- Executes during object creation, after setting all properties.
function plot_type_popup_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
    % hObject    handle to plot_type_popup (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    settingsXMLObj = GUIOperations.GetSettings('/GUI/XML/kinematicsXML.xml');
    plot_str = GUIOperations.XmlObj2StringCellArray(settingsXMLObj.getElementsByTagName('simulator').item(0).getElementsByTagName('plot_functions').item(0).getElementsByTagName('plot_function'),'type');
    set(hObject,'Value', 1);    
    set(hObject, 'String', plot_str);
    setappdata(hObject,'num_plots',1);
end

% --- Executes on selection change in trajectory_popup.
function trajectory_popup_Callback(~, ~, ~) %#ok<DEFNU>
    % hObject    handle to trajectory_popup (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns trajectory_popup contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from trajectory_popup
end

function trajectory_popup_Update(~, handles)
    contents = cellstr(get(handles.model_text,'String'));
    model_type = contents{1};
    if(getappdata(handles.figure1,'toggle'))
        model_config = DevModelConfig(model_type);
    else
        model_config = ModelConfig(model_type);
    end
    setappdata(handles.trajectory_popup,'model_config',model_config);
    % Determine the trajectories
    trajectories_str = model_config.getJointTrajectoriesList();    
    if(~isempty(trajectories_str))
        set(handles.trajectory_popup, 'Value', 1);   set(handles.trajectory_popup, 'String', trajectories_str);
    else
        set(handles.trajectory_popup, 'Value', 1);   set(handles.trajectory_popup, 'String', 'No trajectories are defined for this robot');
    end
end

% --- Executes during object creation, after setting all properties.
function trajectory_popup_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
    % hObject    handle to trajectory_popup (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% --- Executes on selection change in solver_class_popup.
function solver_class_popup_Callback(~, ~, handles) %#ok<DEFNU>
    % hObject    handle to solver_class_popup (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns solver_class_popup contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from solver_class_popup
    
end

function solver_class_popup_Update(hObject,handles)
    contents = cellstr(get(handles.kinematics_popup,'String'));
    kinematics_id = contents{get(handles.kinematics_popup,'Value')};
    if(strcmp(kinematics_id,'Inverse Kinematics'))
        set(hObject,'value',1);
        set(hObject, 'String', {' '});
        set(hObject,'Visible','off');
        set(handles.solver_class_text,'Visible','off');
    else
        set(hObject,'Visible','on');
        set(handles.solver_class_text,'Visible','on');
        settings = getappdata(handles.figure1,'settings');
        solver_str = GUIOperations.XmlObj2StringCellArray(settings.getElementsByTagName('simulator').item(0).getElementsByTagName('solver_class'),'id');
        set(hObject,'Value',1);
        set(hObject, 'String', solver_str);
    end
end

% --- Executes during object creation, after setting all properties.
function solver_class_popup_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
    % hObject    handle to solver_class_popup (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end 
end

% --- Executes on selection change in kinematics_popup.
function kinematics_popup_Callback(hObject, ~, handles) %#ok<DEFNU>
    % hObject    handle to kinematics_popup (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns kinematics_popup contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from kinematics_popup
    contents = cellstr(get(hObject,'String'));
    toggle_visibility(contents{get(hObject,'Value')},handles);
end

% --- Executes during object creation, after setting all properties.
function kinematics_popup_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
    % hObject    handle to kinematics_popup (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    kinematics_str = {'Forward Kinematics','Inverse Kinematics'};
    set(hObject, 'String', kinematics_str);
end
    
%--------------------------------------------------------------------------
%% Push Buttons

function plot_button_Callback(~, ~, handles) %#ok<DEFNU>
    % hObject    handle to plot_button (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    sim = getappdata(handles.figure1,'sim');
    if(isempty(sim))
        warning('No simulator has been generated. Please press run first'); %#ok<WNTAG>
    else
        contents = cellstr(get(handles.plot_type_popup,'String'));
        plot_type = contents{get(handles.plot_type_popup,'Value')};
        GUIOperations.GUIPlot(plot_type,sim,handles,str2double(getappdata(handles.plot_type_popup,'num_plots')),get(handles.undock_box,'Value'))
    end
end


% --- Executes on button press in run_button.
function run_button_Callback(~, ~, handles) %#ok<DEFNU>
    % hObject    handle to run_button (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    clc; 
    contents = cellstr(get(handles.trajectory_popup,'String'));
    trajectory_id = contents{get(handles.trajectory_popup,'Value')};
    % Then read the form of dynamics
    contents = cellstr(get(handles.kinematics_popup,'String'));
    kinematics_id = contents{get(handles.kinematics_popup,'Value')};
    modObj = getappdata(handles.cable_text,'modObj');
    if(strcmp(kinematics_id,'Inverse Kinematics'))
        run_inverse_kinematics(handles,modObj,trajectory_id);
    else
        run_forward_kinematics(handles,modObj,trajectory_id);
    end
end

% --- Executes on button press in plot_movie_button.
function plot_movie_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to plot_movie_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Things that this should do
    sim = getappdata(handles.figure1,'sim');
    if(isempty(sim))
        warning('No simulator has been generated. Please press run first'); %#ok<WNTAG>
    else
        path_string = CASPR_configuration.LoadHomePath();
        % Check if the log folder exists
        if((exist([path_string,'/data'],'dir')~=7)||(exist([path_string,'/data/videos'],'dir')~=7))
            if((exist([path_string,'/data'],'dir')~=7))
                mkdir([path_string,'/data']);
            end
            mkdir([path_string,'/data/videos']);
        end      
        model_config = getappdata(handles.trajectory_popup,'model_config');
        file_name = [path_string,'/data/videos/kinematics_gui_output.avi'];
        [file,path] = uiputfile(file_name,'Save file name');
        sim.plotMovie(model_config.displayRange, model_config.viewAngle, [path,file], sim.timeVector(length(sim.timeVector)), false, 700, 700);
    end
end

% --- Executes on button press in script_button.
function script_button_Callback(~, ~, handles) %#ok<DEFNU>
    % hObject    handle to script_button (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Determine if forward or inverse kinematics
    contents = cellstr(get(handles.kinematics_popup,'String'));
    kinematics_id = contents{get(handles.kinematics_popup,'Value')};
    % Strings for model editing
    model_str = cellstr(get(handles.model_text,'String'));
    cable_str = cellstr(get(handles.cable_text,'String'));
    contents = cellstr(get(handles.trajectory_popup,'String'));
    trajectory_str = contents{get(handles.trajectory_popup,'Value')};
    if(strcmp(kinematics_id,'Inverse Kinematics'))
        base_folder = CASPR_configuration.LoadHomePath();
        r_string = [base_folder,'/GUI/template_scripts/kinematics/script_IK_template.m'];
    else
        contents = cellstr(get(handles.solver_class_popup,'String'));
        solver_class_id = contents{get(handles.solver_class_popup,'Value')};
        base_folder = CASPR_configuration.LoadHomePath();
        if(strcmp(solver_class_id,'FKLeastSquares'))
            r_string = [base_folder,'/GUI/template_scripts/kinematics/script_FK_least_squares_template.m'];
        elseif(strcmp(solver_class_id,'FKDifferential'))
            r_string = [base_folder,'/GUI/template_scripts/kinematics/script_FK_differential_template.m'];
        elseif(strcmp(solver_class_id,'FKHybridLeastSquaresDifferential'))
            r_string = [base_folder,'/GUI/template_scripts/kinematics/script_FK_hybrid_leastsquares_differential_template.m'];
        end
        
    end
    w_string = [base_folder,'/scripts/local/GUI_script_autogenerated.m'];
    r_fid = fopen(r_string,'r');
    w_fid = fopen(w_string,'w');
    while(~feof(r_fid))
        s = fgetl(r_fid);
        % Determine if comment
        new_s = regexprep(s,'%','%%');
        % Replace all references to the model
        new_s = regexprep(new_s,'Example planar XY',model_str);
        new_s = regexprep(new_s,'basic',cable_str);
        new_s = regexprep(new_s,'example_linear',trajectory_str);
        if (CASPR_configuration.LoadDevModelConfig())
            new_s = regexprep(new_s, 'ModelConfig', 'DevModelConfig');
        end
        fprintf(w_fid,[new_s,'\n']);
    end
    fclose(r_fid);
    fclose(w_fid);
    edit(w_string)
end

function Rviz_pushbutton_Callback(~, ~, handles) %#ok<DEFNU>
    % hObject    handle to Rviz_pushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    sim = getappdata(handles.figure1,'sim');
    modObj = getappdata(handles.cable_text,'modObj');
    if(isempty(sim))
        warning('No simulator has been generated. Please press run first'); %#ok<WNTAG>
    else
        try 
            load('CARDSFlowConfig.mat');
            MotionSimulatorBase.PlotCARDSFlow(modObj, sim.trajectory);       
        catch 
            MotionSimulatorBase.PlotRviz(modObj, sim.trajectory);     
        end       
    end
end


% --- Executes on button press in update_button.
function update_button_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to update_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    trajectory_popup_Update(hObject, handles);
end

%--------------------------------------------------------------------------
% Checkboxes
%--------------------------------------------------------------------------
% --- Executes on button press in undock_box.
function undock_box_Callback(~, ~, ~) %#ok<DEFNU>
    % hObject    handle to undock_box (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of undock_box
    % --- Executes on button press in plot_button.
end

%--------------------------------------------------------------------------
%% Toolbar buttons
%--------------------------------------------------------------------------
% --------------------------------------------------------------------
function save_figure_tool_ClickedCallback(~, ~, handles) %#ok<DEFNU>
    % hObject    handle to save_figure_tool (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    tabgp = getappdata(handles.figure1,'tabgp');
    s_tab = get(tabgp,'SelectedTab');
    ax = get(s_tab,'Children');
    f1 = figure; % Open a new figure with handle f1
    copyobj(ax,f1); % Copy axes object h into figure f1
    set(gca,'ActivePositionProperty','outerposition')
    set(gca,'Units','normalized')
    set(gca,'OuterPosition',[0 0 1 1])
    set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
    [file,path] = uiputfile({'*.fig';'*.bmp';'*.eps';'*.emf';'*.jpg';'*.pcx';...
        '*.pbm';'*.pdf';'*.pgm';'*.png';'*.ppm';'*.svg';'*.tif'},'Save file name');
    if(path ~= 0)
        saveas(gcf,[path,file]);
    end
    close(f1);
end

% --------------------------------------------------------------------
function undock_figure_tool_ClickedCallback(~, ~, handles) %#ok<DEFNU>
    % hObject    handle to undock_figure_tool (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    tabgp = getappdata(handles.figure1,'tabgp');
    s_tab = get(tabgp,'SelectedTab');
    ax = get(s_tab,'Children');
    f1 = figure; % Open a new figure with handle f1
    copyobj(ax,f1); % Copy axes object h into figure f1
    set(gca,'ActivePositionProperty','outerposition')
    set(gca,'Units','normalized')
    set(gca,'OuterPosition',[0 0 1 1])
    set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
end

% --------------------------------------------------------------------
function delete_figure_tool_ClickedCallback(~, ~, handles) %#ok<DEFNU>
    % hObject    handle to delete_figure_tool (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    tabgp = getappdata(handles.figure1,'tabgp');
    s_tab = get(tabgp,'SelectedTab');
    if(strcmp('0',get(s_tab,'Title')))
        % Do nothing
    else
        delete(s_tab);
    end
end

%--------------------------------------------------------------------------
% Additional Functions
%--------------------------------------------------------------------------
function run_inverse_kinematics(handles,modObj,trajectory_id)
    % plotting information
    contents = cellstr(get(handles.plot_type_popup,'String'));
    plot_type = contents{get(handles.plot_type_popup,'Value')};
    % Setup the inverse kinematics simulator with the SystemKinematics object
    CASPR_log.Info('Start Setup Simulation');
    set(handles.status_text,'String','Setting up simulation');
    drawnow;
    start_tic = tic;
    sim = InverseKinematicsSimulator(modObj);
    model_config = getappdata(handles.trajectory_popup,'model_config');
    trajectory = model_config.getJointTrajectory(trajectory_id);
    time_elapsed = toc(start_tic);
    CASPR_log.Info(sprintf('End Setup Simulation : %f seconds', time_elapsed));

    % Run the kinematics on the desired trajectory
    CASPR_log.Info('Start Running Simulation');
    set(handles.status_text,'String','Simulation running');
    drawnow;
    start_tic = tic;
    sim.run(trajectory);
    time_elapsed = toc(start_tic);
    CASPR_log.Info(sprintf('End Running Simulation : %f seconds', time_elapsed));

    % After running the simulator the data can be plotted
    % Refer to the simulator classes to see what can be plotted.
    CASPR_log.Info('Start Plotting Simulation');
    start_tic = tic;
    set(handles.status_text,'String','Simulation plotting');
    drawnow;
    GUIOperations.GUIPlot(plot_type,sim,handles,str2double(getappdata(handles.plot_type_popup,'num_plots')),get(handles.undock_box,'Value'));
    time_elapsed = toc(start_tic);
    CASPR_log.Info(sprintf('End Plotting Simulation : %f seconds', time_elapsed));
    set(handles.status_text,'String','No simulation running');
    setappdata(handles.figure1,'sim',sim);
    assignin('base','inverse_kinematic_simulator',sim);
end

function run_forward_kinematics(handles,modObj,trajectory_id)
    % First load all the necessary informatino
    contents = cellstr(get(handles.solver_class_popup,'String'));
    solver_class = contents{get(handles.solver_class_popup,'Value')};
    
    % An inverse kinematics and forward kinematics simulator will be run to
    % show that the results are consistent.
    CASPR_log.Info('Start Setup Simulation');
    set(handles.status_text,'String','Setting up simulation');
    drawnow;
    start_tic = tic;
    % Initialise the least squares solver for the forward kinematics
    if(strcmp(solver_class,'FKLeastSquares'))
        fksolver = FKLeastSquares(modObj, FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_PSEUDOINV, FK_LS_QdotOptionType.FIRST_ORDER_DERIV);
    elseif(strcmp(solver_class,'FKDifferential'))
        fksolver = FKDifferential(modObj);
    elseif(strcmp(solver_class,'FKHybridLeastSquaresDifferential'))
        fksolver = FKHybridLeastSquaresDifferential(modObj, FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_PSEUDOINV, FK_LS_QdotOptionType.FIRST_ORDER_DERIV, 10);
    else
        CASPR_log.Print(CASPR_log.Error,'Unknown forward kinematic solver used');
    end
      
    
    % Initialise the three inverse/forward kinematics solvers
    iksim = InverseKinematicsSimulator(modObj);
    fksim = ForwardKinematicsSimulator(modObj, fksolver);
    model_config = getappdata(handles.trajectory_popup,'model_config');
    trajectory = model_config.getJointTrajectory(trajectory_id);
    time_elapsed = toc(start_tic);
    CASPR_log.Info(sprintf('End Setup Simulation : %f seconds', time_elapsed));

    % Run inverse kinematics
    CASPR_log.Info('Start Running Inverse Kinematics Simulation');
    set(handles.status_text,'String','Running inverse kinematics');
    drawnow;
    start_tic = tic;
    iksim.run(trajectory);
    time_elapsed = toc(start_tic);
    CASPR_log.Info(sprintf('End Running Inverse Kinematics Simulation : %f seconds', time_elapsed));

    % Run forward kinematics
    CASPR_log.Info('Start Running Forward Kinematics Simulation');
    set(handles.status_text,'String','Running forward kinematics');
    drawnow;
    start_tic = tic;
    fksim.run(iksim.cableLengths, iksim.cableLengthsDot, iksim.timeVector, iksim.trajectory.q{1}, iksim.trajectory.q_dot{1});
    time_elapsed = toc(start_tic);
    CASPR_log.Info(sprintf('End Running Forward Kinematics Simulation : %f seconds', time_elapsed));

    % It is expected that iksim and fksim should have the same joint space (the
    % result of fksim)
    CASPR_log.Info('Start Plotting Simulation');
    set(handles.status_text,'String','Simulation plotting');
    drawnow;
    start_tic = tic;
    GUIOperations.GUIPlot('guiPlotJointPose',iksim,handles,1,get(handles.undock_box,'Value'));
    GUIOperations.GUIPlot('guiPlotJointPose',fksim,handles,1,get(handles.undock_box,'Value'));
    GUIOperations.GUIPlot('guiPlotCableLengthError',fksim,handles,1,get(handles.undock_box,'Value'));
    time_elapsed = toc(start_tic);
    CASPR_log.Info(sprintf('End Plotting Simulation : %f seconds', time_elapsed));
    set(handles.status_text,'String','No simulation running');
    assignin('base','forward_kinematic_simulator',fksim);
    % Save the sim
    setappdata(handles.figure1,'sim',fksim);
end

function loadState(handles)
    % load all of the settings and initialise the values to match
    path_string = CASPR_configuration.LoadHomePath();
    settingsXMLObj =  XmlOperations.XmlReadRemoveIndents([path_string,'/GUI/XML/kinematicsXML.xml']);
    setappdata(handles.figure1,'settings',settingsXMLObj);
    file_name = [path_string,'/GUI/config/caspr_gui_state.mat'];
    set(handles.status_text,'String','No simulation running');
    if(exist(file_name,'file'))
        load(file_name);
        set(handles.model_text,'String',state.model_text);
        set(handles.cable_text,'String',state.cable_text);
        setappdata(handles.figure1,'toggle',CASPR_configuration.LoadDevModelConfig());
        setappdata(handles.cable_text,'modObj',state.modObj);
        trajectory_popup_Update([], handles);
        file_name = [path_string,'/GUI/config/kinematics_gui_state.mat'];
        if(exist(file_name,'file'))
            load(file_name);
            mp_text = get(handles.model_text,'String');
            cs_text = get(handles.cable_text,'String');
            if(strcmp(mp_text,state.model_text)&&strcmp(cs_text,state.cable_text))
                set(handles.trajectory_popup,'value',state.trajectory_popup_vale);
                set(handles.kinematics_popup,'value',state.kinematics_popup_vale);
                solver_class_popup_Update(handles.solver_class_popup,handles);
                set(handles.solver_class_popup,'value',state.solver_class_popup);
                set(handles.plot_type_popup,'value',state.plot_type_popup);
                plot_type_popup_Callback(handles.plot_type_popup,[],handles);
            else
                initialise_popups(handles);
            end
        else
            initialise_popups(handles);
        end
    end
end

function saveState(handles,file_path)
    % Save all of the settings
    state.simulator                         =   'kinematics';
    state.model_text                        =   get(handles.model_text,'String');
    state.cable_text                        =   get(handles.cable_text,'String');
    state.trajectory_popup_vale             =   get(handles.trajectory_popup,'value');
    state.kinematics_popup_vale             =   get(handles.kinematics_popup,'value');
    state.solver_class_popup                =   get(handles.solver_class_popup,'value');
    state.plot_type_popup                   =   get(handles.plot_type_popup,'value');
    if(nargin>1)
        save(file_path,'state');
    else
        path_string                         =   CASPR_configuration.LoadHomePath();
        save([path_string,'/GUI/config/kinematics_gui_state.mat'],'state');
    end
end

function toggle_visibility(dynamics_method,handles)
    if(strcmp(dynamics_method,'Inverse Kinematics'))
        % Forward kinematics so hide all of the options
        set(handles.solver_class_text,'Visible','off');
        set(handles.solver_class_popup,'Visible','off');
        set(handles.plot_type_text,'Visible','on');
        set(handles.plot_type_popup,'Visible','on');
    else
        % Inverse kinematics so let all of the options be viewed
        set(handles.solver_class_text,'Visible','on');
        set(handles.solver_class_popup,'Visible','on');
        set(handles.plot_type_text,'Visible','off');
        set(handles.plot_type_popup,'Visible','off');
        solver_class_popup_Update(handles.solver_class_popup,handles);
    end
end

function initialise_popups(handles)
    % Updates
    solver_class_popup_Update(handles.solver_class_popup,handles);
    % Needed callbacks
    plot_type_popup_Callback(handles.plot_type_popup,[],handles);
end


