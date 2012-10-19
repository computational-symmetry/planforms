function varargout = symmetry_sim(varargin)
% SYMMETRY_SIM M-file for symmetry_sim.fig
%      SYMMETRY_SIM, by itself, creates a new SYMMETRY_SIM or raises the existing
%      singleton*.
%
%      H = SYMMETRY_SIM returns the handle to a new SYMMETRY_SIM or the handle to
%      the existing singleton*.
%
%      SYMMETRY_SIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYMMETRY_SIM.M with the given input arguments.
%
%      SYMMETRY_SIM('Property','Value',...) creates a new SYMMETRY_SIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before symmetry_sim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to symmetry_sim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help symmetry_sim

% Last Modified by GUIDE v2.5 09-Dec-2010 10:58:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @symmetry_sim_OpeningFcn, ...
                   'gui_OutputFcn',  @symmetry_sim_OutputFcn, ...
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


% --- Executes just before symmetry_sim is made visible.
function symmetry_sim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to symmetry_sim (see VARARGIN)

% Choose default command line output for symmetry_sim
handles.output = hObject;

%---- Get defaults from gui fig
% alpha = get(handles.alpha_slider, 'Value');
% beta  = get(handles.beta_slider, 'Value');
alphabeta = eval(get( handles.alphabeta_edit, 'String' ) );
angle_rad = atan(alphabeta);
set( handles.theta_edit, 'String', num2str( angle_rad ) );

cyc_per_img = str2double( get( handles.cyc_per_img_edit, 'String') );
base_angle_offset = str2double( get( handles.base_angle_offset_edit, 'String') );
gaussian_space_constant = str2double( get( handles.gaussian_mask_edit, 'String') );

contents = cellstr(get(handles.Display_popupmenu,'String')); 
selection = contents{get(handles.Display_popupmenu,'Value')};
handles.plot_square_type = selection;

contents = cellstr(get(handles.Class_4_6_popupmenu,'String')); 
selection = contents{get(handles.Class_4_6_popupmenu,'Value')};
handles.class_type = selection;

ampl = 1;
phase_rad = 0;
img_pix = 512;
gray_scale = 255;

% Colormap
gm = linspace(0,1, 256)';
gray_map = [gm gm gm];


% handles.alpha = alpha;
% handles.beta  = beta;
handles.angle_rad = angle_rad;
handles.ampl  = ampl;
handles.phase_rad = phase_rad;
handles.cyc_per_img = cyc_per_img;
handles.img_pix = img_pix;
handles.gray_scale = 255;
handles.gray_map = gray_map;
handles.base_angle_offset = base_angle_offset;
handles.gaussian_space_constant = gaussian_space_constant;
handles.sqsupsq_phase_rad = pi;
handles.pair_angle_offset = 0;
handles.plotaxis_h = [];
handles.plotfig_h = [];
handles.parentfig_h = gcbf;
handles.theta_mode = 0; % Default is to change alpha, beta

handles = compute_dependencies( handles );

handles = make_planforms( handles );
handles = plot_planforms( handles );

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes symmetry_sim wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = symmetry_sim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function img_pix_edit_Callback(hObject, eventdata, handles)
% hObject    handle to img_pix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_val = str2double( get(hObject,'String') );
% if new_val > get(hObject, 'Max')
%     new_val = get(hObject, 'Max');
% end
% if new_val < get(hObject, 'Min');
%     new_val = get(hObject', 'Min');
% end
handles.img_pix = new_val;

handles = compute_dependencies( handles );
handles = make_planforms( handles );
handles = plot_planforms( handles );

guidata( hObject, handles );
%--------------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function img_pix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img_pix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cyc_per_img_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cyc_per_img_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cyc_per_img_edit as text
%        str2double(get(hObject,'String')) returns contents of cyc_per_img_edit as a double
new_val = str2double( get(hObject,'String') );
% if new_val > get(hObject, 'Max')
%     new_val = get(hObject, 'Max');
% end
% if new_val < get(hObject, 'Min');
%     new_val = get(hObject', 'Min');
% end
handles.cyc_per_img = new_val;

handles = compute_dependencies( handles );
handles = make_planforms( handles );
handles = plot_planforms( handles );

guidata( hObject, handles );
%--------------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function cyc_per_img_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cyc_per_img_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = make_planforms( handles );

%---    Extract alpha, beta from handles struct
% alpha = handles.alpha;
% beta = handles.beta;

%---- Reparameterized grating
X = handles.X;
Y = handles.Y;
ampl = handles.ampl;

angle_rad = handles.angle_rad; % + handles.base_angle_offset;
cyc_pix = handles.cyc_pix;
gray_scale = handles.gray_scale;
phase_rad = handles.phase_rad;
sqsupsq_phase_rad = handles.sqsupsq_phase_rad;
class_type  = handles.class_type;    

base_angle_offset = handles.base_angle_offset;

r6 = 2*rand(1,6)-ones(1,6);
a_noise = r6 .* ones(1,6)*0;
handles.a_noise = a_noise;

% C1 = gs_grating( X, Y, ampl, phase_rad, angle_rad + a_noise(1) + base_angle_offset, cyc_pix, gray_scale );
% C2 = gs_grating( X, Y, ampl, phase_rad, angle_rad + pi/2 + a_noise(2)+ base_angle_offset, cyc_pix, gray_scale );
% 
% C3 = gs_grating( X, Y, ampl, phase_rad, pi/2 - angle_rad + a_noise(3)+ base_angle_offset, cyc_pix, gray_scale );
% C4 = gs_grating( X, Y, ampl, phase_rad, -angle_rad + a_noise(4)+ base_angle_offset, cyc_pix, gray_scale );
% 
% C5 = gs_grating( X, Y, ampl, phase_rad + sqsupsq_phase_rad, pi/2 - angle_rad + a_noise(3)+ base_angle_offset, cyc_pix, gray_scale ); % c3 with pi phase shift
% C6 = gs_grating( X, Y, ampl, phase_rad + sqsupsq_phase_rad, -angle_rad + a_noise(4)+ base_angle_offset, cyc_pix, gray_scale );       % c4 with pi phase shift

switch class_type
    case 'Class 4'       
        pair_angle = pi/2 + handles.pair_angle_offset;
        plot_n_components = 6;
        plot_components = {'C1', 'C2', 'C3', 'C4' };
        base_angle = ones(1, plot_n_components )*base_angle_offset;
        comp_axis  = [ 0 pair_angle pair_angle 0 ];
        comp_offset = [ angle_rad angle_rad -angle_rad -angle_rad ];
        
        C1 = gs_grating( X, Y, ampl, phase_rad, base_angle(1) + comp_axis(1) + comp_offset(1) + a_noise(1), cyc_pix, 1 );
        C2 = gs_grating( X, Y, ampl, phase_rad, base_angle(2) + comp_axis(2) + comp_offset(2) + a_noise(2), cyc_pix, 1 );
        C3 = gs_grating( X, Y, ampl, phase_rad, base_angle(3) + comp_axis(3) + comp_offset(3) + a_noise(3), cyc_pix, 1 );
        C4 = gs_grating( X, Y, ampl, phase_rad, base_angle(4) + comp_axis(4) + comp_offset(4) + a_noise(4), cyc_pix, 1 );
        C5 = gs_grating( X, Y, ampl, phase_rad + sqsupsq_phase_rad, base_angle(3) + comp_axis(3) + comp_offset(3) + a_noise(3), cyc_pix, 1 );
        C6 = gs_grating( X, Y, ampl, phase_rad + sqsupsq_phase_rad, base_angle(4) + comp_axis(4) + comp_offset(4) + a_noise(4), cyc_pix, 1 );
       
%         C3 = gs_grating( X, Y, ampl, phase_rad, base_angle(3) + comp_axis(3) + comp_offset(3) + a_noise(3), cyc_pix, 1 );
%         C4 = gs_grating( X, Y, ampl, phase_rad, base_angle(4) + comp_axis(4) + comp_offset(4) + a_noise(4), cyc_pix, 1 );
% 
%         C1 = gs_grating( X, Y, ampl, phase_rad + sqsupsq_phase_rad, angle_rad + a_noise(1) + base_angle_offset, cyc_pix, 1 );
%         C2 = gs_grating( X, Y, ampl, phase_rad + sqsupsq_phase_rad, angle_rad + pair_angle + a_noise(2)+ base_angle_offset, cyc_pix, 1 );
%         
%         C3 = gs_grating( X, Y, ampl, phase_rad, pair_angle - angle_rad + a_noise(3)+ base_angle_offset, cyc_pix, 1 );
%         C4 = gs_grating( X, Y, ampl, phase_rad, -angle_rad + a_noise(4)+ base_angle_offset, cyc_pix, 1 );
        
%         C5 = gs_grating( X, Y, ampl, phase_rad + sqsupsq_phase_rad, pi/2 - angle_rad + a_noise(3)+ base_angle_offset, cyc_pix, 1 ); % c3 with pi phase shift
%         C6 = gs_grating( X, Y, ampl, phase_rad + sqsupsq_phase_rad, -angle_rad + a_noise(4)+ base_angle_offset, cyc_pix, 1 );       % c4 with pi phase shift
    case 'Class 6'
        pair_angle = 2*pi/3 + handles.pair_angle_offset;
        plot_n_components = 6;
        plot_components = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6' };
        
        base_angle = ones(1, 6 )*base_angle_offset;
        comp_axis  = [ 0 pair_angle pair_angle 0 2*pair_angle 2*pair_angle ];
        comp_offset = [ angle_rad angle_rad -angle_rad -angle_rad angle_rad -angle_rad];
                
        C1 = gs_grating( X, Y, ampl, phase_rad, base_angle(1) + comp_axis(1) + comp_offset(1) + a_noise(1), cyc_pix, 1 );
        C2 = gs_grating( X, Y, ampl, phase_rad, base_angle(2) + comp_axis(2) + comp_offset(2) + a_noise(2), cyc_pix, 1 );
        
        C3 = gs_grating( X, Y, ampl, phase_rad, base_angle(3) + comp_axis(3) + comp_offset(3) + a_noise(3), cyc_pix, 1 );
        C4 = gs_grating( X, Y, ampl, phase_rad, base_angle(4) + comp_axis(4) + comp_offset(4) + a_noise(4), cyc_pix, 1 );
        
        C5 = gs_grating( X, Y, ampl, phase_rad, base_angle(5) + comp_axis(5) + comp_offset(5) + a_noise(5), cyc_pix, 1 ); 
        C6 = gs_grating( X, Y, ampl, phase_rad, base_angle(6) + comp_axis(6) + comp_offset(6) + a_noise(6), cyc_pix, 1 );      
end


%----   Compute plaids/planforms
scale=max(max(handles.gaussian_mask));
mask = handles.gaussian_mask;
imgstruct.mask = gray_scale*mask;

P1234=((C1+C2+C3+C4).*mask)/4;
P1234=(P1234*gray_scale/2) + ones( size(P1234) )*gray_scale/2;
sq = P1234;
imgstruct.P1234 = P1234;

P1256=(C1+C2+C5+C6).*mask/4;
P1256=(P1256*gray_scale/2) + ones( size(P1256) )*gray_scale/2;
supsq = P1256;
imgstruct.P1256 = P1256;

P12 = (C1 + C2).*mask/2;
scale=max(max(P12));
P12=(P12*gray_scale/2) + ones( size(P12) )*gray_scale/2;
imgstruct.P12 = P12;

P34 = (C3 + C4).*mask/2;
scale=max(max(P34));
P34=(P34*gray_scale/2) + ones( size(P34) )*gray_scale/2;
imgstruct.P34 = P34;

P56 = (C5 + C6).*mask/2;
scale=max(max(P56));
P56=(P56*gray_scale/2) + ones( size(P56) )*gray_scale/2;
imgstruct.P56 = P56;

P123456 = ((C1+C2+C3+C4+C5+C6).*mask)/6;
P123456=(P123456*gray_scale/2) + ones( size(P123456) )*gray_scale/2;
imgstruct.P123456 = P123456;

C1=C1.*mask;
C1=(C1*gray_scale/2) + ones( size(C1) )*gray_scale/2;
imgstruct.C1 = C1;

C2=C2.*mask;
C2=(C2*gray_scale/2) + ones( size(C2) )*gray_scale/2;
imgstruct.C2 = C2;

C3=C3.*mask;
C3=(C3*gray_scale/2) + ones( size(C3) )*gray_scale/2;
imgstruct.C3 = C3;

C4=C4.*mask;
C4=(C4*gray_scale/2) + ones( size(C4) )*gray_scale/2;
imgstruct.C4 = C4;

C5=C3.*mask;
C5=(C5*gray_scale/2) + ones( size(C5) )*gray_scale/2;
imgstruct.C5 = C5;

C6=C4.*mask;
C6=(C6*gray_scale/2) + ones( size(C6) )*gray_scale/2;
imgstruct.C6 = C6;

handles.imgstruct = imgstruct;

return;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function handles = plot_planforms( handles )

switch handles.plot_square_type
    case 'Square'
        switch( handles.class_type )
            case 'Class 4'
                plotimg = handles.imgstruct.P1234; % Make more general later;
            case 'Class 6'
                 plotimg = handles.imgstruct.P123456; % Make more general later;
        end
    case 'Super Square'
        plotimg = handles.imgstruct.P1256;
    case 'C1+C2'
        plotimg = handles.imgstruct.P12;
    case 'C3+C4'
        plotimg = handles.imgstruct.P34;
    case 'C5+C6'
        plotimg = handles.imgstruct.P56;
    case 'C1'
        plotimg = handles.imgstruct.C1;
    case 'C2'
        plotimg = handles.imgstruct.C2;
    case 'C3'
        plotimg = handles.imgstruct.C3;
    case 'C4'
        plotimg = handles.imgstruct.C4;
    case 'C5'
        plotimg = handles.imgstruct.C5;
    case 'C6'
        plotimg = handles.imgstruct.C6;
    case 'Mask'
        plotimg = handles.imgstruct.mask;
end

% If separate printable fig & axis exist, show them
if ishandle( handles.plotfig_h )
    figure( handles.plotfig_h );
    a2 = image( plotimg, 'Parent', handles.plotaxis_h );
    axis(handles.plotaxis_h, 'square', 'off' );
    colormap( handles.gray_map );
end

% Show inline plot
% handles.parentfig_h = gcbf;
% figure( handles.parentfig_h );

a1 = image( plotimg, 'Parent', handles.axes1 );
colormap( handles.gray_map );
axis(handles.axes1, 'square', 'off' );
handles.fig_handle = a1;

return;
%--------------------------------------------------------------------------


function handles = compute_dependencies( handles )

[handles.X, handles.Y] = meshgrid(0:handles.img_pix-1);
handles.cyc_pix = handles.cyc_per_img/handles.img_pix;
handles.gaussian_mask = make_gaussian_mask( handles );

return;

%--------------------------------------------------------------------------
function c = gs_grating( X, Y, ampl, phase_rad, angle_rad, cyc_pix, gray_scale )
% Generates gray scale grating

f  = cyc_pix*2*pi;
aa = cos( angle_rad )*f;
bb = sin( angle_rad )*f;

c = ampl*cos( aa*X + bb*Y +  phase_rad );

% c = gray_scale/2 * cc + ones( size( cc ) ) * gray_scale/2;

return
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function cc = grating( X, Y, ampl, phase_rad, angle_rad, cyc_pix )
% Generates gray scale grating

f  = cyc_pix*2*pi;
aa = cos( angle_rad )*f;
bb = sin( angle_rad )*f;

cc = ampl*cos( aa*X + bb*Y +  phase_rad );

% c = gray_scale/2 * cc + ones( size( cc ) ) * gray_scale/2;

return
%--------------------------------------------------------------------------


% --- Executes on selection change in Display_popupmenu.
function Display_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Display_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Display_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Display_popupmenu
contents = cellstr(get(hObject,'String')); 
selection = contents{get(hObject,'Value')};

handles.plot_square_type = selection;

handles = plot_planforms( handles );

guidata( hObject, handles );
%--------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function Display_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Display_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function angle_noise_slider_Callback(hObject, eventdata, handles)
% hObject    handle to angle_noise_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_val = get(hObject,'Value');
if new_val > get(hObject, 'Max')
    new_val = get(hObject, 'Max');
end
if new_val < get(hObject, 'Min');
    new_val = get(hObject', 'Min');
end
handles.angle_noise = new_val;
set( handles.angle_noise_edit, 'String', num2str( new_val ) );

handles = compute_dependencies( handles );
handles = make_planforms( handles );
handles = plot_planforms( handles );
guidata(hObject, handles);
%-------------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function angle_noise_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_noise_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function angle_noise_edit_Callback(hObject, eventdata, handles)
% hObject    handle to angle_noise_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_noise_edit as text
%        str2double(get(hObject,'String')) returns contents of angle_noise_edit as a double


% --- Executes during object creation, after setting all properties.
function angle_noise_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_noise_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in new_inst_pushbutton.
function new_inst_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to new_inst_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = make_planforms( handles );
handles = plot_planforms( handles );
%--------------------------------------------------------------------------


% --- Executes on button press in base_angle_noise_checkbox.
function base_angle_noise_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to base_angle_noise_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of
% base_angle_noise_checkbox


function base_angle_offset_edit_Callback(hObject, eventdata, handles)
% hObject    handle to base_angle_offset_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_val = eval( get(hObject,'String') );
% if new_val > get(hObject, 'Max')
%     new_val = get(hObject, 'Max');
% end
% if new_val < get(hObject, 'Min');
%     new_val = get(hObject', 'Min');
% end
handles.base_angle_offset = new_val;

handles = compute_dependencies( handles );
handles = make_planforms( handles );
handles = plot_planforms( handles );
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function base_angle_offset_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to base_angle_offset_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function g = make_gaussian_mask( handles )

[x, y] = meshgrid( linspace(-handles.img_pix/2, handles.img_pix/2, handles.img_pix), linspace(-handles.img_pix/2, handles.img_pix/2, handles.img_pix) );

g = exp(-((x .^ 2) + (y .^ 2)) / (handles.gaussian_space_constant ^ 2));
return

	 



function gaussian_mask_edit_Callback(hObject, eventdata, handles)
% hObject    handle to gaussian_mask_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_val = str2double( get(hObject,'String') );
% if new_val > get(hObject, 'Max')
%     new_val = get(hObject, 'Max');
% end
% if new_val < get(hObject, 'Min');
%     new_val = get(hObject', 'Min');
% end
handles.gaussian_space_constant = new_val;
handles = compute_dependencies( handles );
handles = make_planforms( handles );
handles = plot_planforms( handles );
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function gaussian_mask_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gaussian_mask_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sqsupsq_phase_edit_Callback(hObject, eventdata, handles)
% hObject    handle to sqsupsq_phase_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_val = eval( get(hObject,'String') );
handles.sqsupsq_phase_rad = new_val;

set( handles.sqsupsq_slider, 'Value', new_val );

handles = compute_dependencies( handles );
handles = make_planforms( handles );
handles = plot_planforms( handles );
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sqsupsq_phase_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sqsupsq_phase_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Class_4_6_popupmenu.
function Class_4_6_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Class_4_6_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String')); 
selection = contents{get(hObject,'Value')};

handles.class_type = selection;
handles = compute_dependencies( handles );
handles = make_planforms( handles );
handles = plot_planforms( handles );
guidata(hObject, handles);
%-------------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function Class_4_6_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Class_4_6_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Print2File_pushbutton.
function Print2File_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Print2File_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% curr_dir = pwd();
% outptfn_prefix = [ datestr( now, 'yymmdd_HHMM') '_'];
% outptfn = [ outptfn_prefix '_' handles.plot_square_type '_' num2str( handles.angle_rad ) '_' handles.class_type '.pdf' ];
% 
% newfigure = figure; 
% newaxes = copyobj(handles.axes1,newfigure); 
% print( newfigure, '-dpdf', outptfn );
% close(newfigure) 
% fprintf('Saved %s to %s\n', outptfn, curr_dir );



function pair_angle_offset_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pair_angle_offset_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_val = eval( get(hObject,'String') );
handles.pair_angle_offset = new_val;
handles = compute_dependencies( handles );
handles = make_planforms( handles );
handles = plot_planforms( handles );
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pair_angle_offset_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pair_angle_offset_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in printable_fig_pushbutton.
function printable_fig_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to printable_fig_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotfig_h = figure;
handles.plotaxis_h = gca;
handles = compute_dependencies( handles );
handles = make_planforms( handles );
handles = plot_planforms( handles );
guidata(hObject, handles);


% --- Executes on slider movement.
function theta_slider_Callback(hObject, eventdata, handles)
% hObject    handle to theta_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_val = get(hObject,'Value');
if new_val > get(hObject, 'Max')
    new_val = get(hObject, 'Max');
end
if new_val < get(hObject, 'Min');
    new_val = get(hObject', 'Min');
end

handles.angle_rad = new_val;

% Update alphabeta edit box to show new beta/alpha ratio
alphabeta = [ num2str( sin( new_val ) ) '/' num2str( cos( new_val ) ) ];
set( handles.alphabeta_edit, 'String', alphabeta );

% Update theat edit box
set( handles.theta_edit, 'String', num2str( handles.angle_rad ) );

% set( handles.angle_noise_edit, 'String', num2str( new_val ) );

handles = compute_dependencies( handles );
handles = make_planforms( handles );
handles = plot_planforms( handles );
guidata(hObject, handles);
%-------------------------------------------------------------------------





% --- Executes during object creation, after setting all properties.
function theta_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function alphabeta_edit_Callback(hObject, eventdata, handles)
% hObject    handle to alphabeta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphabeta_edit as text
%        str2double(get(hObject,'String')) returns contents of alphabeta_edit as a double
new_val = eval( get(hObject,'String') );

handles.angle_rad = atan(new_val);
set( handles.theta_slider, 'Value', handles.angle_rad );
set( handles.theta_edit, 'String', handles.angle_rad );

handles = compute_dependencies( handles );
handles = make_planforms( handles );
handles = plot_planforms( handles );
guidata(hObject, handles);
%-------------------------------------------------------------------------



% --- Executes during object creation, after setting all properties.
function alphabeta_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphabeta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function theta_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pair_angle_offset_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function theta_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphabeta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sqsupsq_slider_Callback(hObject, eventdata, handles)
% hObject    handle to sqsupsq_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_val = get(hObject,'Value');
if new_val > get(hObject, 'Max')
    new_val = get(hObject, 'Max');
end
if new_val < get(hObject, 'Min');
    new_val = get(hObject', 'Min');
end

handles.sqsupsq_phase_rad = new_val;

set( handles.sqsupsq_phase_edit, 'String', num2str(new_val ) );

handles = compute_dependencies( handles );
handles = make_planforms( handles );
handles = plot_planforms( handles );
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sqsupsq_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sqsupsq_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
