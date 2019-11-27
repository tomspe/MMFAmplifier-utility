function varargout = GainInMMF_sim2(varargin)
% GainInMMF_sim2 MATLAB code for GainInMMF_sim2.fig
%      GainInMMF_sim2, by itself, creates a new GainInMMF_sim2 or raises the existing
%      singleton*.
%
%      H = GainInMMF_sim2 returns the handle to a new GainInMMF_sim2 or the handle to
%      the existing singleton*.
%
%      GainInMMF_sim2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GainInMMF_sim2.M with the given input arguments.
%
%      GainInMMF_sim2('Property','Value',...) creates a new GainInMMF_sim2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GainInMMF_sim2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GainInMMF_sim2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GainInMMF_sim2

% Last Modified by GUIDE v2.5 02-Sep-2019 18:07:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GainInMMF_sim2_OpeningFcn, ...
                   'gui_OutputFcn',  @GainInMMF_sim2_OutputFcn, ...
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


% --- Executes just before GainInMMF_sim2 is made visible.
function GainInMMF_sim2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GainInMMF_sim2 (see VARARGIN)

Thumb1 = imread('ModeThumbnail1.jpg') ;
Thumb2 = imread('ModeThumbnail2.jpg') ;
Thumb3 = imread('ModeThumbnail3.jpg') ;
Thumb_l = imread('ModeThumbnail L.jpg') ;
Thumb_m = imread('ModeThumbnail M.jpg') ;
imshow(Thumb1, 'parent', handles.ModeThumbnail1_axes) ;
imshow(Thumb2, 'parent', handles.ModeThumbnail2_axes) ;
imshow(Thumb3, 'parent', handles.ModeThumbnail3_axes) ;
imshow(Thumb_l, 'parent', handles.IndexThumbnail1_axes) ;
imshow(Thumb_m, 'parent', handles.IndexThumbnail2_axes) ;


% Choose default command line output for GainInMMF_sim2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GainInMMF_sim2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GainInMMF_sim2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes on button press in Find_u_values_pb.
function Find_u_values_pb_Callback(hObject, eventdata, handles)
% hObject    handle to Find_u_values_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GainMedia = cellstr(get(handles.GainMediumChoice_lb,'String')) ;
handles.GainMediumString = GainMedia{get(handles.GainMediumChoice_lb,'Value')} ;
lambda_s = 1e-9*str2double(get(handles.lambda_s_et, 'String')) ;
lambda_p = 1e-9*str2double(get(handles.lambda_p_et, 'String')) ;
NA = str2double(get(handles.NA_et, 'String')) ;
a = 1e-6*str2double(get(handles.CoreRadius_et, 'String')) ;
n_core = str2double(get(handles.n_core_et, 'String')) ;

n_clad = sqrt( n_core^2 - NA^2 ) ; V_s = 2*pi*a*NA/lambda_s  ; V_p = 2*pi*a*NA/lambda_p ;
    % % u values and propagation constants - following textbook by Buck
    % % (approximation)
    
u_matrix_s = Find_u_Solutions_by_Mitschke(V_s, 1e-14)   %handles.u_matrix_s = Find_u_Solutions_by_Buck(V_s) ;
u_matrix_p = Find_u_Solutions_by_Mitschke(V_p, 1e-14)     %handles.u_matrix_p = Find_u_Solutions_by_Buck(V_p) ;
beta_matrix_s = CalculateBetaFor_u_Solutions(u_matrix_s, n_core, a, lambda_s) ;
beta_matrix_p = CalculateBetaFor_u_Solutions(u_matrix_p, n_core, a, lambda_p) ;
handles.Betas_s = ArrangeModeMatrixAsList(beta_matrix_s, u_matrix_s) ;
handles.Betas_p = ArrangeModeMatrixAsList(beta_matrix_p, u_matrix_p) ;
DrawModeNumberingTable(u_matrix_p, u_matrix_s) ;
handles.NumOfModes_s = 2*length(find(u_matrix_s)) - size(u_matrix_s, 2) ;
handles.NumOfModes_p = 2*length(find(u_matrix_p)) - size(u_matrix_p, 2) ;

handles.lambda_s = lambda_s ; handles.lambda_p = lambda_p ; handles.NA = NA ; handles.a = a ; handles.n_core = n_core ; handles.n_clad = n_clad ; handles.V_s = V_s ; handles.V_p = V_p ;
handles.u_matrix_s = u_matrix_s ; handles.u_matrix_p = u_matrix_p ;


%% Default values for the modal composition of the pump
%  for the amplitudes:
AmpMat_p = ones(size(handles.u_matrix_p)) ;
%AmpMat_p(1,1) = 1 ;   % default is to input just the first LP01 mode
handles.AmplitudeTable_p_set1 = WriteAmplitudesToCellArray(AmpMat_p, handles.u_matrix_p) ;   % for the degenerqcy set with cosines
AmpMat_p(1,1) = 0 ;     % LP01 has no degeneracy anyway...
handles.AmplitudeTable_p_set2 = WriteAmplitudesToCellArray(AmpMat_p, handles.u_matrix_p) ;  % for the degenerqcy set with sines
% and display this in the GUI (default is the amplitudes rather than the phases, for the cos set)
set(handles.AmpOrPhase_rb, 'Value', 1) ;
set(handles.CosOrSinModeSet_rb, 'Value', 1) ;
set(handles.PumpModes_table, 'Data', handles.AmplitudeTable_p_set1) ;
set(handles.AmpOrPhase_rb, 'String',  'Toggle component: displaying Amplitudes') ;
set(handles.CosOrSinModeSet_rb, 'String',  'Toggle theta degeneracy: displaying cos set') ;
%  for the phases:
PhaseMat_p = zeros(size(handles.u_matrix_p)) ;      % all phases are 0 as default
handles.PhaseTable_p_set1 = WriteAmplitudesToCellArray(PhaseMat_p, handles.u_matrix_p) ;
handles.PhaseTable_p_set2 = WriteAmplitudesToCellArray(PhaseMat_p, handles.u_matrix_p) ;
% enabling buttons of next steps
set(handles.PumpTableMark_pb, 'Enable', 'on') ; 
set(handles.PumpTableClear_pb, 'Enable', 'on') ; 
set(handles.PumpTableRandomize_pb, 'Enable', 'on') ; 
set(handles.AmpOrPhase_rb, 'Enable', 'on') ; 
set(handles.CosOrSinModeSet_rb, 'Enable', 'on') ; 
set(handles.GeneratePump_pb, 'Enable', 'on') ;
set(handles.LoadPumpConfig_pb, 'Enable', 'on') ;


guidata(hObject, handles);


% --- Executes on button press in GeneratePump_pb.
function GeneratePump_pb_Callback(hObject, eventdata, handles)
% hObject    handle to GeneratePump_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
        % the pump total power is normalized to 1
        v1 = ReadTable(handles.AmplitudeTable_p_set1, handles.u_matrix_p, 0) ;   % degeneracy group of cosines
        v2 = ReadTable(handles.AmplitudeTable_p_set2, handles.u_matrix_p, 1) ;   % degeneracy group of sines
        v = [ v1 v2 ] ;
        handles.List_Amplitudes_p = v ;
        % the pump phases:
        v3 = (pi/180)*ReadTable(handles.PhaseTable_p_set1, handles.u_matrix_p, 0) ;   
        v4 = (pi/180)*ReadTable(handles.PhaseTable_p_set2, handles.u_matrix_p, 1) ; 
        handles.List_Phases_p = [ v3 v4 ] ;

%                         %load('PumpConfigFromHugo_element.10.13.mat') ;
%                         load('D:\Users\Tom\Documents\MATLAB\MMFA sim for cluster\HugoResults\Element15-30\Simulation2\PumpConfigAsList.mat') ;
%                         handles.List_Amplitudes_p = abs(PumpConfigFromHugo) ;
%                         handles.List_Phases_p = angle(PumpConfigFromHugo) ;


% This factor sets the finesse with which we sample the z axis in the
% propagation of the speckle. If the factor is N, then we set dz to be such
% that the speckle dephases no more than pi/N between slices
handles.DephasingSamplingFactorAlongZ = str2double(get(handles.DephasingSamplingFactor_et, 'String')) ;
% Calculating the optimal dz for propagation:
handles.dz = 1e6*CalculateOptimal_dz(handles.Betas_p, handles.List_Amplitudes_p, handles.DephasingSamplingFactorAlongZ) ;

% These factors set the excursion outside the core area with which we place our spatial grid over the fiber crossection.
% The big excursion is for larger grid used only in the intialization stage where the mode profiles are normalized.
% The small excursion is for the grid over which the bulk of the calculations (all throughout the length of the fiber) are performed
NumOfPointsLargerAxis = str2double(get(handles.NumOfSpatialPointsForModeNormalization_et, 'String')) ;    % 360 ;
handles.NumOfSpatialPoints = str2double(get(handles.NumOfSpatialPoints_et, 'String')) ;
NumOfPointsInsideCore = str2double(get(handles.NumOfPointsInsideCore_et, 'String')) ;
%handles.BeyondCoreFactor_Large = 5/3 ;
%handles.BeyondCoreFactor_Main = 1.05 ;
CropMar = (NumOfPointsLargerAxis - handles.NumOfSpatialPoints)/2 ; 
dx = 2/(NumOfPointsInsideCore - 1) ;   % because the diameter is 2 in pure numbers (radius of 1)
XaxisLarger = -((NumOfPointsLargerAxis - 1)/2)*dx :dx:((NumOfPointsLargerAxis - 1)/2)*dx ;
handles.Xaxis = XaxisLarger((CropMar+1):(end-CropMar)) ;
handles.Yaxis = handles.Xaxis ;

AreaUsedForNormalizing = (2*handles.a*XaxisLarger(end))^2 ; % correction - multiply by 2 (diameter not radius)  %the length of the side of this square is the fiber diameter multiplied by the enlargenment factor 
handles.dA = AreaUsedForNormalizing/NumOfPointsLargerAxis/NumOfPointsLargerAxis ;  % in m^2

% Preparing the 3D arrays containing the spatial profiles of all exisitng
% modes for later use
handles.FieldsStack_s  = GenerateStackOfModeFields(handles.u_matrix_s, handles.V_s, XaxisLarger, CropMar) ;
handles.FieldsStack_p = GenerateStackOfModeFields(handles.u_matrix_p, handles.V_p, XaxisLarger, CropMar) ;

% for the signal, pre-calculate all multiplications of mode(i)*mode(j)
handles.SignalCouplingArray = GenerateModeOverlapArray(handles.u_matrix_s, handles.FieldsStack_s) ;


% the summed field can be found
Sx = length(handles.Xaxis) ;
Sy = length(handles.Yaxis) ;
ComplexAmplitudes = (handles.List_Amplitudes_p).*exp(1i*handles.List_Phases_p) ;
ComplexAmplitudeStack = DuplicateValuesInListToStack(ComplexAmplitudes, Sy, Sx) ;
BetaStack = DuplicateValuesInListToStack(handles.Betas_p, Sy, Sx) ;
SummedField = SuperposeModeFields(handles.FieldsStack_p, BetaStack, ComplexAmplitudeStack, 0) ;
% display and enabling buttons for next steps
imagesc((abs(SummedField)).^2, 'parent', handles.InputPump_axes) ;  
set(handles.InputPump_axes, 'Visible', 'off') ;
set(handles.GeneratePump_pb, 'Enable', 'off') ;
set(handles.InputPumpAnimation_pb, 'Visible', 'on') ;
handles.PumpAnimationPosition = 0 ;
set(handles.CalculateTM_pb, 'Enable', 'on') ;

guidata(hObject, handles);


% --- Executes on button press in RecommendPumpPin_pb.
function RecommendPumpPin_pb_Callback(hObject, eventdata, handles)
% hObject    handle to RecommendPumpPin_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[  Zaxis, DopingFactor, ~ ] = SetUpFiberLengthAndGain(handles) ;

switch handles.GainMediumString
    case 'Erbium'
        [ PumpAbsorptionCoefficient, PumpSaturationIntensity, ~, ~, ~ ] = ErbiumAbsorptionAndGainModel(DopingFactor) ;
    case 'Ytterbium'
        [ PumpAbsorptionCoefficient, PumpSaturationIntensity, ~, ~, ~ ] = YtterbiumAbsorptionAndGainModel(DopingFactor) ;
end
alpha = PumpAbsorptionCoefficient ;   %in 1/m

IntensityAlongFiber = Propagator1DForIntensityInReverse(Zaxis, alpha, PumpSaturationIntensity, 0.5*PumpSaturationIntensity) ;
OutputPowerAtCenter = handles.dA*IntensityAlongFiber(end) ;

set(handles.PumpPower_et, 'String', num2str(3.57*OutputPowerAtCenter)) ;

guidata(hObject, handles);



%%----------
%% All functions actually running some calculation of propagation along Z:



% --- Executes on button press in CalculateTM_pb.
function CalculateTM_pb_Callback(hObject, eventdata, handles)
% hObject    handle to CalculateTM_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% retrieve the fiber length and gain parameters, and also the pump configuration (truly normalized to requested power)
[  Zaxis, DopingFactor, PumpPower ] = SetUpFiberLengthAndGain(handles) ;
% normalizing the complex amplitudes to have the requested total power
TotalAmpsSquared = sum((abs(handles.List_Amplitudes_p)).^2) ;    
InjectedPumpAmplitudesInSqrtW = sqrt(PumpPower/TotalAmpsSquared)*(handles.List_Amplitudes_p).*exp(1i*handles.List_Phases_p) ;
% run pump through the fiber
[ PumpAmplitudesAlongZ, Level2PopulationThroughoutVolume ] = PumpPropagator_Corrected(handles.FieldsStack_p, handles.Betas_p, InjectedPumpAmplitudesInSqrtW, ...
            handles.Xaxis, Zaxis, handles.dA, handles.GainMediumString, DopingFactor, 1) ;
% Pump runs "backward" so the output is at the 1st column:
Pout = sum(PumpAmplitudesAlongZ(:,1).*conj(PumpAmplitudesAlongZ(:,1))) ;
Pin = sum(PumpAmplitudesAlongZ(:,end).*conj(PumpAmplitudesAlongZ(:,end))) ;
disp([ 'Pump absorption was ' num2str(10*log10(Pin/Pout)) 'dB' ])

[ TransmissionMatrix ] = CalculateTransmissionMatrix_Corrected(handles.SignalCouplingArray, handles.Betas_s, ...
                        Level2PopulationThroughoutVolume, handles.Xaxis, Zaxis, handles.GainMediumString, DopingFactor, 1) ;

zman = now ;
eval([ 'TM_'  datestr(zman, 30) ' = TransmissionMatrix ;' ]) ;
save([ 'TM_'  datestr(zman, 30) ], [ 'TM_'  datestr(zman, 30) ]) ;
figure('color', 'w') ;
subplot(1,2,1) ; imagesc(abs(TransmissionMatrix)) ; colorbar ; title('TM absolute values')
subplot(1,2,2) ; imagesc(angle(TransmissionMatrix)) ; colorbar ; title('TM phase values')

handles.CurrentTM = TransmissionMatrix ;
guidata(hObject, handles) ;
                             



%%%---------




%%----------
%% GUI functions:







% --- Executes on button press in InputPumpAnimation_pb.
function InputPumpAnimation_pb_Callback(hObject, eventdata, handles)
% hObject    handle to InputPumpAnimation_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MakeMovie = 0 ;


if strcmp(get(hObject, 'String'), 'Play')
    handles.PlayAnimation = true ;
    set(hObject, 'String', 'Stop') ;
else 
    handles.PlayAnimation = false ;
    set(hObject, 'String', 'Play') ;
end
guidata(hObject, handles) ;

if handles.PlayAnimation
    Sx = length(handles.Xaxis) ;
    Sy = length(handles.Yaxis) ;
    ComplexAmplitudes = (handles.List_Amplitudes_p).*exp(1i*handles.List_Phases_p) ;
    ComplexAmplitudeStack = DuplicateValuesInListToStack(ComplexAmplitudes,  Sy, Sx) ;
    BetaStack = DuplicateValuesInListToStack(handles.Betas_p,  Sy, Sx) ;
    c = 0 ;
    z_in_um = handles.PumpAnimationPosition ;
end

if MakeMovie
    dz = 7.5 ;  % um
    outputVideo = VideoWriter([ 'Movie' datestr(now, 30) '.avi' ]) ;
    open(outputVideo) ;
else
    dz = handles.dz ;
end
while handles.PlayAnimation   &&   ( ~MakeMovie   ||   (z_in_um-handles.PumpAnimationPosition) < 2000 )
    z_in_um = handles.PumpAnimationPosition + c*dz ;
    SummedField = SuperposeModeFields(handles.FieldsStack_p, BetaStack, ComplexAmplitudeStack, z_in_um) ;
    imagesc((abs(SummedField)).^2, 'parent', handles.InputPump_axes) ;
    set(handles.InputPump_axes, 'Visible', 'off') ;
    c = c + 1 ;
    pause(5e-3*dz) ;   %we normalize the movie animation speed to the length of our dz increments
    if strcmp(get(hObject, 'String'), 'Play') 
        handles.PumpAnimationPosition = z_in_um ;
        break ;
    end
    
    if MakeMovie
        tmpf=getframe(handles.InputPump_axes) ;
        %image(tmpf.cdata);
        imwrite(tmpf.cdata,'tmpim.jpg');
        img = imread('tmpim.jpg') ;
        writeVideo(outputVideo,img) ;
    end
end

if MakeMovie
    close(outputVideo)
end

guidata(hObject, handles) ;



%%%---------


%%%------------
% --- Executes when entered data in editable cell(s) in PumpModes_table.
function handles = PumpModes_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to PumpModes_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
set(handles.GeneratePump_pb, 'Enable', 'on') ;
set(handles.InputPumpAnimation_pb, 'Visible', 'off') ;
cla(handles.InputPump_axes) ;
% in case pump has been generated before, we have fields which are lists
% for amplitudes and phases, should be removed
if isfield(handles, 'List_Amplitudes_p') 
    handles = rmfield(handles, 'List_Amplitudes_p') ;
    handles = rmfield(handles, 'List_Phases_p') ;
end
set(handles.CalculateTM_pb, 'Enable', 'off') ;


if get(handles.AmpOrPhase_rb, 'Value') == 1
    if get(handles.CosOrSinModeSet_rb, 'Value') == 1
        handles.AmplitudeTable_p_set1 = get(handles.PumpModes_table, 'Data') ;
    else
        handles.AmplitudeTable_p_set2 = get(handles.PumpModes_table, 'Data') ;
    end
else
    if get(handles.CosOrSinModeSet_rb, 'Value') == 1
        handles.PhaseTable_p_set1 = get(handles.PumpModes_table, 'Data') ;
    else
        handles.PhaseTable_p_set2 = get(handles.PumpModes_table, 'Data') ;
    end
end
ResetGeneratedPump(hObject, handles) ;

guidata(hObject, handles);

% --- Executes on button press in CosOrSinModeSet_rb.
function CosOrSinModeSet_rb_Callback(hObject, eventdata, handles)
% hObject    handle to CosOrSinModeSet_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CosOrSinModeSet_rb

if get(hObject, 'Value') == 1
    set(hObject, 'String', 'Toggle theta degeneracy: displaying cos set') ;
    if get(handles.AmpOrPhase_rb, 'Value') == 1
        set(handles.PumpModes_table, 'Data', handles.AmplitudeTable_p_set1) ;
    else
        set(handles.PumpModes_table, 'Data', handles.PhaseTable_p_set1) ;
    end
else
    set(hObject, 'String', 'Toggle theta degeneracy: displaying sin set') ;
    if get(handles.AmpOrPhase_rb, 'Value') == 1
        set(handles.PumpModes_table, 'Data', handles.AmplitudeTable_p_set2) ;
    else
        set(handles.PumpModes_table, 'Data', handles.PhaseTable_p_set2) ;
    end
end
guidata(hObject, handles);

% --- Executes on button press in AmpOrPhase_rb.
function AmpOrPhase_rb_Callback(hObject, eventdata, handles)
% hObject    handle to AmpOrPhase_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AmpOrPhase_rb

if get(hObject, 'Value') == 1
    set(hObject, 'String', 'Toggle component: displaying Amplitudes') ;
    if get(handles.CosOrSinModeSet_rb, 'Value') == 1
        set(handles.PumpModes_table, 'Data', handles.AmplitudeTable_p_set1) ;
    else
        set(handles.PumpModes_table, 'Data', handles.AmplitudeTable_p_set2) ;
    end

else
    set(hObject, 'String', 'Toggle component: displaying Phases[ in deg ]') ;
    if get(handles.CosOrSinModeSet_rb, 'Value') == 1
        set(handles.PumpModes_table, 'Data', handles.PhaseTable_p_set1) ;
    else
        set(handles.PumpModes_table, 'Data', handles.PhaseTable_p_set2) ;
    end

end
guidata(hObject, handles);


%%%%


function AmplitudeTable = WriteAmplitudesToCellArray(AmplitudeMatrix, ValidityMatrix)
% receive a matrix, write it to a table
% we need the validity matrix to distinguish between legit modes with 0
% Amplitude, and illegal modes
% each non-zero entry in the validity matrix represents an existing mode (this can be the u values matrix)
[ m, n ] = size(AmplitudeMatrix) ;
AmplitudeTable = cell(m, n) ;
for k = 1:m*n
    if ValidityMatrix(k) ~= 0 
        AmplitudeTable{k} = [ '   '  num2str(AmplitudeMatrix(k)) ] ;
    else
        AmplitudeTable{k} = '   --' ;
    end
end


function List = ReadTable(Table, ValidityMatrix, DiscardFirstRowFlag)
% receive a cell array, write it to a list
% we need the validity matrix to distinguish between legit modes and
% illegal ones, overwriting user input to the table
% each non-zero entry in the validity matrix represents an existing mode (this can be the u values matrix)
[ m, n ] = size(Table) ;
if DiscardFirstRowFlag
    FirstRowIndex = 2 ;
else
    FirstRowIndex = 1 ;
end
NumOfModes = length(find(ValidityMatrix(FirstRowIndex:end,:))) ;
List = zeros(1, NumOfModes) ;
c = 1 ;
for k = 1:m*n
    [ row, ~] = ind2sub([ m, n ], k) ;
    if (ValidityMatrix(k) ~= 0)   &&   (row >= FirstRowIndex)      % even if the user has written an amplitude value 
                                                                    % to a forbidden mode, we overrule and ignore it 
        List(c) = sscanf(Table{k}, '%f') ;
        c = c + 1  ;
    end
end


% --- Executes on button press in PumpTableClear_pb.
function PumpTableClear_pb_Callback(hObject, eventdata, handles)
% hObject    handle to PumpTableClear_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% we generate the correct-sized all-zeros table, write it to the GUI display
% (without caring if it's amplitude or phase), than call the CellEdit
% function so that it will deal with updating the amp/phase data in handles
Table = WriteAmplitudesToCellArray(zeros(size(handles.u_matrix_p)), handles.u_matrix_p) ;
set(handles.PumpModes_table, 'Data', Table) ;
handles = PumpModes_table_CellEditCallback(hObject, eventdata, handles) ;
ResetGeneratedPump(hObject, handles) ;
guidata(hObject, handles);


% --- Executes on button press in PumpTableMark_pb.
function PumpTableMark_pb_Callback(hObject, eventdata, handles)
% hObject    handle to PumpTableMark_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% we generate the correct-sized all-ones table, write it to the GUI display
% (without caring if it's amplitude or phase), than call the CellEdit
% function so that it will deal with updating the amp/phase data in handles
Table = WriteAmplitudesToCellArray(ones(size(handles.u_matrix_p)), handles.u_matrix_p) ;
set(handles.PumpModes_table, 'Data', Table) ;
handles = PumpModes_table_CellEditCallback(hObject, eventdata, handles) ;
ResetGeneratedPump(hObject, handles) ;
guidata(hObject, handles);


% --- Executes on button press in PumpTableRandomize_pb.
function PumpTableRandomize_pb_Callback(hObject, eventdata, handles)
% hObject    handle to PumpTableRandomize_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.AmpOrPhase_rb, 'Value')     % we are on amplitudes
    RandAmps = round(100*rand(size(handles.u_matrix_p)))/10 ;
    if ~get(handles.CosOrSinModeSet_rb, 'Value')
        RandAmps(1,:) = 0 ;
    end
    Table = WriteAmplitudesToCellArray(RandAmps, handles.u_matrix_p) ;
else            % we are on phases 
    RandPhases = round(360*rand(size(handles.u_matrix_p))) ;
    if ~get(handles.CosOrSinModeSet_rb, 'Value')
        RandPhases(1,:) = 0 ;
    end
    Table = WriteAmplitudesToCellArray(RandPhases, handles.u_matrix_p) ;
end
set(handles.PumpModes_table, 'Data', Table) ;
handles = PumpModes_table_CellEditCallback(hObject, eventdata, handles) ;
ResetGeneratedPump(hObject, handles) ;
guidata(hObject, handles);


% --- Executes on button press in LoadPumpConfig_pb.
function LoadPumpConfig_pb_Callback(hObject, eventdata, handles)
% hObject    handle to LoadPumpConfig_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,~] = uigetfile('*.mat') ;
DataStruct = load([PathName '\' FileName]) ;
fieldname = fieldnames(DataStruct) ;
fieldname = fieldname{1} ;
Data = eval(['DataStruct.' fieldname]) ;
[ h, w ] = size(Data) ;
if h == 1 && w == handles.NumOfModes_p
    PumpConfig = Data ;
elseif w == 1 && h == handles.NumOfModes_p
    PumpConfig = Data.' ;
else
    errordlg('Mismatch to the expected number of Pump Modes!', 'Error in loading pump configuration', 'modal') ;
    return
end

handles = UpdatePumpConfigIntoTables(hObject, handles, PumpConfig) ;

guidata(hObject, handles);
    

function handles = UpdatePumpConfigIntoTables(hObject, handles, PumpConfig)
N = length(find(handles.u_matrix_p)) ;
Amps = abs(PumpConfig) ;
Phases = (180/pi)*angle(PumpConfig) ;

handles.AmplitudeTable_p_set1 = WriteAmplitudesToCellArray(ArrangeListAsModeMatrix(Amps(1:N), handles.u_matrix_p, 0), handles.u_matrix_p) ;
handles.AmplitudeTable_p_set2 = WriteAmplitudesToCellArray(ArrangeListAsModeMatrix(Amps(N+1:end), handles.u_matrix_p, 1), handles.u_matrix_p) ;
handles.PhaseTable_p_set1 = WriteAmplitudesToCellArray(ArrangeListAsModeMatrix(Phases(1:N), handles.u_matrix_p, 0), handles.u_matrix_p) ;
handles.PhaseTable_p_set2 = WriteAmplitudesToCellArray(ArrangeListAsModeMatrix(Phases(N+1:end), handles.u_matrix_p, 1), handles.u_matrix_p) ;
% force the toggle buttons back to the default setting, now we know which
% of the 4 tables we should display 
set(handles.CosOrSinModeSet_rb, 'Value', 1) ;
set(handles.AmpOrPhase_rb, 'Value', 1) ;
set(handles.PumpModes_table, 'Data', handles.AmplitudeTable_p_set1) ;
ResetGeneratedPump(hObject, handles) ;

guidata(hObject, handles);


function ResetGeneratedPump(hObject, handles)

handles.PlayAnimation = false ;
set(handles.InputPumpAnimation_pb, 'String', 'Play') ;
set(handles.InputPumpAnimation_pb, 'Visible', 'Off') ;
cla(handles.InputPump_axes) ;
set(handles.GeneratePump_pb, 'Enable', 'on') ;
set(handles.CalculateTM_pb, 'Enable', 'off') ;


guidata(hObject, handles);

function ResetUtility(hObject, handles)

set(handles.PumpTableMark_pb, 'Enable', 'off') ; 
set(handles.PumpTableClear_pb, 'Enable', 'off') ; 
set(handles.PumpTableRandomize_pb, 'Enable', 'off') ; 
set(handles.AmpOrPhase_rb, 'Enable', 'off') ;
set(handles.AmpOrPhase_rb, 'String', 'Toggle component: displaying Amplitudes') ;
set(handles.CosOrSinModeSet_rb, 'Enable', 'off') ;
set(handles.CosOrSinModeSet_rb, 'String', 'Toggle theta degeneracy: displaying cos set') ;
set(handles.PumpModes_table, 'Data', []) ;

set(handles.GeneratePump_pb, 'Enable', 'off') ;
set(handles.InputPumpAnimation_pb, 'Visible', 'off') ;
cla(handles.InputPump_axes) ;
if isfield(handles, 'ListAmplitudes_p')
    handles = rmfield(handles, 'ListAmplitudes_p') ;
    handles = rmfield(handles, 'ListPhases_p') ;
end
if isfield(handles, 'PumpAnimationPosition')
    handles = rmfield(handles, 'PumpAnimationPosition') ;
end

set(handles.CalculateTM_pb, 'Enable', 'off') ;


guidata(hObject, handles);


function GainFactor_et_Callback(hObject, eventdata, handles)
% hObject    handle to GainFactor_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GainFactor_et as text
%        str2double(get(hObject,'String')) returns contents of GainFactor_et as a double

% --- Executes during object creation, after setting all properties.
function GainFactor_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GainFactor_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function FiberLength_et_Callback(hObject, eventdata, handles)
% hObject    handle to FiberLength_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FiberLength_et as text
%        str2double(get(hObject,'String')) returns contents of FiberLength_et as a double

% --- Executes during object creation, after setting all properties.

function FiberLength_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FiberLength_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%------------------------------------------------------------------------------


% --- Executes on selection change in GainMediumChoice_lb.
function GainMediumChoice_lb_Callback(hObject, eventdata, handles)
% hObject    handle to GainMediumChoice_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GainMediumChoice_lb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GainMediumChoice_lb
contents = cellstr(get(hObject,'String')) ;
switch contents{get(hObject,'Value')}
case 'Erbium'
    set(handles.lambda_s_et, 'String', '1550')
    set(handles.lambda_p_et, 'String', '980')
case 'Ytterbium'
    set(handles.lambda_s_et, 'String', '1030')
    set(handles.lambda_p_et, 'String', '976')
end
ResetUtility(hObject, handles) ;

    

% --- Executes during object creation, after setting all properties.
function GainMediumChoice_lb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GainMediumChoice_lb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function lambda_s_et_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_s_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_s_et as text
%        str2double(get(hObject,'String')) returns contents of lambda_s_et as a double
ResetUtility(hObject, handles) ;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function lambda_s_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_s_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function lambda_p_et_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_p_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_p_et as text
%        str2double(get(hObject,'String')) returns contents of lambda_p_et as a double
ResetUtility(hObject, handles) ;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function lambda_p_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_p_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function NA_et_Callback(hObject, eventdata, handles)
% hObject    handle to NA_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NA_et as text
%        str2double(get(hObject,'String')) returns contents of NA_et as a double
ResetUtility(hObject, handles) ;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function NA_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NA_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CoreRadius_et_Callback(hObject, eventdata, handles)
% hObject    handle to CoreRadius_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CoreRadius_et as text
%        str2double(get(hObject,'String')) returns contents of CoreRadius_et as a double
ResetUtility(hObject, handles) ;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function CoreRadius_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CoreRadius_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function n_core_et_Callback(hObject, eventdata, handles)
% hObject    handle to n_core_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_core_et as text
%        str2double(get(hObject,'String')) returns contents of n_core_et as a double
ResetUtility(hObject, handles) ;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function n_core_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_core_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function DephasingSamplingFactor_et_Callback(hObject, eventdata, handles)
% hObject    handle to DephasingSamplingFactor_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DephasingSamplingFactor_et as text
%        str2double(get(hObject,'String')) returns contents of DephasingSamplingFactor_et as a double

v1 = ReadTable(handles.AmplitudeTable_p_set1, handles.u_matrix_p, 0) ;   % degeneracy group of cosines
v2 = ReadTable(handles.AmplitudeTable_p_set2, handles.u_matrix_p, 1) ;   % degeneracy group of sines
v = [ v1 v2 ] ;
handles.List_Amplitudes_p = v ;
% This factor sets the finesse with which we sample the z axis in the
% propagation of the speckle. If the factor is N, then we set dz to be such
% that the speckle dephases no more than pi/N between slices
handles.DephasingSamplingFactorAlongZ = str2double(get(handles.DephasingSamplingFactor_et, 'String')) ;
% Calculating the optimal dz for propagation:
handles.dz = 1e6*CalculateOptimal_dz(handles.Betas_p, handles.List_Amplitudes_p, handles.DephasingSamplingFactorAlongZ) ;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function DephasingSamplingFactor_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DephasingSamplingFactor_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NumOfSpatialPoints_et_Callback(hObject, eventdata, handles)
% hObject    handle to NumOfSpatialPoints_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumOfSpatialPoints_et as text
%        str2double(get(hObject,'String')) returns contents of NumOfSpatialPoints_et as a double
ResetGeneratedPump(hObject, handles) ;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function NumOfSpatialPoints_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumOfSpatialPoints_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function NumOfSpatialPointsForModeNormalization_et_Callback(hObject, eventdata, handles)
% hObject    handle to NumOfSpatialPointsForModeNormalization_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumOfSpatialPointsForModeNormalization_et as text
%        str2double(get(hObject,'String')) returns contents of NumOfSpatialPointsForModeNormalization_et as a double
ResetGeneratedPump(hObject, handles) ;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function NumOfSpatialPointsForModeNormalization_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumOfSpatialPointsForModeNormalization_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function NumOfPointsInsideCore_et_Callback(hObject, eventdata, handles)
% hObject    handle to NumOfPointsInsideCore_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumOfPointsInsideCore_et as text
%        str2double(get(hObject,'String')) returns contents of NumOfPointsInsideCore_et as a double
ResetGeneratedPump(hObject, handles) ;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function NumOfPointsInsideCore_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumOfPointsInsideCore_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected cell(s) is changed in PumpModes_table.
function PumpModes_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to PumpModes_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected cell(s) is changed in SignalModes_table.
function SignalModes_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to SignalModes_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function InputPump_axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to InputPump_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  

% --- Executes on mouse press over axes background.
function InputSignal_axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to InputSignal_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




function PumpPower_et_Callback(hObject, eventdata, handles)
% hObject    handle to PumpPower_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PumpPower_et as text
%        str2double(get(hObject,'String')) returns contents of PumpPower_et as a double


% --- Executes during object creation, after setting all properties.
function PumpPower_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PumpPower_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




        
