clearvars;
close all;
% Add path to stlTools
% You can download the package for free from: 
% https://es.mathworks.com/matlabcentral/fileexchange/51200-stltools

%addpath('./stlTools');

% Set the name of the mat file containing all the info of the 3D model

MatFileName = 'leonardo_m346_3d_model.mat';


% Define the list of parts which will be part of the rigid aircraft body
%rigid_body_list   = {     'Body_m346.stl'  };
rigid_body_list   = { 'm346_body.stl' 'Rear_canopy.stl' 'Center_canopy.stl' 'Forward_canopy.stl' 'ail_left.stl' 'ail_right.stl' 'hst_left.stl' 'hst_right.stl' 'lef_left.stl' 'lef_right.stl' 'Rudder.stl'  };
%rigid_body_list   = { 'M346_original_rotated.stl'  };


% Define the color of each part
rigid_body_colors = {0.8 * [1, 1, 1], 0.1 * [1, 1, 1], 0.1 * [1, 1, 1], 0.1 * [1, 1, 1], 0.6 * [1, 1, 1], 0.6 * [1, 1, 1], 0.6 * [1, 1, 1], 0.6 * [1, 1, 1], 0.6 * [1, 1, 1], 0.6 * [1, 1, 1], 0.6 * [1, 1, 1]};
% Define the transparency of each part
alphas            = [1, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9];
% Define the model ofset to center the A/C Center of Gravity
offset_3d_model   = [0,0,0];
% Define the controls 
ControlsFieldNames = {...
'model'              'label',             'color',                             'rot_point',            'rot_vect', 'max_deflection'};
Controls = {                                                                                                       
%'hst_right.stl',       'HST_R',   0.3*[0.8, 0.8, 1], [47.0,    2.07,  1.78]            [-0.7, 0.3 , 0], [-30, +30];
%'hst_left.stl',        'HST_L',   0.3*[0.8, 0.8, 1], [47.0,   -2.07,  1.78],           [-0.7,-0.3, 0], [-30, +30];
% 'LE_Left.stl',        'LE_L',   0.3*[0.8, 0.8, 1], [9267, -2493.0,  380.2]-offset_3d_model,  [0.7849, -0.6197, 0], [ -1, +30];
% 'LE_right.stl',       'LE_R',   0.3*[0.8, 0.8, 1], [9267, +2493.0,  380.2]-offset_3d_model, [-0.7849, -0.7197, 0], [ -1, +30];
% 'Rudder.stl',          'RUD',   0.3*[0.8, 0.8, 1], [12930,    0.0, 1387.0]-offset_3d_model,            [0, 0, -1], [-30, +30];
% 'Elevon_Left.stl',  'FLAP_L',   0.3*[0.8, 0.8, 1], [11260, -860.9,  368.3]-offset_3d_model,  [+0.0034, 0.9999, 0], [-30, +30];
% 'Elevon_Right.stl', 'FLAP_R',   0.3*[0.8, 0.8, 1], [11260, +860.9,  368.3]-offset_3d_model,  [-0.0034, 0.9999, 0], [-30, +30];
};
% Definition of the Model3D data structure
% Rigid body parts
for i = 1:length(rigid_body_list)
    Model3D.Aircraft(i).model = rigid_body_list{i};
    Model3D.Aircraft(i).color = rigid_body_colors{i};
    Model3D.Aircraft(i).alpha = alphas(i);
    % Read the *.stl file
   [Model3D.Aircraft(i).stl_data.vertices, Model3D.Aircraft(i).stl_data.faces, ~, Model3D.Aircraft(i).label] = stlRead(rigid_body_list{i});
    Model3D.Aircraft(i).stl_data.vertices  = Model3D.Aircraft(i).stl_data.vertices - offset_3d_model;
end
% Controls parts
for i = 1:size(Controls, 1)
    for j = 1:size(Controls, 2)
        Model3D.Control(i).(ControlsFieldNames{j}) = Controls{i, j};
    end
    % Read the *.stl file
    [Model3D.Control(i).stl_data.vertices, Model3D.Control(i).stl_data.faces, ~, ~] = stlRead( Model3D.Control(i).model);
    Model3D.Control(i).stl_data.vertices = Model3D.Control(i).stl_data.vertices - offset_3d_model;
end

%% Save mat file
save(MatFileName, 'Model3D');

%% Check the results
% Get maximum dimension to plot the circles afterwards
AC_DIMENSION = max(max(sqrt(sum(Model3D.Aircraft(1).stl_data.vertices.^2,2))));
% for i=1:length(Model3D.Control)
%     AC_DIMENSION = max(AC_DIMENSION,max(max(sqrt(sum(Model3D.Control(i).stl_data.vertices.^2,2)))));
% end
% Define the figure properties
AX = axes('position',[0.0 0.0 1 1]);
axis off
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[scrsz(3)/40 scrsz(4)/12 scrsz(3)/2*1.0 scrsz(3)/2.2*1.0],'Visible','on');
set(AX,'color','none');
axis('equal')
hold on;
cameratoolbar('Show')
% Circles around the aircraft transformation group handles
euler_hgt(1)  = hgtransform('Parent',           AX, 'tag', 'OriginAxes');
euler_hgt(2)  = hgtransform('Parent', euler_hgt(1), 'tag', 'roll_disc');
euler_hgt(3)  = hgtransform('Parent', euler_hgt(1), 'tag', 'pitch_disc');
euler_hgt(4)  = hgtransform('Parent', euler_hgt(1), 'tag', 'heading_disc');
euler_hgt(5)  = hgtransform('Parent', euler_hgt(2), 'tag', 'roll_line');
euler_hgt(6)  = hgtransform('Parent', euler_hgt(3), 'tag', 'pitch_line');
euler_hgt(7)  = hgtransform('Parent', euler_hgt(4), 'tag', 'heading_line');
% Plot objects
% -------------------------------------------------------------------------
% Plot airframe
AV_hg         = hgtransform;

% controls_deflection_deg transformation group handles
% CONT_hg       = zeros(1,length(Model3D.Control));
% for i=1:length(Model3D.Control)
%     CONT_hg(i) = hgtransform('Parent', AV_hg, 'tag', Model3D.Control(i).label);
% end

for i = 1:length(Model3D.Aircraft)
    AV = patch(Model3D.Aircraft(i).stl_data, ...
        'FaceColor',        Model3D.Aircraft(i).color, ...
        'EdgeColor',        'none',        ...
        'FaceLighting',     'gouraud',     ...
        'Parent',            AV_hg,       ...
        'AmbientStrength',   0.15);
end


% AIRCRAFT BODY
 M1 = makehgtform('zrotate', 0*pi / 180);  % Heading rotation
 M2 = makehgtform('yrotate', 0*pi / 180);  % Pitch rotation
 M3 = makehgtform('xrotate', 0*pi / 180);  % bank_deg rotation
 set(AV_hg, 'Matrix',M1 * M2 * M3)


%CONT(length(Model3D.Control))=0;

% %Plot controls
% for i=1:length(Model3D.Control)
%     CONT(i) = patch(Model3D.Control(i).stl_data, ... 
%         'FaceColor',        Model3D.Control(i).color, ...
%         'EdgeColor',        'none',        ...
%         'FaceLighting',     'gouraud',     ...
%         'Parent',            AV_hg,       ...
%         'AmbientStrength',  0.15);
%     % Plot the rotation point and the rotation axis of each control
%     % (double-check correct implementation and rotation direction of each
%     % control surface)
%     p = Model3D.Control(i).rot_point;
%     vect = Model3D.Control(i).rot_vect;
%     plot3(p(1)+[0, AC_DIMENSION*vect(1)/2], p(2)+[0, AC_DIMENSION*vect(2)/2], p(3)+[0, AC_DIMENSION*vect(3)/2], 'b-o', 'MarkerSize', 10, 'LineWidth', 2);
% 
%     M1 = makehgtform('translate', -Model3D.Control(i).rot_point);   % Heading
%     M2 = makehgtform('axisrotate', Model3D.Control(i).rot_vect, 30 * pi / 180);  % Pitch
%     M3 = makehgtform('translate', Model3D.Control(i).rot_point);  % bank_deg
%     set(CONT_hg(i), 'Matrix', M3 * M2 * M1);
% 
% end




% Fixing the axes scaling and setting a nice view angle
axis('equal');
% axis([-1 1 -1 1 -1 1] * 2.0 * AC_DIMENSION)
set(gcf,'Color',[1 1 1])
axis off
camlight('left');
%camlight(20, 50)   %creates a light at AZ, EL from camera.
 
material('DULL');
view([30 10])
%zoom(0.5);
% Add a camera light, and tone down the specular highlighting


% --------------------------------------------------------------------
% Define the radius of the sphere
R = 1.4 * AC_DIMENSION;

% Outer circles
phi = (-pi:pi/36:pi)';
D1 = [sin(phi) cos(phi) zeros(size(phi))];
HP(1) = plot3(R*D1(:,1),R*D1(:,2),+R*D1(:,3),'Color','b','tag','Zplane','Parent',euler_hgt(4));
HP(2) = plot3(R*D1(:,2),R*D1(:,3),+R*D1(:,1),'Color',[0 0.8 0],'tag','Yplane','Parent',euler_hgt(3));
HP(3) = plot3(R*D1(:,3),R*D1(:,1),+R*D1(:,2),'Color','r','tag','Xplane','Parent',euler_hgt(2));

% +0,+90,+180,+270 Marks
S = 0.95;
phi = -pi+pi/2:pi/2:pi;
D1 = [sin(phi); cos(phi); zeros(size(phi))];
plot3([S*R*D1(1,:); R*D1(1,:)],[S*R*D1(2,:); R*D1(2,:)],[S*R*D1(3,:); R*D1(3,:)],'Color','b','tag','Zplane','Parent',euler_hgt(4));
plot3([S*R*D1(2,:); R*D1(2,:)],[S*R*D1(3,:); R*D1(3,:)],[S*R*D1(1,:); R*D1(1,:)],'Color',[0 0.8 0],'tag','Yplane','Parent',euler_hgt(3));
plot3([S*R*D1(3,:); R*D1(3,:)],[S*R*D1(1,:); R*D1(1,:)],[S*R*D1(2,:); R*D1(2,:)],'Color','r','tag','Xplane','Parent',euler_hgt(2));
text(R*1.05*D1(1,:),R*1.05*D1(2,:),R*1.05*D1(3,:),{'N','E','S','W'},'Fontsize',9,'color',[0 0 0],'HorizontalAlign','center','VerticalAlign','middle');

% +45,+135,+180,+225,+315 Marks
S = 0.95;
phi = -pi+pi/4:2*pi/4:pi;
D1 = [sin(phi); cos(phi); zeros(size(phi))];
plot3([S*R*D1(1,:); R*D1(1,:)],[S*R*D1(2,:); R*D1(2,:)],[S*R*D1(3,:); R*D1(3,:)],'Color','b','tag','Zplane','Parent',euler_hgt(4));
HT = text(R*1.05*D1(1,:),R*1.05*D1(2,:),R*1.05*D1(3,:),{'NW','NE','SE','SW'},'Fontsize',8,'color',[0 0 0],'HorizontalAlign','center','VerticalAlign','middle');

% 10 deg sub-division marks
S = 0.98;
phi = -[0:10:90 80:-10:0 -10:-10:-90 -80:10:0];
PHI_TEXT{length(phi)}='';
for i=1:length(phi)
    PHI_TEXT{i} = num2str(phi(i));
end
theta_t = -[0:10:90 80:-10:0 -10:-10:-90 -80:10:0];
THETA_TEXT{length(theta_t)}='';
for i=1:length(theta_t)
    THETA_TEXT{i} = num2str(theta_t(i));
end
phi = -180:10:180;
phi = phi*pi/180;
D1 = [sin(phi); cos(phi); zeros(size(phi))];
plot3([S*R*D1(1,:); R*D1(1,:)],[S*R*D1(2,:); R*D1(2,:)],[S*R*D1(3,:); R*D1(3,:)],'Color','b','tag','Zplane','Parent',euler_hgt(4));
plot3([S*R*D1(2,:); R*D1(2,:)],[S*R*D1(3,:); R*D1(3,:)],[S*R*D1(1,:); R*D1(1,:)],'Color',[0 0.8 0],'tag','Yplane','Parent',euler_hgt(3));
plot3([S*R*D1(3,:); R*D1(3,:)],[S*R*D1(1,:); R*D1(1,:)],[S*R*D1(2,:); R*D1(2,:)],'Color','r','tag','Xplane','Parent',euler_hgt(2));

% Plot guide lines
HL(1) = plot3([-R R],[0 0],[0 0],'b-','tag','heading_line','parent',euler_hgt(7));
HL(2) = plot3([-R R],[0 0],[0 0],'g-','tag','pitch_line','parent',euler_hgt(6),'color',[0 0.8 0]);
HL(3) = plot3([0 0],[-R R],[0 0],'r-','tag','roll_line','parent',euler_hgt(5));

