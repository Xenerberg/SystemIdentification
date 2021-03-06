%% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\Users\Hrishik\Documents\MATLAB\data_trajectory_T6_45\Results_EPOS_IC_H.dat
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2017/03/07 10:32:17

%% Initialize variables.
filename = '\data_trajectory_T6_45\Results_EPOS_IC_H.dat';
t_0 = 200;
t_f = 400;
%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
%   column29: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
t = dataArray{:, 1};
q_0 = dataArray{:, 2};
q_1 = dataArray{:, 3};
q_2 = dataArray{:, 4};
q_3 = dataArray{:, 5};
om_x = dataArray{:, 6};
om_y = dataArray{:, 7};
om_z = dataArray{:, 8};
rho_t_x = dataArray{:, 9};
rho_t_y = dataArray{:, 10};
rho_t_z = dataArray{:, 11};
r_x = dataArray{:, 12};
r_y = dataArray{:, 13};
r_z = dataArray{:, 14};
r_dot_x = dataArray{:, 15};
r_dot_y = dataArray{:, 16};
r_dot_z = dataArray{:, 17};
om_b_x = dataArray{:, 18};
om_b_y = dataArray{:, 19};
om_b_z = dataArray{:, 20};
r_11 = dataArray{:, 21};
r_12 = dataArray{:, 22};
r_13 = dataArray{:, 23};
r_21 = dataArray{:, 24};
r_22 = dataArray{:, 25};
r_23 = dataArray{:, 26};
r_31 = dataArray{:, 27};
r_32 = dataArray{:, 28};
r_33 = dataArray{:, 29};

q_tag = [circshift(q_1,195),circshift(q_2,195),circshift(q_3,195),circshift(q_0,195)];
q_tag = q_tag(find(t==t_0):find(t==t_f),:);
save('data_trajectory_T6_45\q_tag.mat','q_tag');

%% Clear temporary variables
clearvars filename formatSpec fileID dataArray ans;