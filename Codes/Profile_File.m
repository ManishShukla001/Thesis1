% ts_000_create_profile_file.m
%
% This script takes coordinates and writes a profile .nc file with all the
% model ouput data at that location. This information is passed to various
% timeseries visualisation scripts for quick and efficent plotting
%
% Usage: 
% Input a TUFLOW FV results file and list of coordinates and names to
% generate the profile file. Please ensure that coordinated used are in the
% same projection as the TUFLOW FV model. 
%
% SDE 2018

%clear
%close 

% ------------------- User Input ------------------
tfv_res_files = {'C:\Manish\May\Simulation\2\r1\HYD_0S.nc';
                 };
Sites=readmatrix("center.xlsx") ;   

for i = 1:length(Sites)
 Names{i} = ['p' num2str(i)];
    end          
             
             
output_sites        = [ 652650	1081750 ];

%						159.0845  -31.3727
%						159.0906  -31.3814
%						159.1001  -31.3948
%						159.1154  -31.4032
%						159.1266  -31.4105
%						159.1202  -31.4165
%						159.1178  -31.4236
%						159.1277  -31.4402];
                    
site_names          = {'Point_1'}; %cell(1, 2984);%{ 'Point_1'
                      %  'Point_2'
                      %  'Point_3'
                      %  'Point_4'
                       % 'Point_5'
                       % 'Point_6'
                       % 'Point_7'
                       % 'Point_8'
                       % 'Point_9'};
%for i=1:2984
%    site_names{i}=sprintf('Point_%d', i);
%end

                    
% ----------------- Generate Profiles ----------------
for i=1:length(tfv_res_files)
    tfv_resfile = tfv_res_files{i};
    % create a meaningful output name
    out_name = strrep(tfv_resfile,'.nc','_PROFILESK.nc');

    % call the profile generation script
    fv_create_profiles(tfv_resfile,out_name,Sites,Names);
end