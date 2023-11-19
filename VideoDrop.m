%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file goes through the videodrop files in a folder and saves the file
% name, median diameter, concentration, and histogram values (bin width and
% height) in a table. These values can then be used in prism or other
% software to plot figures and do statistical analysis. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Setup
clear all;

% 1.1 Paths
dataset_dir = "C:\Users\ball5790\OneDrive - Nexus365\2023\January 2023\Videodrop\January 16 2023\";
datadir = dir(dataset_dir);                                         % list of folders
datadir = datadir(~ismember({datadir.name}, {'.', '..'}));          % remove aliases for parent and current folder
datadir = datadir([datadir(:).isdir]);                              % only list folders and not files in dataset_dir

file_name_save = '16012023.csv';                                    % name of file you will save

% 1.2 Initialise struct for data values
group = struct('name',{},'median_diameter',{},'concentration',{},'hist_values',{});

for i = 1:length(datadir) 
    group(i).name = datadir(i).name;                                
end

%% 2. Loop over all folders
for i = 1:length(datadir)
    
    folder_name = datadir(i).name;                                          % folder name (i.e. BSA1_PCaN_0h_16012023)
    
    file_dir = dir(fullfile(strcat(dataset_dir,folder_name,'\'),'*.csv'));  % access the csv file (data file) in the current folder
    file_name = file_dir.name;                                              % file name
    full_file_dir = strcat(file_dir.folder,'\',file_name);                  % full file path
    
    current_group = strcat(folder_name(1:3),folder_name(5:end));
    
    [hist_edges, hist_val, hist_midpoints,median_diam] = readVD(full_file_dir);
    conc = getConc(full_file_dir);
    
    group(i).median_diameter = median_diam;
    group(i).concentration = conc;
    group(i).hist_values = hist_val;
    
end

results_table = struct2table(group);
hist_values = array2table(results_table(:,4).hist_values, 'VariableNames', hist_midpoints);
results_table = [results_table(:,1:3) hist_values];
writetable(results_table, fullfile(dataset_dir, file_name_save));

%% 3. Define functions
function [hist_edges, hist_val, hist_midpoints, median_diam] = readVD(full_file_dir)
    
    % Define histogram variables
    size_min = 0;                                                           % nm
    size_max = 1000;                                                        % nm
    size_hist_step = 50;
    hist_step = size_min:size_hist_step:size_max;
    
    first_midpoint = (size_hist_step-size_min)/2;
    hist_midpoints = first_midpoint:size_hist_step:size_max;
    hist_midpoints = arrayfun(@num2str, hist_midpoints, 'UniformOutput',false);

    % Get histogram values
    full_data = readtable(full_file_dir,'VariableNamingRule','preserve');   % read csv file
    size_table = [full_data(:,6),full_data(:,1),full_data(:,4)];             
    size_table = sortrows(size_table, 'diameter','ascend');
    size_array = table2array(size_table);
    diameter = size_array(:,3);
    [hist_val, hist_edges] = histcounts(diameter,hist_step);
    
    % Calculate median and std of particle diameter
    median_diam = median(diameter);
    std_diameter = std(diameter);
end

function [conc] = getConc(full_file_dir)
    fid = fopen(full_file_dir);
    
    for i = 1:10 % the concentration value should be at the top of the csv file
        tline = fgetl(fid);
        
        if contains(tline, 'Concentration')
            conc = str2double(tline(15:21));
        end
    end
end
