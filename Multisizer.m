%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to process Multisizer data which is exported as a
% .#m4 file. It assumes that the blank and sample measurement have the same 
% name except that the blank ends with an odd number and the sample measurement
% ends with an even number: blank BSA1_16112__01.#m4, sample BSA1_16112__02.#m4
% 
% It will output the size histogram as well as raw histogram data,
% concentration, mean size, median size, and gas volume fraction. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Setup
%1.1 Paths
basedir = 'C:\Users\ball5790\OneDrive - Nexus365\2022\November 2022\Multisizer\November 16 2022'; % top level location of data, to speed selection
[file,path] = uigetfile('*.#m4','Select file(s)','Multiselect','on',basedir);

% 1.2 Load first file
if ischar(file) 
    disp('Both the sample and blank files need to be selected')% one file was selected
    return
else
    fid = fopen(fullfile(path,file{1}));                       % multiple files were selected
end

% 1.3 Initialise variables
index = 0;
pattern = '[-+]?\d*\.\d+|\d+';  % To extract number from file
vol_sample = [];                % mL (volume of sample added to cuvette)
vol_cuvette = [];               % mL (volume of the cuvette)
vol_analysed = [];              % mL (total volume sent through the aperture)

% 1.4 Read lines from the first file to get parameters
while ~feof(fid)
    index = index + 1;
    tline = fgetl(fid);
        
    if contains(tline, 'SampleVol=')
        vol_sample = str2num(cell2mat(regexp(tline, pattern, 'match')));  
        
    elseif contains(tline, 'ElectroVol=')
        vol_cuvette = str2num(cell2mat(regexp(tline, pattern, 'match'))); 
        
    elseif contains(tline, 'AnalyticVol=')
        vol_analysed = str2num(cell2mat(regexp(tline, pattern, 'match')))/1000;
    
    end
end

% Close the file
fclose(fid);

%% 2. Extract data
% 2.1 Initialise variables
section = '';
bins = [];                      % bin width
n = [];                         % bin height

% 2.2 Get bin widths and heights
for i = 1:length(file)
    fid = fopen(fullfile(path,file{i}));
    tline = fgetl(fid);
    
    while ischar(tline)
        tline = fgetl(fid);
        if strcmp(tline([1 end]),'[]') % a 'section' starts and ends with []
            section = tline;
            index = 1;
            continue
        end
        
        switch section
            case '[#Bindiam]'
                bins(index,i)=str2num(tline);
            case '[#Binheight]'
                n(index,i)=str2num(tline);
        end
        index = index+1;
    end
    fclose(fid);
end

%% 3. Subtract blanks from sample measurements
% 3.1 Initialise variables
[row,col] = size(n); %col is the number of files selected
n_subtracted = zeros(length(bins),col/2);
index = 1;

% 3.2 Subtract blank from sample
for i = 1:2:col
    n_subtracted(:,index) = n(:,i+1)-n(:,i);
    index = index + 1;
end

% 3.3 Set negative values equal to 0
n_subtracted(n_subtracted<0) = 0; 

% 3.4 Find the mean bin heights of technical replicates
n_subtracted = round(mean(n_subtracted,2));

%% 4. Calculate stats
% 4.1 Particle concentration
n_tot = sum(n_subtracted);                              % total number of particles
conc = (n_tot/vol_analysed)*(vol_cuvette/vol_sample);   % concentration

%4.2 Mean and median size of particles 
sizes = [];                   

for i = 1:length(bins)
    sizes(end+1:end+n_subtracted(i),1) = bins(i);
end

mean_diam = mean(sizes);
median_diam = median(sizes);

% 4.3 Gas volume fraction 
gas_volume_analysed = sum((4/3)*pi.*(sizes./2).^3)*1e-9;                    % uL
gas_volume = (vol_cuvette/vol_sample)*(gas_volume_analysed/vol_analysed);   % uL/mL
 

%% 5. Save data
% % 5.1 Plot histogram
% figure(1);
% hold on;
% histogram('BinEdges',bins(:,1)','BinCounts',n_subtracted(1:end-1)');
% set(gca, 'XScale', 'log')
% xL=xlim;
% yL=ylim;
% 
% fig_name = strcat(file{1}(1:end-6),'.png');
% saveas(gcf, fullfile(path,fig_name), 'png');

% 5.2 Export raw data 
data_array = {bins(:,1), n_subtracted, conc, mean_diam, median_diam, gas_volume}; 
data_table = cell2table(data_array, 'VariableNames', {'Diameter', 'Counts','Concentration', 'Mean Diameter', 'Median Diameter', 'Gas Volume'});
file_name = strcat(file{1}(1:end-6),'.csv');
writetable(data_table, fullfile(path,file_name));