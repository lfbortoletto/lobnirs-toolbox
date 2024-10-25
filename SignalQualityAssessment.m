function [sqiTable] = SignalQualityAssessment(dataFolder, nirsExtension, sqiTable)
% SignalQualityAssessment Manually assess and evaluate the quality of NIRS signals across multiple channels.
%
% This function provides an interactive environment for manually reviewing the 
% signal quality of channels from Near-Infrared Spectroscopy (NIRS) data files. 
% Additionally, it helps evaluate the performance of automated signal quality 
% metrics.
%
% INPUTS:
% 1) dataFolder: (String) Path to the folder containing one or more NIRS data files. 
%                The function will search for files with the specified extension.
%
% 2) nirsExtension: (String) File extension used to identify NIRS data files 
%                   (e.g., '.lob' or '.nirs').
%
% 3) sqiTable: (Optional) A previously saved table containing the results of 
%              a past signal quality assessment. Providing this input allows 
%              you to resume from where you left off, avoiding redundant assessments.
%
% OUTPUT:
% sqiTable: A struct summarizing the manual assessment results and the evaluated 
%           signal quality metrics for each NIRS channel processed.
%
% USAGE:
% sqiTable = SignalQualityAssessment('C:\data\NIRS', '.nirs'); 
% sqiTable = SignalQualityAssessment('C:\data\NIRS', '.nirs', savedSqiTable);
%
% This function supports iterative assessments, enabling you to stop and continue 
% as needed by reloading a saved 'sqiTable'.
%
% Created by: L. F. Bortoletto (2024/10/22)
%
% -------------------------------------------------------------------------

% Get files inside folder.
[filePathList, ~] = getAllFiles(dataFolder, nirsExtension);

% Get the number of subjects in folder.
nSubjects = numel(filePathList);

% Load all nirs files.
for participant = 1:nSubjects
    nirs{participant} = load(filePathList{participant},'-mat'); 
end

% Get the number of channels in the probe.
nChannels = size(nirs{1}.data.SD.MeasList, 1);

% Initialize sqiTable if input isn't given and valid elements list.
if exist('sqiTable', 'var')
    
    % Initialize table valid elements to choose from.
    validElementsList = 1:nChannels*nSubjects;
    
    % Remove elements already processed.
    removeElement = find(~isnan(sqiTable.manual(:)));
    validElementsList(removeElement) = []; %#ok<FNDSB>
    
else
    
    % Initialize signal quality index tables.
    sqiTable.manual = nan(nChannels, nSubjects);
    sqiTable.snr = nan(nChannels, nSubjects);
    sqiTable.sci = nan(nChannels, nSubjects);
    
    % Initialize table valid elements to choose from.
    validElementsList = 1:nChannels*nSubjects;
    
end

% Initialize button_value.
button_value = 0;

while ~isnan(button_value)
        
    % Clear previous button_value and snrValue.
    clear button_value snrValue sciValue;
    
    % Choose a random element for analysis.
    currentElement = randsample(validElementsList, 1);
    
    % Find participant and channel number chosen.
    [channelIndex, subjectIndex] = ind2sub(size(sqiTable.manual), currentElement);
    
    % Remove from valid list.
    removeElement = find(validElementsList == currentElement);
    validElementsList(removeElement) = []; %#ok<FNDSB>
    
    % Get chosen time series timings and samples.
    timings = nirs{subjectIndex}.data.t;
    samples = nirs{subjectIndex}.data.d(:, channelIndex);
    
    % Compute SNR.
    snrValue = getSNR(samples);
    
    % Compute SCI.
    sciValue = getSCI(nirs{subjectIndex}.data, nChannels, channelIndex);
    
    % Compute PSD values.
    [F, I] = getPSD(timings, samples);
    
    % Display time series and ask user to define the time series quality.
    button_value = plot_and_get_button(F, I, timings, samples);
     
    % Store button value in the true table.
    sqiTable.manual(channelIndex, subjectIndex) = button_value;
    
    % Store computed snr in the snr table.
    sqiTable.snr(channelIndex, subjectIndex) = snrValue;
    
    % Store computed sci in the sci table.
    sqiTable.sci(channelIndex, subjectIndex) = sciValue;

end

end

function [filePathList, fileNameList] = getAllFiles(dirName, extension)
    % Initialize the output file lists
    filePathList = {};
    fileNameList = {};

    % Check if the extension is provided, otherwise get all files
    if nargin < 2 || isempty(extension)
        extension = '';  % No extension filtering
    end

    % List all files and folders in the current directory
    dirData = dir(dirName);

    % Filter out '.' and '..' directories
    dirIndex = [dirData.isdir];
    validNames = {dirData(dirIndex).name};
    validNames = validNames(~ismember(validNames, {'.', '..'}));

    % Recursively call getAllFiles on subdirectories
    for i = 1:length(validNames)
        nextDir = fullfile(dirName, validNames{i});
        [subFilePaths, subFileNames] = getAllFiles(nextDir, extension);
        filePathList = [filePathList; subFilePaths];
        fileNameList = [fileNameList; subFileNames];
    end

    % Get a list of file names (not directories)
    fileIndex = ~[dirData.isdir];
    fileName = {dirData(fileIndex).name};

    % Filter files based on the given extension if provided
    if ~isempty(extension)
        fileName = fileName(endsWith(upper(fileName), upper(extension)));
    end

    % Add full path to files and update file lists
    if ~isempty(fileName)
        fullFileName = cellfun(@(x) fullfile(dirName, x), ...
                               fileName, 'UniformOutput', false);
        filePathList = [filePathList; fullFileName'];
        fileNameList = [fileNameList; fileName'];
    end
end

function button_value = plot_and_get_button(X1, Y1, X2, Y2)
    % Create a fullscreen figure
    fig = figure('Name', 'Custom Plot Layout with Buttons', 'NumberTitle', 'off', ...
                 'Position', get(0, 'ScreenSize'));  % Fullscreen window

    % Create a tiled layout with 2 rows and 3 columns
    t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % First plot: occupying row 1, column 1 only
    nexttile([1, 1]);
    plot(X1, Y1, 'LineWidth', 1.5);
    xlabel('X1-axis');
    ylabel('Y1-axis');
    title('Plot 1');
    grid on;

    % Second plot: occupying row 1, columns 2 and 3
    nexttile([1, 2]);
    plot(X2, Y2, 'LineWidth', 1.5);
    xlabel('X2-axis');
    ylabel('Y2-axis');
    title('Plot 2');
    grid on;

    % Create a button panel occupying the second row
    button_panel = uipanel('Parent', fig, 'Position', [0.1, 0.2, 0.8, 0.1]);

    % Create buttons for Stop, Bad, and Good
    stop_button = uicontrol('Parent', button_panel, 'Style', 'pushbutton', ...
        'String', 'Stop', 'Units', 'normalized', 'Position', [0.05, 0.1, 0.25, 0.8], ...
        'Callback', @(src, event) buttonCallback(nan));

    bad_button = uicontrol('Parent', button_panel, 'Style', 'pushbutton', ...
        'String', 'Bad', 'Units', 'normalized', 'Position', [0.37, 0.1, 0.25, 0.8], ...
        'Callback', @(src, event) buttonCallback(0));

    good_button = uicontrol('Parent', button_panel, 'Style', 'pushbutton', ...
        'String', 'Good', 'Units', 'normalized', 'Position', [0.69, 0.1, 0.25, 0.8], ...
        'Callback', @(src, event) buttonCallback(1));

    % Initialize the button value
    button_value = [];

    % Wait until a button is pressed by blocking MATLAB execution
    uiwait(fig);

    % Callback function to store the button value and close the window
    function buttonCallback(value)
        button_value = value;  % Store the button value
        uiresume(fig);         % Resume execution
        close(fig);            % Close the figure
    end
end


function [F, I] = getPSD(X, Y)
    
% Define variables.
fSample = median(1./diff(X));
sampleSize = length(Y);

% Compute PSD through Welch method.
[Pxx,F] = pwelch(Y, sampleSize, [], [], fSample);

% Plot PSD.
I = 10*log10(Pxx);
   
end

function snrValue = getSNR(samples)

    M = mean(samples);
    STD = std(samples);
    
    snrValue = M/STD;

end

function sciValue = getSCI(nirs, nChannels, channelIndex)
    
    % Get WL 1 and WL 2 channel indices.
    if channelIndex <= nChannels/2
        indexWL_1 = channelIndex;
        indexWL_2 = channelIndex + nChannels/2;
    else
        indexWL_2 = channelIndex;
        indexWL_1 = channelIndex - nChannels/2;
    end
    
    % Define bandpass frequency range.
    bandpassInterval = [0.5, 2.5];

    % Filter data in the desired frequencies
    dFilteredWl_1 = hmrBandpassFiltLOB(nirs.d(:, indexWL_1), nirs.SD.f, bandpassInterval(1), bandpassInterval(2));
    dFilteredWl_2 = hmrBandpassFiltLOB(nirs.d(:, indexWL_2), nirs.SD.f, bandpassInterval(1), bandpassInterval(2));        
    
    % Compute SCI.
    C = corrcoef(dFilteredWl_1, dFilteredWl_2); 
    sciValue = C(1,2);    

end