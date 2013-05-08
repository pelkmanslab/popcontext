function handles = MeasurePopulationContext(handles)

% Help for the Measure Population Context module:
% Category: Measurement
%
%
% SHORT DESCRIPTION:
%
% This module calculates the population context parameters (local cell
% density, cell islet edge, cell area and population size), as described by
% Snijder et al., Nature, 2009 and in Snijder & Sacher, et al, in revision,
% 2011
%
% The population context is a set of measurements that describe the spatial
% (and cellular growth) organization of your cells in their images. The
% measurements include: (1) the local cell density; (2) the population
% size; (3) the object (cell or nucleus) area; (4) wether a cell is located at the
% cell islet edge or not; (5) the distance of each cell to the closest edge
% of a cell colony.; and (6) the distance of a cell to its closest
% neighbor. 
%
% If population context is measured per well, the module will also store
% the results of the image filename parsing: (1) Well Row; (2) Well Column;
% (3) Image Site; (4) Microscope type; (5) Timepoint; and (6) the Image
% Snake. 
%
%
% *************************************************************************
% Here is more detailed information on the different options:
%
% Which object would you like to use for the cell population context
% measurement?
%   --> Select the object to be used for population context analysis.
%   Typically a nucleus or cell object would be used.
%
% What is the typical pixel diameter of your cells?
%   --> This parameter should reflect the typical diameter (approximately)
%   of cells. It determines the size of the smoothing function applied for
%   the local cell density calculation, and the distance of object
%   expansion for detection of edge cells. The default value of 140 works
%   for typical (HeLa) cells in 10x images that can be detected with a
%   nuclear diameter between 10 and 40 pixels (default settings of
%   IdentifyPrimAutomatic).
%   An optimal setting should result in a smoothed LCD figure with an
%   average of above 1 if you see different densities in your image.
%   Additionally if your image clearly contains edge and non-edge
%   cells the % of edge cells should not be either 0% or 100%. If you get
%   0% edge cells, lower the value of this setting, if you get 100% edge
%   cells, you'd want to increase it.  
%   The Nearest Neighbor Distance calculation gives a suggestion for what
%   might be a reasonable cell diameter setting for your data, which is
%   shown in the lower right panel in the display figure. 
%   Depending on your cell type, approximately 5 x your typical nucleus
%   diameter might be a good first estimate. 
%
% How many times do you want to shrink to speed up the analysis?
% 'Automatic' or enter a value not smaller than 1. 
%   --> Leave 'Automatic' if the function is running fast enough. Automatic
%   shrinks your images 2x if you analyze pop.context per image, and
%   shrinks it such that one dimension of your combined well image is still
%   larger than 1000px. Value should be an integer equal to or larger than
%   1. 
%
% Do you want to measure the population context per image or per well? Per
% well only works for "known" filenames.  
%   --> Population context measurements can be made per image or per well.
%   If you have aquired multiple adjacent images per well, the module can
%   try to infer from your filenames the well and site information, and
%   virtually stitch back together the individual imaged sites. Overall
%   this increases the accuracy of the population context measurements.
%   Ideally primary object detection is set to discard objects that touch
%   the image edge, and the stitching algorithm will slightly overlap
%   neighboring image sites to overcome the negative bias due to the object
%   discarding. See below for examples of known file names.
%
% Per well: What is the amount of overlap in pixels between images of the
% same well? Enter 'Automatic' for automatic detection. 
%   --> Leave 'Automatic' if you're analyzing popcontext per well, and if
%   the algorithm finds a well defined optimal overlap value at a minimal
%   error. However, if your images contain very few cells (in which case
%   popcontext analysis might not make sense), or if your cells are very
%   spread out, the automatic algorithm might fail and this allows you to
%   override the automatc detection. The overlap optimization is performed
%   integrated over all wells present in your dataset.
%
% Per well: Which image filename contains well and site information?
%   --> Determines which filename to parse to get image well and site
%   information from. See below for a list of known filename examples.
%
% Per well: What image-snake should be used?
%   --> The image snake represent the path with which your automatic
%   microscope acquires multiple imaged sites per well, and subsequently
%   the image-site number assigned to each location. Known snakes are
%   implemented for the ImageXpress Micro and Yokagawa CV7000, cellWoRx and
%   BD Pathway. When set to 'Automatic' the snake is inferred from the
%   naming style of your images. 
%
% Per well: What are the dimensions of your image grid? Enter 'Automatic'
% or for example '3,5' for 3 rows and 5 columns. 
%   --> Determines the grid size of imaged sites per well. Leave to
%   'Automatic' for a best guess depending on the number of sites imaged
%   per well, but this option allows the default to be overridden in case
%   automatic inference failed.
%
% Per well: Do you want to store PDF files of the overlap optimization and
% of the population context measurement of the first well?
%   --> Set to yes if you want automatic PDF creation of the image overlap
%   search and of the population context measurement of the first well. The
%   PDF will close the output figure during creation, but if multiple wells
%   are analyzed subsequent wells will be shown in a new figure. The PDF
%   files will be stored with a timestamp in the default output directory.
%
%
% Note: MeasurePopulationContext per well is not compatible with
% Batch/cluster analysis, as it requires all the measurements over the
% entire dataset (images from a single multiwell plate) to be present in
% the handles structure. But the code was originally written to be run on a
% handles.Measurements structure after CP analysis was complete, so it can
% realtively easily be reverted.
%
% On the interpretation of the LCD: LCD is a relative unit that depends on
% your image acquisition and microscope, as well as on the cell width
% parameter used to calculate it. If you want to compare absolute LCD
% values between different datasets make sure the images are of the same
% type and the cell diameter parameter has the same value.
%
% ******************** known file format examples *************************
% Known filename formats include the standard output formats of several
% automated high content microscopes including the API cellWoRx, MD
% ImageXpress, Yokagawa CV7000 and BD Pathway. The image snake (what image
% site per well corresponds to which position in the image grid) is
% automatically inferred from the naming style, but can be overwritten in
% the module settings.
%
% Here is a example list of understood filename formats:
%  strImageName = 'YD_50K_KY_P1_1_1_H07f00d0.png';
%  strImageName = 'SV40_FAKKO_GLY_GLYGGPP_GM1_noGM1_A01_01_w460_SegmentedCells.png';
%  strImageName = 'TDS_VSV_50K_P2_1_B02f0d0_SegmentedCells.png';
%  strImageName = 'Tf_KY_P2_20x_C10_19_w460.tif'
%  strImageName = 'Dapi - n000000.tif'
%  strImageName = 'DAPI_p53SLS - n000000.png';
%  strImageName = 'Frank_si_VSV_B02_25_w460.tif';
%  strImageName = '\Well H12\Dapi_p53SLS - n000000.tif';
%  strImageName = 'Olivia-1_A10_s1_w12E22EFEB-B167-43E0-A05F-997CCA19728A.tif'
%  strImageName = 'Frank-VSV-pH-dyn_B02_s1_w11BBF4034-97B9-4912-9DA5-6FBAF05BA7E4.tif'
%  strImageName = 'VACV_rescreen_CP078-1aa_K04_8_w530.png';
%  strImageName = 'HPV16_batch1_CP001-1ea_A20_9_w530.tif'
%  strImageName = 'BAC-siRNA-Optimisation_C01_s3CB0B5EFE-CA88-49D1-B8B8-2115D7B91A6F.png'
%  strImageName = 'RNABlockCourseSamplePlate01_B05_T0005F451L01A01Z01C02.png';
%
%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% MeasurePopulationContext was developed & written by Berend Snijder,
% University of Zurich, (c) 2011. 
%
% Website: http://www.cellprofiler.org & http://www.infectome.org
%
% $Revision: 8629 $

% %%%%%%%%%%%%%%%%%
% %%% VARIABLES %%%
% %%%%%%%%%%%%%%%%%

drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = Which object would you like to use for the cell population context measurement?
%choiceVAR01 = Nuclei
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
strObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What is the typical pixel diameter of your cells?
%defaultVAR02 = 140
%inputtypeVAR02 = popupmenu scale
intCellDiameter = str2double(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = How many times do you want to shrink to speed up the analysis? 'Automatic' or enter a value not smaller than 1.
%defaultVAR03 = Automatic
%inputtypeVAR03 = popupmenu scale
strShrinkFactor = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = Do you want to measure the population context per image or per well? Per well only works for "known" filenames (see help).
%choiceVAR04 = Well
%choiceVAR04 = Image
%inputtypeVAR04 = popupmenu
strImageOrWellMeasurement = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Per well: What is the amount of overlap in pixels between images of the same well? Enter 'Automatic' for automatic detection.
%defaultVAR05 = Automatic
%inputtypeVAR05 = popupmenu scale
strImageOverlap = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Per well: Which image filename contains well and site information?
%infotypeVAR06 = imagegroup
%inputtypeVAR06 = popupmenu
strImageName = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Per well: What image-snake should be used?
%choiceVAR07 = Automatic
%choiceVAR07 = ImageXpress Micro and Yokagawa CV7000
%choiceVAR07 = cellWoRx and BD Pathway
%inputtypeVAR07 = popupmenu
strImageSnake = char(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = Per well: What are the dimensions of your image grid? Enter 'Automatic' or for example '3,5' for 3 rows and 5 columns.
%choiceVAR08 = Automatic
%inputtypeVAR08 = popupmenu scale
strImageSiteGridDimensions = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = Per well: Do you want to store PDF files of the overlap optimization and of the population context measurement of the first well?
%choiceVAR09 = No
%choiceVAR09 = Yes
%inputtypeVAR09 = popupmenu
strStorePDF = char(handles.Settings.VariableValues{CurrentModuleNum,9});



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% MAKE MEASUREMENTS & SAVE TO HANDLES STRUCTURE %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

% store some current variable values
SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);

% temporarily overwrite the colormap for this figure?
strOldCMap = handles.Preferences.IntensityColorMap;
handles.Preferences.IntensityColorMap = 'jet';

% overwrite colormap and clear figure
if any(findobj == ThisModuleFigureNumber)    
    currentfig = CPfigure(handles,'Text',ThisModuleFigureNumber);
    colormap('jet')    
end

% add settings to the population context module output so that we can have
% multiple modules in the same pipeline
strSettingDescription = sprintf('Per%s_Cell%d',strImageOrWellMeasurement,intCellDiameter);


% %%%%%%%%%%%%%%%%%%%%%
% % STORE OBJECT AREA %
% independent of wether we work per well or per image, we also want to
% store the object area (so that people do not have to additionally run
% MeasureObjectAreaShape, which measures a lot more and is therefore
% slower). This measurement needs to be run per image as it requires the
% segmentation to be present. Retrieves the label matrix image that
% contains the segmented objects which will be measured with this module.
LabelMatrixImage =  CPretrieveimage(handles,['Segmented', strObjectName],ModuleName,'MustBeGray','DontCheckScale');
% get object count
NumObjects = max(LabelMatrixImage(:));
% default to empty matrix
matObjectAreas = [];
if  NumObjects > 0
    % get pixel size
    PixelSize = str2double(handles.Settings.PixelSize);
    % use regionprops to measure the object area
    warning('off', 'MATLAB:divideByZero');
    props = regionprops(LabelMatrixImage, 'Area');
    warning('on', 'MATLAB:divideByZero');
    matObjectAreas = cat(1,props.Area)*PixelSize^2;
end
% Save object area measurement
handles = CPaddmeasurements(handles, strObjectName, ['PopContext_',strSettingDescription,'_Area'], matObjectAreas);
% %%%%%%%%%%%%%%%%%%%%%



% if we do per-image analysis, draw a figure on each cellprofiler cycle
if strcmpi(strImageOrWellMeasurement,'image')

        % init output as empty cell array
        if (SetBeingAnalyzed == handles.Current.StartingImageSet)
            handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_LocalCellDensity']) = cell(1,handles.Current.NumberOfImageSets);
            handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_SingleCell']) = cell(1,handles.Current.NumberOfImageSets);
            handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_Edge']) = cell(1,handles.Current.NumberOfImageSets);
            handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_DistanceToEdge']) = cell(1,handles.Current.NumberOfImageSets);
            handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_PopulationSize']) = cell(1,handles.Current.NumberOfImageSets);
            handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_NearestNeighbor']) = cell(1,handles.Current.NumberOfImageSets);
        end
    
        % get object count and nuclei positions
        matObjectCount = handles.Measurements.Image.(['Count_',strObjectName]){SetBeingAnalyzed};
        matNucleiPositions = [handles.Measurements.(strObjectName).Location_Center_X{SetBeingAnalyzed},handles.Measurements.(strObjectName).Location_Center_Y{SetBeingAnalyzed}];
    
        % skip empty wells or images
        %if matObjectCount>0

            % shrink factor
            if strncmpi(strShrinkFactor,'auto',4)
                intShrinkFactor = 2;
            else
                intShrinkFactor = str2double(strShrinkFactor);
                if isnan(intShrinkFactor) || intShrinkFactor<1
                    CPmsgbox(sprintf('The shrink factor should be either ''Automatic'' or not smaller than 1; Now using ''Automatic'' settings.'),'Invalid shrink factor','help')
                    intShrinkFactor = 2;
                    % error('bs:InvalidShrinkFactor','The shrink factor should be either ''Automatic'' or not smaller than 1')
                end
            end
            % other params
            intFilterSigma = intCellDiameter / 6;
            intFilterSize = intCellDiameter;
            PSF = fspecial('gaussian',intFilterSize,intFilterSigma);

            % resize if necessary
            if intShrinkFactor~=1
                PSF = imresize(PSF,1/intShrinkFactor);
                PSF = PSF - min(PSF(:));
                PSF = PSF / max(PSF(:));
            end

            % do we parse wells? no :)
            boolWeCanParseWells = false;
            
            % get the image and corresponding image dimensions.
            matImage = CPretrieveimage(handles,strImageName,ModuleName);
            [intMaxWelPosY,intMaxWelPosX] = size(matImage);        

            % Do image analysis & draw figure
            [matEdgePerCell,matClosestDistanceToEdgePerCell,matLCDsForCurrentCells,matSingleForCurrentCells,matPopSizeForCurrentCells,matNearestNeighbourDistancePerCell] = doPopContextImageAnalysis(handles,matNucleiPositions,intShrinkFactor,intFilterSize,intFilterSigma,PSF,intMaxWelPosX,intMaxWelPosY,SetBeingAnalyzed,boolWeCanParseWells);
            

            % now store the calculated LCD, EDGE, Distance2Edge values, per
            % original image-index and object-index.
            handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_LocalCellDensity']){SetBeingAnalyzed} = matLCDsForCurrentCells;
            handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_SingleCell']){SetBeingAnalyzed} = matSingleForCurrentCells;
            handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_Edge']){SetBeingAnalyzed} = matEdgePerCell;
            handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_DistanceToEdge']){SetBeingAnalyzed} = matClosestDistanceToEdgePerCell;
            handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_PopulationSize']){SetBeingAnalyzed} = matPopSizeForCurrentCells;
            handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_NearestNeighbor']){SetBeingAnalyzed} = matNearestNeighbourDistancePerCell;
            

        %end
    
elseif strcmpi(strImageOrWellMeasurement,'well')
    % we're attempting per well analysis

    % on the first run, display a message on the output figure
    if (SetBeingAnalyzed == handles.Current.StartingImageSet)

        % draw/activate figure if user did not close the window
        if any(findobj == ThisModuleFigureNumber)

            % leave a message
            currentfig = CPfigure(handles,'Text',ThisModuleFigureNumber);
            TextString = 'The MeasurePopulationContext module is set to analyze wells from multiwell plates and will therefore start only at the last image cycle.';
            uicontrol(currentfig,'style','text','units','normalized','fontsize',handles.Preferences.FontSize,'HorizontalAlignment','center','string',TextString,'position',[.05 .45 .9 .1],'BackgroundColor',[.7 .7 .9],'tag','TextUIControl');

        end

    end

    % only analyze population context on the last image cycle
    if SetBeingAnalyzed == handles.Current.NumberOfImageSets


        % all user selected object count per image
        matObjectCount = cat(1,handles.Measurements.Image.(['Count_',strObjectName]){:});

        % initialize output over the full range
        handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_LocalCellDensity']) = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)';
        handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_Edge']) = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)';
        handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_DistanceToEdge']) = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)';
        handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_SingleCell']) = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)'; 
        handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_PopulationSize']) = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)'; 
        handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_NearestNeighbor']) = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)'; 
        
        % all user selected image file names
        cellFileNames = handles.Measurements.Image.(sprintf('FileName_%s',strImageName));

        % set filter sigma equal to the typical cell diameter / 6
        intFilterSigma = intCellDiameter / 6;
        intFilterSize = intCellDiameter;
        intImageBorderOverlap = intFilterSigma;

        % filter all file names to read out well-row, well-column and timepoint
        % information. In my opinion CP should have this covered :)
        boolWeCanParseWells = false;
        boolFigureWasClosed = CPwasfigureclosed(handles,ThisModuleFigureNumber);
        
        try
            % parse image filenames
            [cellRows, cellColumns, ~, cellTimepoints] = cellfun(@filterimagenamedata, cellFileNames,'UniformOutput',false);

            % all present wells & timepoints combinations
            matWellPositionsPerImage = [cell2mat(cellRows)',cell2mat(cellColumns)',cell2mat(cellTimepoints)'];

            % lookup image positions and attempt automatic microscope
            % type detection
            [cellImagePosition,cellstrMicroscopeType] = cellfun(@check_image_position,cellFileNames,'UniformOutput',false);
            matImagePosition = cell2mat(cellImagePosition);
            strMicroscopeType = unique(cellstrMicroscopeType);

            % check user setting for what type of image snake is
            % preferred (automatic, cwx/bd or md/cv7k).
            if strcmpi(strImageSnake(1),'c')% cellworx/bd
                strMicroscopeType = 'CW';
            elseif strcmpi(strImageSnake(1),'i')% md/cv7k
                strMicroscopeType = 'MD';
            end

            % get image snake
            if strncmpi(strImageSiteGridDimensions,'auto',4)
                % using automatic grid site detection
                matImageSnake = get_image_snake(max(matImagePosition),strMicroscopeType);
            elseif regexpi(strImageSiteGridDimensions,'\d{1,},\d{1,}')
                C = regexpi(strImageSiteGridDimensions,'(\d{1,}),(\d{1,})','tokens');
                matGridDimensions = [str2double(C{1}{1}),str2double(C{1}{2})];
                matImageSnake = get_image_snake(max(matImagePosition),strMicroscopeType,matGridDimensions);
            else
                CPmsgbox(sprintf('Your image site grid dimensions should be either ''Automatic'' or of the type ''number,number''. Switching to automatic best guess.'),'Incorrect grid dimension input','help')
                matImageSnake = get_image_snake(max(matImagePosition),strMicroscopeType);
            end
            
            % we should add these features as image measurements so people
            % can see what we inferred (This feature inference could
            % actually form a separate CP module)  
            handles.Measurements.Image.(['PopContext_',strImageName,'_WellRow']) = cellRows;
            handles.Measurements.Image.(['PopContext_',strImageName,'_WellColumn']) = cellColumns;
            handles.Measurements.Image.(['PopContext_',strImageName,'_ImageSite']) = cellImagePosition;
            handles.Measurements.Image.(['PopContext_',strImageName,'_Microscope']) = cellstrMicroscopeType;
            handles.Measurements.Image.(['PopContext_',strImageName,'_Timepoints']) = cellTimepoints;
            handles.Measurements.Image.(['PopContext_',strImageName,'_ImageSnake']) = matImageSnake; %note, we store a matrix in this feature..

            % set boolWeCanParseWells to true;
            boolWeCanParseWells = true;            
        catch  %#ok<CTCH>

            % catch problems and flag we can't parse filenames
            boolWeCanParseWells = false;
        end

        % the presence of NaNs indicates failed filename parsing as well.
        if ~boolWeCanParseWells || any(isnan(matWellPositionsPerImage(:))) || any(isnan(matImagePosition(:)))
            boolWeCanParseWells = false;
            CPmsgbox(sprintf('The MeasurePopulationContext module failed to successfully read your filenames. Switching to analyzing the population context per image.'),'Sorry, I do not understand your file names...','help')
        end

        % if we can't or shouldn't parse wells, set these variables to default
        % values that effectively result in per-image analysis.
        if ~boolWeCanParseWells
            % revert to default (no-well) image positions
            matWellPositionsPerImage = [(1:numel(cellFileNames))',ones(numel(cellFileNames),1),ones(numel(cellFileNames),1)];
            matImagePosition = ones(numel(cellFileNames),1);
            matImageSnake = zeros(2,100);
        end

        % get the image and corresponding image dimensions.
        matImage = CPretrieveimage(handles,strImageName,ModuleName);
        matImageSize = size(matImage);

        % generate explicit meta data, i.e. image-index and object-index
        cellNucleiMetaData = arrayfun(@(x,y) cat(2,repmat(y,x,1), (1:x)'), matObjectCount,(1:numel(matObjectCount))', 'UniformOutput',false);

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%% LOOP TO SEE OPTIMAL IMAGE OVERLAP... 
%         % We should probably also test negative image overlap values (i.e.
%         % image-spacing values). This might catch cases where nuclei at the
%         % edge of the image were not discarded.
%         %
%         % The automatic method is fast and robust, so it makes sense to do it
%         % in an automated fashion... but we might want to let the user
%         % overwrite the value, i.e. 'Automatic' is defualt, or enter number.


        % store unmodified nuclei positions in a cellArray
        cellOrigNucleiPositions = cellfun(@(x,y) [x,y],handles.Measurements.(strObjectName).Location_Center_X,handles.Measurements.(strObjectName).Location_Center_Y,'UniformOutput',false);

        % calculate new origins for each image, use these as offsets.
        matOrigNucleusOffsetX = matImageSnake(1,matImagePosition) * matImageSize(1,2);% width
        matOrigNucleusOffsetY = matImageSnake(2,matImagePosition) * matImageSize(1,1);% height

        % see if user wants automatic overlap detection
        if strncmpi(strImageOverlap,'auto',4) && boolWeCanParseWells

            iScoreCount = 0;
            matOverlapValues = -90:2:1200;
            matOverlapScore = [];
            boolOptimumFound = false;
            numOfHigherValues = 50;

            % loop over all possible overlap values, and break out the loop if
            % optimum is found.
            for iOverlap = matOverlapValues

                if ~boolOptimumFound
                    intImageBorderOverlap = iOverlap;
                elseif boolOptimumFound
                    % if the optimal overlap value is found, recalculate values using
                    % this setting, then break out of loop.
                    [foo,minIX]=nanmin(matOverlapScore); %#ok<ASGLU>
                    intImageBorderOverlap = matOverlapValues(minIX);
                end

                % calculate a rough 2D histogram of object locations
                % integrated
                % over all wells per plate, for varying image overlap values
                [~,~,~,matComplete2DBinCount] = calculateCombinedNucleiPositions(cellOrigNucleiPositions,matOrigNucleusOffsetX,matOrigNucleusOffsetY,matImageSnake,matImagePosition,matImageSize,intImageBorderOverlap);

                % Minimize the deviation of actual points from a smooth curve of
                % the X and Y projections. Smoothness is preferred / expected.
                matSumY = sum(matComplete2DBinCount,1);
                matSumX = sum(matComplete2DBinCount,2);
                matSmoothSumY = BSsmooth(matSumY(:)');
                matSmoothSumX = BSsmooth(matSumX(:)');

                if any(findobj == ThisModuleFigureNumber) || strcmpi(strStorePDF(1),'y')

                    % activate figure
                    currentfig = CPfigure(handles,'Text',ThisModuleFigureNumber);

                    % silly CP keeps drawing a menu label... try and
                    % prevent it like this.
                    FigUserData = get(currentfig,'UserData');
                    FigUserData.ImageFlag = 1;
                    set(currentfig,'UserData',FigUserData);                    
                    
                    % Remove uicontrols from last cycle
                    delete(findobj(ThisModuleFigureNumber,'tag','TextUIControl'));                    
                    if iOverlap==matOverlapValues(1)
                        clf(currentfig);
                    end

                    
%                     %%% A subplot of the figure window is set to display the original image.
                    hAx=subplot(2,2,1,'Parent',ThisModuleFigureNumber);
                    CPimagesc(matComplete2DBinCount,handles,hAx);
                    title(hAx,'2D cell density','FontSize',8);

%                     %%% A subplot of the figure window is set to display the colored label
%                     %%% matrix image.
                    hAx = subplot(2,2,2,'Parent',ThisModuleFigureNumber);
                    cla
                    hold on
                    plot(matSumY,'*b')
                    plot(matSmoothSumY,':r')
                    title(hAx,'Vertical sum projection and smooth fit','FontSize',8);                        
                    hold off

                    hAx = subplot(2,2,3,'Parent',ThisModuleFigureNumber);
                    cla
                    hold on
                    plot(matSumX,'*b')
                    plot(matSmoothSumX,':r')
                    title(hAx,'Horizontal sum projection and smooth fit','FontSize',8);
                    hold off


%                     %%% Error plot
                    hAx=subplot(2,2,4,'Parent',ThisModuleFigureNumber);
                    cla
                    hold on
                    plot(matOverlapValues(~isnan(matOverlapScore)),matOverlapScore(~isnan(matOverlapScore)))
                    axis tight
                    if boolOptimumFound
                       line([matOverlapValues(minIX),matOverlapValues(minIX)],[min(matOverlapScore),max(matOverlapScore)],'linestyle',':','color','r')
                       text(matOverlapValues(minIX),mean([min(matOverlapScore),max(matOverlapScore)]),sprintf('overlap = %d',matOverlapValues(minIX)),'fontsize',7)
                    end
                    title(hAx,'Error function','FontSize',8);                                                
                    hold off

                    drawnow

                end

                if boolOptimumFound
                    fprintf('%s: found optimal overlap value of %d in %d iterations\n',mfilename,intImageBorderOverlap,iScoreCount)
                    break
                end

                % add current score to score overview
                iScoreCount = iScoreCount + 1;

                matOverlapScore(iScoreCount) = ...
                    sum(abs(matSmoothSumY(:)-matSumY(:))) + ...
                    sum(abs(matSmoothSumX(:)-matSumX(:))); %#ok<AGROW>

                % minimum found if last X overlap-scores measurements are bigger
                % the -Xth (where X = numOfHigherValues).
                if length(matOverlapScore)>numOfHigherValues
                    if all(matOverlapScore(end-(numOfHigherValues-1):end)>matOverlapScore(end-numOfHigherValues))
                        boolOptimumFound = true;
                    end
                end

            end % end optimization loop

            % throw a message if we did not find an optimal overlap value, and
            % revert to a best guess.
            if ~boolOptimumFound
                CPmsgbox(sprintf('The MeasurePopulationContext module failed to find an optimum image overlap automatically. Setting image overlap to 1/6th the width of a typical cell (%d), but you can consider overwriting this value in the module settings.',intFilterSigma),'Population context warning','help')
                intImageBorderOverlap = intFilterSigma;
            end

        elseif ~strncmpi(strImageOverlap,'auto',4) && ~isnan(str2double(strImageOverlap))
            % user passed image overlap
            intImageBorderOverlap = str2double(strImageOverlap);
        else 
            % any other case, revert to default/automatic.
            intImageBorderOverlap = intFilterSigma;        
        end

        % (re)calculate nuclei positions with the final values for
        % intImageBorderOverlap
        [cellNucleiPositions,intMaxWelPosX,intMaxWelPosY,matComplete2DBinCount] = calculateCombinedNucleiPositions(cellOrigNucleiPositions,matOrigNucleusOffsetX,matOrigNucleusOffsetY,matImageSnake,matImagePosition,matImageSize,intImageBorderOverlap);
        
        % a shrink-factor to speed up the calculations, that does the final
        % calculation on an image that is at least 1000 pixels in one
        % dimension.
        if strncmpi(strShrinkFactor,'auto',4)
            intShrinkFactor = ceil(max(intMaxWelPosY,intMaxWelPosX)/1000);
        else
            intShrinkFactor = str2double(strShrinkFactor);
            if isnan(intShrinkFactor) || intShrinkFactor<1
                CPmsgbox(sprintf('The shrink factor should be either ''Automatic'' or not smaller than 1; Now using ''Automatic'' settings.'),'Invalid shrink factor','help')
                intShrinkFactor = ceil(max(intMaxWelPosY,intMaxWelPosX)/1000);                
            end
        end
        

%         %%% END OF IMAGE OVERLAP CALCULATIONS
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%% CALCULATE LCD, EDGE, DISTANCE-TO-EDGE, AND SINGLE-CELL PER WELL %%%
%         % loop over each well, & recalculate LCD & EDGE measurements

        % calculate the point spread function (PSF) for image dilution
        PSF = fspecial('gaussian',intFilterSize,intFilterSigma);


        % do resizing of PSF, after calculating the reference LCD values at full
        % resolution. I've tested and shrunk LCDs and corresponding cell counts per
        % pixel are approx. equal to their full resolution versions.
        %
        % Also, only calculate the PSF once, don't repeat per well/image.
        PSF = imresize(PSF,1/intShrinkFactor);
        PSF = PSF - min(PSF(:));
        PSF = PSF / max(PSF(:));

        % print figure if requested, than closre figure. if no print
        % requested, just clear figure if present 
        if any(findobj == ThisModuleFigureNumber)
            
            currentfig = CPfigure(handles,'text',ThisModuleFigureNumber);
            % Remove uicontrols from last cycle
            delete(findobj(ThisModuleFigureNumber,'tag','TextUIControl'));
            if strcmpi(strStorePDF(1),'y')
                strPrintName = sprintf('%s_%s_OverlapOptimization.pdf',datestr(now,30),ModuleName);
                h = CPmsgbox(sprintf('Please hold while printing figure to PDF...'),'Please hold...','help');
                try %#ok<TRYNC>
                    gcf2pdf(currentfig, handles.Current.DefaultOutputDirectory, strPrintName)
                catch objFoo
                    objFoo
                end
                close(h)                
                % close it, as gcf2pdf places it out of visible area of
                % screen...
                close(currentfig)

                % and create a new output figure for showing (at least) the
                % first well popcontext measurements
                CPfigure(handles,'Text',ThisModuleFigureNumber);
                % Remove uicontrols from last cycle
                delete(findobj(ThisModuleFigureNumber,'tag','TextUIControl'));                
            end
            clf(currentfig);
            
        end        
        
        % loop over each individual well and timepoint (if found according to
        % our parsing...)
        iWellCounter = 0;
        for iPos = unique(matWellPositionsPerImage,'rows')'

            fprintf('%s: processing well: row %d, col %d (t=%d)\n',mfilename,iPos(1),iPos(2),iPos(3))

            % look up images corresponding to current well
            matImageIX = ismember(matWellPositionsPerImage,iPos','rows');

            % get all nuclei positions and meta data for current well
            matNucleiPositions = cat(1,cellNucleiPositions{matImageIX});
            matNucleiMetaData = cat(1,cellNucleiMetaData{matImageIX});

            % skip empty wells or images
            if size(matNucleiPositions,1)==0
                fprintf('%s: skipping well: no cells\n',mfilename)
                continue
            end

            % count how many wells were processed
            iWellCounter = iWellCounter + 1;   
            
            % Do image analysis & plot figure (force figure plotting if PDf
            % is requested, otherwise, leave it to if the figure was closed
            % or not)  
            if strcmpi(strStorePDF(1),'y') && iWellCounter==1
                [matEdgePerCell,matClosestDistanceToEdgePerCell,matLCDsForCurrentCells,matSingleForCurrentCells,matPopSizeForCurrentCells,matNearestNeighbourDistancePerCell] = doPopContextImageAnalysis(handles,matNucleiPositions,intShrinkFactor,intFilterSize,intFilterSigma,PSF,intMaxWelPosX,intMaxWelPosY,iPos,boolWeCanParseWells,false);
            else
                [matEdgePerCell,matClosestDistanceToEdgePerCell,matLCDsForCurrentCells,matSingleForCurrentCells,matPopSizeForCurrentCells,matNearestNeighbourDistancePerCell] = doPopContextImageAnalysis(handles,matNucleiPositions,intShrinkFactor,intFilterSize,intFilterSigma,PSF,intMaxWelPosX,intMaxWelPosY,iPos,boolWeCanParseWells,boolFigureWasClosed);
            end
            
            % if requested to print the first figure as pdf plot, do so.
            if strcmpi(strStorePDF(1),'y') && iWellCounter==1
                currentfig = CPfigure(handles,'Text',ThisModuleFigureNumber);
                % Remove uicontrols from last cycle
                delete(findobj(ThisModuleFigureNumber,'tag','TextUIControl'));                
                strPrintName = sprintf('%s_%s_FirstWellExample.pdf',datestr(now,30),ModuleName);
                h = CPmsgbox(sprintf('Please hold while printing figure to PDF...'),'Please hold...','help');
                try %#ok<TRYNC>
                    gcf2pdf(currentfig, handles.Current.DefaultOutputDirectory, strPrintName)
                catch objFoo
                    objFoo
                end                
                close(h)
                % close it, as gcf2pdf places it out of visible area of
                % screen...
                close(currentfig)
                
                % and create a new output figure, if the user did not close
                % the figure before...
                if ~boolFigureWasClosed
                    CPfigure(handles,'Text',ThisModuleFigureNumber);
                    % Remove uicontrols from last cycle
                    delete(findobj(ThisModuleFigureNumber,'tag','TextUIControl'));                    
                end
            end            
            

            % now store the calculated LCD, EDGE, Distance2Edge values, per
            % original image-index and object-index.
            for iImage = unique(matNucleiMetaData(:,1))'
                matImageIX = ismember(matNucleiMetaData(:,1),iImage);
                matObjectIX = matNucleiMetaData(matImageIX,2);

                handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_LocalCellDensity']){iImage}(matObjectIX) = matLCDsForCurrentCells(matImageIX);
                handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_SingleCell']){iImage}(matObjectIX) = matSingleForCurrentCells(matImageIX);
                handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_Edge']){iImage}(matObjectIX) = matEdgePerCell(matImageIX);
                handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_DistanceToEdge']){iImage}(matObjectIX) = matClosestDistanceToEdgePerCell(matImageIX);
                handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_PopulationSize']){iImage}(matObjectIX) = matPopSizeForCurrentCells(matImageIX);
                handles.Measurements.(strObjectName).(['PopContext_',strSettingDescription,'_NearestNeighbor']){iImage}(matObjectIX) = matNearestNeighbourDistancePerCell(matImageIX);
            end

        end

    end% check for last image cycle
    
end% check if we do per well or per image analysis

% restore default colormap (we overwrite it in the beginning)
handles.Preferences.IntensityColorMap = strOldCMap;    

end% function



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%                         HERE BE SUBFUNCTIONS                        %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function boolFigureWasClosed = CPwasfigureclosed(handles,ThisModuleFigureNumber)
    boolFigureWasClosed = false;%default
    matAllModuleFigureHandles = NaN(1,handles.Current.NumberOfModules);
    for i = 1:handles.Current.NumberOfModules
        matAllModuleFigureHandles(i) = handles.Current.(sprintf('FigureNumberForModule%02d',i));
    end
    % this seems to happen when no figures were opened at the start of CP.
    if numel(unique(matAllModuleFigureHandles))==1
        boolFigureWasClosed = true;
    end
    % if figures were opened at the beginning of CP, check if this figure
    % is present now...
    if ~boolFigureWasClosed
        boolFigureWasClosed = ~any(findobj == ThisModuleFigureNumber);
    end
end


function [matEdgePerCell,matClosestDistanceToEdgePerCell,matLCDsForCurrentCells,matSingleForCurrentCells,matPopSizeForCurrentCells,matNearestNeighbourDistancePerCell] = doPopContextImageAnalysis(handles,matNucleiPositions,intShrinkFactor,intFilterSize,intFilterSigma,PSF,intMaxWelPosX,intMaxWelPosY,iPos,boolWeCanParseWells,boolFigureWasClosed)

    if ~exist('boolFigureWasClosed','var')
        % if not passed, assume it was open. drawing will check for its
        % presence anyway.
        boolFigureWasClosed=false;
    end
    
    % measure the neirest neighbour distance per cell
    matNearestNeighbourDistancePerCell = nan(size(matNucleiPositions,1),1);
    if size(matNucleiPositions,1)>1
        matNotSelfIX = ~eye(size(matNucleiPositions,1),size(matNucleiPositions,1));
        for iCell = 1:size(matNucleiPositions,1)
            matNearestNeighbourDistancePerCell(iCell,1) = min( sqrt( ...
                    (matNucleiPositions(matNotSelfIX(:,iCell),1) - matNucleiPositions(iCell,1)) .^2 + ...
                    (matNucleiPositions(matNotSelfIX(:,iCell),2) - matNucleiPositions(iCell,2)) .^2 ...
                ) ); 
        end    
    end
    
    
    % round off nuclei positions (shrunk if necessary) to be able to index
    % a matrix with them 
    matNucleiPositions = ceil(matNucleiPositions / intShrinkFactor);
    
    % what is the LCD margin to be considered single cells (i.e. if LCD < minLCD * factor) you are single. 
    intSingleCellFactor = 1.1;
    
    % create map with dots for each cell
    % perhaps work with dots, gaussian blurred...
    matImageMapWithDots = zeros(ceil(intMaxWelPosY / intShrinkFactor),ceil(intMaxWelPosX / intShrinkFactor));
    matImageMapWithDots(sub2ind(size(matImageMapWithDots),matNucleiPositions(:,2),matNucleiPositions(:,1))) = 1;

    % gaussian blur mask. note: mask size determines radius for LCD & EDGE
    % measurement, which is crucial
    matImageMap = imfilter(matImageMapWithDots,PSF,'symmetric','conv');
    % this gives LCD values that are independent of the shrink factor
    % (tested between 2 - 7)
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDGE DETECTION %%%%%%%%%%%
%     %%% Calculate edges: make binary map, expand dots with sigma from mask,
%     %%% look for empty areas (excluding too small empty areas), mark those
%     %%% cells bordering empty areas as edge cells.

    % look for empty space
    
    % strategy 1: expand dots, look for edges.
    se = strel('disk',ceil(intFilterSigma/intShrinkFactor));
    matImageMapWithEmptySpace = imdilate(matImageMapWithDots,se);
    
    %matImageMapWithEmptySpace = bwmorph(matImageMapWithDots,'thicken',ceil(intFilterSigma/intShrinkFactor));
    %matImageMapWithEmptySpace = bwmorph(matImageMapWithEmptySpace,'close');
    %matImageMapWithEmptySpace = bwmorph(matImageMapWithEmptySpace,'fill');
    matImageMapWithEmptySpace(:,1) = 1;
    matImageMapWithEmptySpace(:,end) = 1;
    matImageMapWithEmptySpace(1,:) = 1;
    matImageMapWithEmptySpace(end,:) = 1;

    matImageMapWithEmptySpace = ~matImageMapWithEmptySpace;
    % alternative: threshold gaussian blur...
%     matImageMapWithEmptySpace = matImageMap < 1; % empty = 1/true


%     %%% remove too small empty space regions
    % first, label all objects
    matEmptySpaceLabelMatrix = bwlabel(matImageMapWithEmptySpace);
    % get area counts per object
    props = regionprops(matEmptySpaceLabelMatrix,'Area'); %#ok<MRPBW>
    matEmptySpaceRegionSize = cat(1,props.Area);

    % find too small empty areas, and exclude these (size somehow in
    % relation to mask-size used in image blur)
    matTooSmallEmptyAreasObjectID = find(matEmptySpaceRegionSize < 2*(pi*ceil(intFilterSigma/intShrinkFactor)^2) );
    % remove too small empty-space-objects
    if ~isempty(matTooSmallEmptyAreasObjectID)
        %fprintf('%s: \tdiscarding %d too small empty areas\n',mfilename,length(matTooSmallEmptyAreasObjectID))
        matImageMapWithEmptySpace(ismember(matEmptySpaceLabelMatrix,matTooSmallEmptyAreasObjectID)) = 0;
    end
%     %%%

    % behaviour of edge() changed over different matlab versions... :(
%     matImageMapWithEmptySpace = edge(uint8(matImageMapWithEmptySpace),'sobel');
    matImageMapWithEmptySpace = detectEdge(uint8(matImageMapWithEmptySpace));

    % find edges of scratch-mask
    [matEmptySpaceEdgeX,matEmptySpaceEdgeY] = find(matImageMapWithEmptySpace);

    % create euclidean distance to closest sratch for each cell position
    % (note that for a distance ranking, the sqrt() can be ommitted)
    if size(matNucleiPositions,1)>0
        matIDOfClosestCellPerEdgePixel = NaN(size(matEmptySpaceEdgeX,1),1);        
        for iPixel = 1:size(matEmptySpaceEdgeX,1)
            [foo, matIDOfClosestCellPerEdgePixel(iPixel)] = min( sqrt( ...
                    (matNucleiPositions(:,1) - matEmptySpaceEdgeY(iPixel)) .^2 + ...
                    (matNucleiPositions(:,2) - matEmptySpaceEdgeX(iPixel)) .^2 ...
                ) ); %#ok<ASGLU>
        end
        matEdgePerCell = zeros(size(matNucleiPositions,1),1);
        matEdgePerCell(unique(matIDOfClosestCellPerEdgePixel)) = 1;
    else
        matEdgePerCell = zeros(size(matNucleiPositions,1),1);
    end


    % calculate the closest distance between each cell and an edge pixel
    matClosestDistanceToEdgePerCell = NaN(size(matNucleiPositions,1),1);
    if ~isempty(matEmptySpaceEdgeY) && size(matNucleiPositions,1) > 0
        for iCell = 1:size(matNucleiPositions,1)
            matClosestDistanceToEdgePerCell(iCell) = min( sqrt( ...
                    (matEmptySpaceEdgeY - matNucleiPositions(iCell,1)) .^2 + ...
                    (matEmptySpaceEdgeX - matNucleiPositions(iCell,2)) .^2 ...
                ) );
        end    
    elseif size(matNucleiPositions,1) == 0
        matClosestDistanceToEdgePerCell = [];
    else
        % there are no edges in current image, set to max distance
        % possible. (i.e. diagonal on entire image)
        matClosestDistanceToEdgePerCell(:,1) = round(sqrt(intMaxWelPosX^2 + intMaxWelPosY^2));
    end

    % calculate LCD per cell
    matLCDsForCurrentCells = arrayfun(@(x,y) matImageMap(x,y),matNucleiPositions(:,2),matNucleiPositions(:,1));    

    % Calculate if cells are all alone, if so they are 'single' cells
    intMaxPSFValue = max(PSF(:));
    matSingleForCurrentCells = matLCDsForCurrentCells<=(intMaxPSFValue*intSingleCellFactor);

    % population size = total number of cells (per image or per well)
    matPopSizeForCurrentCells = repmat(size(matNucleiPositions,1),size(matSingleForCurrentCells));
    
    % get figure number for current module number
    ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',handles.Current.CurrentModuleNumber]);
    
    
    % Draw figure per well or per image
    if any(findobj == ThisModuleFigureNumber) && ~boolFigureWasClosed
        % Remove uicontrols from last cycle
        delete(findobj(ThisModuleFigureNumber,'tag','TextUIControl'));        

        CPfigure(handles,'Text',ThisModuleFigureNumber);

        % draw stuff
        hAx = subplot(2,3,1,'Parent',ThisModuleFigureNumber);
%         cla
        CPimagesc(matImageMap,handles,hAx);
        if boolWeCanParseWells
            title(hAx,sprintf('Well row=%d col=%d (t=%d).\nShrunk %dx. Cell diameter %d',iPos(1),iPos(2),iPos(3),intShrinkFactor,intFilterSize),'fontsize',8);
        else
            title(hAx,sprintf('Image cycle %d.\nShrunk %dx. Cell diameter %d',iPos(1),intShrinkFactor,intFilterSize),'fontsize',8);
        end

        hAx = subplot(2,3,2,'Parent',ThisModuleFigureNumber);
        cla
        hold on
        scatter(matNucleiPositions(:,1),-matNucleiPositions(:,2),15,matEdgePerCell,'Marker','*')
        ylim([-ceil(intMaxWelPosY/intShrinkFactor),0]);xlim([0,ceil(intMaxWelPosX/intShrinkFactor)])
        title(hAx,sprintf('Edge Cells (%d%%)',round(mean(matEdgePerCell)*100)),'FontSize',8)
        hold off

        hAx = subplot(2,3,3,'Parent',ThisModuleFigureNumber);
        cla
        hold on        
        scatter(matNucleiPositions(:,1),-matNucleiPositions(:,2),15,matClosestDistanceToEdgePerCell,'Marker','*')
        ylim([-ceil(intMaxWelPosY/intShrinkFactor),0]);xlim([0,ceil(intMaxWelPosX/intShrinkFactor)])
        title(hAx,'Minimal distance to edge','FontSize',8)
        hold off        

        hAx = subplot(2,3,4,'Parent',ThisModuleFigureNumber);
        cla
        hold on        
        scatter(matNucleiPositions(:,1),-matNucleiPositions(:,2),15,matSingleForCurrentCells,'Marker','*')
        ylim([-ceil(intMaxWelPosY/intShrinkFactor),0]);xlim([0,ceil(intMaxWelPosX/intShrinkFactor)])
        title(hAx,sprintf('Single Cells (%d%%)',round(mean(matSingleForCurrentCells)*100)),'FontSize',8)
        hold off        

        hAx = subplot(2,3,5,'Parent',ThisModuleFigureNumber);
        cla
        hold on        
        hist(matLCDsForCurrentCells)
        title(hAx,'LCD histogram','FontSize',8)    
        hold off        

        hAx = subplot(2,3,6,'Parent',ThisModuleFigureNumber);
        cla
        hold on        
        scatter(matNucleiPositions(:,1),-matNucleiPositions(:,2),15,matNearestNeighbourDistancePerCell,'Marker','*')
        ylim([-ceil(intMaxWelPosY/intShrinkFactor),0]);xlim([0,ceil(intMaxWelPosX/intShrinkFactor)])
        title(hAx,sprintf('Nearest Neighbour Distance\n(Suggested diameter = %d)',round(1.1*nanmedian(matNearestNeighbourDistancePerCell*6))),'FontSize',8)
        hold off

        drawnow

    end    
    
    
end


function [cellNucleiPositions,intMaxWelPosX,intMaxWelPosY,matComplete2DBinCount] = calculateCombinedNucleiPositions(cellOrigNucleiPositions,matOrigNucleusOffsetX,matOrigNucleusOffsetY,matImageSnake,matImagePosition,matImageSize,intImageBorderOverlap)


    % calculate nucleus offsets for varying image border overlap values
    matNucleusOffsetX = matOrigNucleusOffsetX - (matImageSnake(1,matImagePosition) * intImageBorderOverlap);
    matNucleusOffsetY = matOrigNucleusOffsetY - (matImageSnake(2,matImagePosition) * intImageBorderOverlap);

    % get max well dimensions, for 2D binning later
    intMaxWelPosX = (max(matImageSnake(1,:))+1) * matImageSize(1,2) - max(matImageSnake(1,:) * intImageBorderOverlap);% max well width
    intMaxWelPosY = (max(matImageSnake(2,:))+1) * matImageSize(1,1) - max(matImageSnake(2,:) * intImageBorderOverlap);% max well height

    % get original nuclei positions
    cellNucleiPositions = cellOrigNucleiPositions;

    % add offsets the original nuclei positions accounting for image
    % overlaps
    for i = 1:length(cellNucleiPositions)
        if ~isempty(cellNucleiPositions{i})
            cellNucleiPositions{i}(:,1) = cellNucleiPositions{i}(:,1) + matNucleusOffsetX(i);
            cellNucleiPositions{i}(:,2) = cellNucleiPositions{i}(:,2) + matNucleusOffsetY(i);
        end
    end

    if nargout > 3
        % plot total LCD map of all wells combined 
        matAllWellPositions = cat(1,cellNucleiPositions{:});
        % get 2D heatmap of cell count
        matDescriptor = [0, intMaxWelPosY, round(intMaxWelPosY / 50);
                         0, intMaxWelPosX, round(intMaxWelPosX / 50)];
        [matComplete2DBinCount]=histogram2(matAllWellPositions(:,2)',matAllWellPositions(:,1)',matDescriptor);
    end




end % function







function [intRow, intColumn, strWellName, intTimepoint] = filterimagenamedata(strImageName)

    intRow = NaN;
    intColumn = NaN;
    strWellName = NaN;
    intTimepoint  = 1;

    if nargin == 0
%         strImageName = '070420_Tf_KY_P2_20x_E01_19_w460.tif'
%         strImageName = '070610_Tfn_MZ_kinasescreen_CP022-1cd_A01f0d0.tif'
%         strImageName = 'Y:\Data\Users\Jean-Philippe\p53hitvalidation\ko_p53Pro72_plate1_triplicate1\TIFF\Well A24'
%         strImageName = '080611-olivia-1_A10_s1_w12E22EFEB-B167-43E0-A05F-997CCA19728A.tif'
%         strImageName = '080815-Frank-VSV-pH-dyn_F02_s1_w1E073FDA7-1105-4532-B633-E9D89C3F23FB.tif'
%         strImageName = '\\nas-biol-imsb-1\share-3-$\Data\Users\HPV16_DG\2008-08-16_HPV16_batch1_CP003-1ec\Well A001\Dapi - n000000.tif'
%         strImageName = 'BAC-siRNA-Optimisation_C01_s3CB0B5EFE-CA88-49D1-B8B8-2115D7B91A6F.png'
%         strImageName = '110420GpiGfpHoechstFakInhib_t_0021_C03_s6_w231366D01-7DCF-4473-A4A4-7A78094ADD3E.png';
%         strImageName = '080815-Frank-VSV-pH-dyn_B02_s1_w11BBF4034-97B9-4912-9DA5-6FBAF05BA7E4.tif'
%         strImageName = '110519_NB_E37_AAVRS1DOX_VS1_H12_12_w460.png';
%       strImageName = 'SettingB_E05_w167678E42-9203-4F07-9F36-EE22FDBE1B90.png';
%          strImageName = '111106-InSitu-MovieTest01_G09_T0003F004L01A01Z01C01.png';
        strImageName = 'bDZ01-1A_wD17_s3_z0_t0_cGFP.tif';   
    end

    strWellName = char(strrep(regexp(strImageName,'_[A-Z]\d\d_','Match'),'_',''));
    strWellName2 = char(strrep(regexp(strImageName,'_[A-Z]\d\df','Match'),'_',''));
    strWellName3 = char(strrep(regexp(strImageName,'Well [A-Z]\d{2,}','Match'),'_',''));
    strWellName4 = regexp(strImageName,'_\w\d\d_s\d{1,}_w\d','Match');
    strWellName5 = regexp(strImageName,'_w\w\d{2,}_','Match');
    


    % look for timepoint information
    cellstrTimepoint = regexp(strImageName,'_[tT]_?(\d{1,})_?','Tokens');
    if not(isempty(cellstrTimepoint))
        intTimepoint = str2double(cellstrTimepoint{1});
    else
        intTimepoint = 1;
    end

    % see what matches and get info from there
    if not(isempty(strWellName))
        intRow=double(strWellName(end,1))-64;
        intColumn=str2double(strWellName(end,2:3));
        strWellName=strWellName(end,1:end);
    elseif not(isempty(strWellName2))
        intRow=double(strWellName2(end,1))-64;
        intColumn=str2double(strWellName2(end,2:3));
        strWellName=strWellName2(end,1:end-1);
    elseif not(isempty(strWellName3))
        strImageData = regexp(strImageName,'Well ([A-Z])(\d{2,})','Tokens');
        intRow=double(strImageData{1}{1})-64;
        intColumn=str2double(strImageData{1}{2});
        strWellName=sprintf('%s%.02d',strImageData{1}{1},intColumn);
    elseif not(isempty(strWellName4))
%         %%% MD
        strImageData = regexpi(strImageName,'_(\w)(\d\d)_s','Tokens');
        intRow=double(strImageData{1})-64;
        intColumn=str2double(strImageData{2});
        strWellName=[strImageData{1},strImageData{2}];
    elseif not(isempty(strWellName5))
% %         %%% iBRAIN2
        strImageData = regexpi(strImageName,'_w(\w)(\d{2,})_','Tokens');
        intRow=double(strImageData{1}{1})-64;
        intColumn=str2double(strImageData{1}{2});
        strWellName=[strImageData{1}{1},strImageData{1}{2}];
    else
        intRow = NaN;
        intColumn = NaN;
        strWellName = NaN;
        warning('bs:UnknownFileName','filterimagenamedata: unable to get well data from image name %s',strImageName)
    end


end



function [matImagePosition,strMicroscopeType] = check_image_position(strImageName)
% help for check_image_position()
% BS, 082015
% usage [matImagePosition,strMicroscopeType] = check_image_position(strImageName)
%
% possible values for strMicroscopeType are 'CW', 'BD', 'MD' and 'CV7K'.

    if nargin == 0
%         strImageName = '061224_YD_50K_KY_P1_1_1_H07f00d0.png';
%         strImageName = '090313_SV40_FAKKO_GLY_GLYGGPP_GM1_noGM1_A01_01_w460_SegmentedCells.png';
%         strImageName = '070314_TDS_VSV_50K_P2_1_B02f0d0_SegmentedCells.png';
%         strImageName = '070420_Tf_KY_P2_20x_C10_19_w460.tif'
%         strImageName = 'Dapi - n000000.tif'
%         strImageName = 'DAPI_p53SLS - n000000.png';
%         strImageName = '040423_frank_si_VSV_B02_25_w460.tif';
%         strImageName = 'Y:\Data\Users\Jean-Philippe\p53hitvalidation\ko_p53Pro72_plate1_triplicate1\TIFF\Well H12\Dapi_p53SLS - n000000.tif';
%         strImageName = '080611-olivia-1_A10_s1_w12E22EFEB-B167-43E0-A05F-997CCA19728A.tif'
%         strImageName = '080815-Frank-VSV-pH-dyn_B02_s1_w11BBF4034-97B9-4912-9DA5-6FBAF05BA7E4.tif'
%         strImageName = '081008_VV_rescreen_CP078-1aa_K04_8_w530.png';
%         strImageName = '2008-08-14_HPV16_batch1_CP001-1ea_A20_9_w530.tif'
%         strImageName = 'BAC-siRNA-Optimisation_C01_s3CB0B5EFE-CA88-49D1-B8B8-2115D7B91A6F.png'
%          strImageName = '110930-RNABlockCourseSamplePlate01_B05_T0005F451L01A01Z01C02.png';
        strImageName = 'bDZ01-1A_wD17_s3_z0_t0_cGFP.tif'; 
    end

    strMicroscopeType = '';
    matImagePosition = NaN;

    % CW
    strNomenclature1 = regexp(strImageName,'f\d\dd\d\.','Match');
    strNomenclature1a = regexp(strImageName,'f\dd\d\.','Match');
    strNomenclature1b = regexp(strImageName,'f\d\dd\d','Match');
    strNomenclature1c = regexp(strImageName,'f\dd\d','Match');
    strNomenclature2 = regexp(strImageName,'_w\d\d\d\.','Match');
    strNomenclature3 = regexp(strImageName,' - n\d{2,}\.','Match');

    % MD MICROEXPRESS
    strNomenclature4 = regexp(strImageName,'_\w\d\d_s\d{1,}_w\d','Match');
    % MD with only one channel matches this but not the previous
    strNomenclature4a = regexp(strImageName,'_\w\d\d_s\d{1,}[A-Z0-9\-]{36}','Match');
    strNomenclature4b = regexp(strImageName,'_\w\d\d_\d{1,}_w\d','Match');

    % CV7K
    strNomenclature5 = regexp(strImageName, '_([^_]{3})_(T\d{4})(F\d{3})(L\d{2})(A\d{2})(Z\d{2})(C\d{2})', 'Match');
    
    % iBRAIN2
    strNomenclature6 = regexp(strImageName, '_s\d{1,}_', 'Match');
    


    if not(isempty(strNomenclature1))
        strMicroscopeType = 'CW';
        strImagePosition = regexp(strImageName,'f\d\dd\d\.','Match');
        strImagePosition = strImagePosition{1,1}(1,1:3);
        strImagePosition = strrep(strImagePosition,'f','');
        matImagePosition = str2double(strImagePosition)+1;
    elseif not(isempty(strNomenclature2))
        strMicroscopeType = 'CW';
        strImagePosition = regexp(strImageName,'_\d\d_w\d\d\d\.','Match');
        if not(isempty(strImagePosition))
            matImagePosition = str2double(strImagePosition{1,1}(1,2:3));
        else
            strImagePosition = regexp(strImageName,'_\d_w\d\d\d\.','Match');
            matImagePosition = str2double(strImagePosition{1,1}(1,2));
        end
    elseif not(isempty(strNomenclature1a))
        strMicroscopeType = 'CW';
        strImagePosition = regexp(strImageName,'f\dd\d\.','Match');
        strImagePosition = strImagePosition{1,1}(1,1:2);
        strImagePosition = strrep(strImagePosition,'f','');
        matImagePosition = str2double(strImagePosition)+1;
    elseif not(isempty(strNomenclature1b))
        strMicroscopeType = 'CW';
        strImagePosition = regexp(strImageName,'f\d\dd\d','Match');
        strImagePosition = strImagePosition{1,1}(1,1:3);
        strImagePosition = strrep(strImagePosition,'f','');
        matImagePosition = str2double(strImagePosition)+1;
    elseif not(isempty(strNomenclature1c))
        strMicroscopeType = 'CW';
        strImagePosition = regexp(strImageName,'f\dd\d','Match');
        strImagePosition = strImagePosition{1,1}(1,2);
        matImagePosition = str2double(strImagePosition)+1;
    elseif not(isempty(strNomenclature3))
        strMicroscopeType = 'BD';
        strImagePosition = regexpi(strImageName,' - n(\d{2,})\.','Tokens');
        if ~isempty(strImagePosition)
            matImagePosition = str2double(strImagePosition{1}{1})+1;
        end
    elseif not(isempty(strNomenclature4))
        strMicroscopeType = 'MD';
        strImagePosition = regexpi(strImageName,'_\w\d\d_s(\d{1,})_w\d','Tokens');
        matImagePosition = str2double(strImagePosition{1});
    elseif not(isempty(strNomenclature4a))
        strMicroscopeType = 'MD';
        strImagePosition = regexpi(strImageName,'_\w\d\d_s(\d{1,})[A-Z0-9\-]{36}','Tokens');
        matImagePosition = str2double(strImagePosition{1});
    elseif not(isempty(strNomenclature4b))
        strMicroscopeType = 'MD';
        strImagePosition = regexpi(strImageName,'_\w\d\d_(\d{1,})_w\d','Tokens');
        matImagePosition = str2double(strImagePosition{1});
    elseif not(isempty(strNomenclature5))
        strMicroscopeType = 'CV7K';
        strImagePosition = regexp(strImageName, '_([^_]{3})_(T\d{4})F(\d{3})(L\d{2})(A\d{2})(Z\d{2})C(\d{2})', 'Tokens');
        matImagePosition = str2double(strImagePosition{1}(3));
    elseif not(isempty(strNomenclature6))
        strMicroscopeType = 'MD'; %iBRAIN2 = MD
        strImagePosition = regexp(strImageName, '_s(\d{1,})_', 'Tokens');
        matImagePosition = str2double(strImagePosition{1}(1));
    else
        warning('bs:UnknownFileName','unknown file name %s',strImageName)
    end
end




function [matImageSnake,matStitchDimensions] = get_image_snake(intMaxImagePosition, strMicroscopeType, matGridDimensions)
% help for get_image_snake()
% BS, 082015
% usage [intImagePosition,strMicroscopeType] = get_image_snake(strImageName)
%
% possible values for strMicroscopeType are 'CW', 'BD', 'MD' and 'CV7K'.

    if nargin == 0
        intMaxImagePosition = 49;
        strMicroscopeType = 'CW';
    end

    if nargin > 0 & ~isnumeric(intMaxImagePosition)
        error('%s: intMaxImagePosition should be a number.',mfilename)
    end

    if nargin == 2
        if ~strcmpi(strMicroscopeType,'CW') && ...
           ~strcmpi(strMicroscopeType,'BD') && ...
           ~strcmpi(strMicroscopeType,'MD') && ...
           ~strcmpi(strMicroscopeType,'CV7K') && ...
            error('%s: unrecognized microscope type ''%s''. Allowed values are ''CW'', ''BD'', ''MD'' and ''CV7K''.',mfilename,strMicroscopeType)
        end
    end

    % Default to a lot of zeros
    matImageSnake = zeros(2,100);
    matStitchDimensions = [1,1];

    % look up the most likely configuration of number of images...
    if intMaxImagePosition == 0
        matImageSnake = [0;0];
        matStitchDimensions = [1,1];
    elseif intMaxImagePosition ==1
        matImageSnake = [0;0];
        matStitchDimensions = [1,1];
    elseif intMaxImagePosition ==2
        matStitchDimensions = [1,2];
    elseif intMaxImagePosition ==4
        matStitchDimensions = [2,2];
    elseif intMaxImagePosition ==5
        matStitchDimensions = [2,3];
    elseif intMaxImagePosition == 9
        matStitchDimensions = [3,3];
    elseif intMaxImagePosition == 12
        matStitchDimensions = [4,3];
    elseif intMaxImagePosition == 15
        matStitchDimensions = [5,3];
    elseif intMaxImagePosition == 16
        matStitchDimensions = [4,4];
    elseif intMaxImagePosition == 20
        matStitchDimensions = [5,4];
    elseif intMaxImagePosition == 25
        matStitchDimensions = [5,5];
    elseif intMaxImagePosition == 36
        matStitchDimensions = [6,6];
    elseif intMaxImagePosition == 192
        matStitchDimensions = [16,12];
    elseif intMaxImagePosition == 180
        matStitchDimensions = [15,12];
    elseif intMaxImagePosition == 42
        matStitchDimensions = [7,6];
    elseif intMaxImagePosition == 30
        matStitchDimensions = [6,5];
    elseif intMaxImagePosition == 48
        matStitchDimensions = [8,6];
    elseif intMaxImagePosition == 49
        matStitchDimensions = [7,7];
    elseif intMaxImagePosition == 63
        matStitchDimensions = [9,7];
    elseif intMaxImagePosition == 80
        matStitchDimensions = [8,10];
    else
        % simple algorithm to guess the correct dimensions for image
        % stitching
        matI = round(sqrt(intMaxImagePosition))-5 : round(sqrt(intMaxImagePosition))+5; % set up search dimension
        matII = matI' * matI; % symmetric search space
        [i,j] = find(triu(matII)==intMaxImagePosition); % find the multiple that matches our image count
        matStitchDimensions = sort([matI(i(1)),matI(j(1))],'descend');% assume more rows than columns
    end

    % if user defined stitch dimensions, use those...
    if nargin==3
        matStitchDimensions = matGridDimensions;
    end
    
    % generate the so-called image snake depending on the microscope type
    if intMaxImagePosition > 1

        % alternative code could be as follows:
        matRows = [];
        matCols = [];

        % for the cellWoRx and BD-Pathway
        if strcmpi(strMicroscopeType,'CW') || strcmpi(strMicroscopeType,'BD')
            for i = 1:matStitchDimensions(1);
                if ~mod(i,2)
                    matRows = [matRows,matStitchDimensions(2)-1:-1:0];
                else
                    matRows = [matRows,0:matStitchDimensions(2)-1];
                end
                matCols = [matCols,repmat((matStitchDimensions(1)-i),1,matStitchDimensions(2))];
            end
            matImageSnake = [matRows;matCols];

        % for the MD and CV7K
        elseif strcmpi(strMicroscopeType,'MD') || strcmpi(strMicroscopeType,'CV7K')
            for i = 1:matStitchDimensions(1);
                matRows = [matRows,0:matStitchDimensions(2)-1];
                matCols = [matCols,repmat(i-1,1,matStitchDimensions(2))];
            end
            matImageSnake = [matRows;matCols];
        end

        fprintf('%s: %d images per well, of type %s. cols = %d, rows = %d\n',mfilename,intMaxImagePosition, char(strMicroscopeType), max(matRows(:)+1),max(matCols(:)+1))
    end

end







function basedir=getbasedir(strRootPath)
% getlastdir retreives the base path excluding the last filename or
% directory of a path

    if nargin==0
        strRootPath = {'Z:\Data\Users\YF_DG\20080309165648_M1_080309_YF_DG_batch1_CP004-1de\BATCH\Measurements_Image_FileNames.mat'};
    end

    if ischar(strRootPath)
        basedir = doit(strRootPath);
    elseif iscell(strRootPath)
        basedir = cellfun(@doit,strRootPath,'UniformOutput',0);
    else
        error('unknown input type for getbasedir')
    end


end

function strOutput = doit(strInput)

    if strncmp(strInput,'\\',2)
        boolStartWithSlashSlash = 1;
        strInput = strInput(3:end);
    else
        boolStartWithSlashSlash = 0;
    end

    strInput = strrep(strInput,strcat(filesep,filesep),filesep);
    
    strOutput = '';
    matFilesepIndices = strfind(strInput, filesep);

    if isempty(matFilesepIndices) && strcmp(filesep,'\')
        matFilesepIndices = strfind(strInput, '/');        
    elseif isempty(matFilesepIndices) && strcmp(filesep,'/')
        matFilesepIndices = strfind(strInput, '\');
    end

    intPathLength = size(strInput,2);
    if matFilesepIndices(end) == intPathLength
        basedir = strInput(1:matFilesepIndices(end-1));
    else
        basedir = strInput(1:matFilesepIndices(end)); 
    end

    if boolStartWithSlashSlash
        strOutput = ['\\',basedir];
    else
        strOutput = basedir;
    end
end


function [result,matBin,descriptor]=histogram2(x,y,descriptor)
%HISTOGRAM2 Computes the two dimensional frequency histogram of two
%           row vectors x and y.
%   [RESULT,DESCRIPTOR] = HISTOGRAM2(X,Y) or
%   [RESULT,DESCRIPTOR] = HISTOGRAM2(X,Y,DESCRIPTOR) or
%where
%   DESCRIPTOR = [LOWERX,UPPERX,NCELLX;
%                 LOWERY,UPPERY,NCELLY]
%
%   RESULT       : A matrix vector containing the histogram
%   DESCRIPTOR   : The used descriptor
%
%   X,Y          : The row vectors to be analyzed
%   DESCRIPTOR   : The descriptor of the histogram
%     LOWER?     : The lowerbound of the ? dimension of the histogram
%     UPPER?     : The upperbound of the ? dimension of the histogram
%     NCELL?     : The number of cells of the ? dimension of the histogram
%
%   See also: http://www.cs.rug.nl/~rudy/matlab/

%   R. Moddemeijer 
%   Copyright (c) by R. Moddemeijer
%   $Revision: 1.2 $  $Date: 2001/02/05 09:54:29 $

    if nargin <1
       disp('Usage: RESULT = HISTOGRAM2(X,Y)')
       disp('       RESULT = HISTOGRAM2(X,Y,DESCRIPTOR)')
       disp('Where: DESCRIPTOR = [LOWERX,UPPERX,NCELLX;')
       disp('                     LOWERY,UPPERY,NCELLY]')
       return
    end

    % Some initial tests on the input arguments

    [NRowX,NColX]=size(x);

    if NRowX~=1
      error('Invalid dimension of X');
    end;

    [NRowY,NColY]=size(y);

    if NRowY~=1
      error('Invalid dimension of Y');
    end;

    if NColX~=NColY
      error('Unequal length of X and Y');
    end;

    if nargin>3
      error('Too many arguments');
    end;

    if nargin==2
      minx=min(x);
      maxx=max(x);
      deltax=(maxx-minx)/(length(x)-1);
    %   ncellx=ceil(length(x)^(1/3));
      ncellx=floor(length(x)^(1/3));
      miny=min(y);
      maxy=max(y);
      deltay=(maxy-miny)/(length(y)-1);
      ncelly=ncellx;
      descriptor=[minx-deltax/2,maxx+deltax/2,ncellx;miny-deltay/2,maxy+deltay/2,ncelly];
    end;

    lowerx=descriptor(1,1);
    upperx=descriptor(1,2);
    ncellx=descriptor(1,3);
    lowery=descriptor(2,1);
    uppery=descriptor(2,2);
    ncelly=descriptor(2,3);

    if ncellx<1 
      error('Invalid number of cells in X dimension')
    end;

    if ncelly<1 
      error('Invalid number of cells in Y dimension')
    end;

    if upperx<=lowerx
      error('Invalid bounds in X dimension')
    end;

    if uppery<=lowery
      error('Invalid bounds in Y dimension')
    end;

    result(1:ncellx,1:ncelly)=0;

    xx=round( (x-lowerx)/(upperx-lowerx)*ncellx + 1/2 );
    yy=round( (y-lowery)/(uppery-lowery)*ncelly + 1/2 );

    % [BS-HACK] as in histc i want a vector with for each object the
    % histogram-bin it is in. let's use vector-indexing for 2D object
    matBin(1:NColX)=NaN;
    matBinID=reshape(1:(ncellx*ncelly),ncellx,ncelly);

    for n=1:NColX
      indexx=xx(n);
      indexy=yy(n);
      if indexx >= 1 & indexx <= ncellx & indexy >= 1 & indexy <= ncelly
        result(indexx,indexy)=result(indexx,indexy)+1;
        % [BS-HACK]
        matBin(n)=matBinID(indexx,indexy);    
      end;
    end;
end



function strPrintName = gcf2pdf(figHandle, strRootPath, strFigureName,varargin)
% help for gcf2pdf 
%
% gcf2pdf converts current figure into a vectorized pdf file.
%
% usage: 
% strPrintName = gcf2pdf(strRootPath,strFigureName,varargin)

    % activate current figure
    figure(figHandle)
    drawnow

    strPrintName = '';

%     %%% allow gcf2pdf to also see hidden windows
    strOrigSetting = get(0,'ShowHiddenHandles');
    if strcmp(strOrigSetting,'off')
        set(0,'ShowHiddenHandles','on');
    end    
    
    set(0,'ShowHiddenHandles','on')
    set(gcf,'RendererMode','manual')    
    if ~isempty(find(strcmpi(varargin,'zbuffer'), 1))    
        disp(sprintf('%s: using renderer: Z-Buffer',mfilename))
        set(gcf,'Renderer','zbuffer')
    else
        disp(sprintf('%s: using renderer: Painters',mfilename))        
        set(gcf,'Renderer','painters')        
    end

    % prepare for pdf printing
    scrsz = [1,1,1920,1200];
    set(gcf, 'Position', [1 scrsz(4) scrsz(3) scrsz(4)]);     
    shading interp
    set(gcf,'PaperPositionMode','auto', 'PaperOrientation','landscape')
    set(gcf, 'PaperUnits', 'normalized'); 
    printposition = [0 .2 1 .8];

    set(gcf,'PaperPosition', printposition)
    
    % check varargin for papersize references
    cellstrPaperSizes = {'A0','A1','A2','A3','A4','A5'};
    if sum(ismember(upper(varargin),cellstrPaperSizes))>0
        % IF NOHEADER CONTAINS PAPER SIZE REFERENCE, USE THIS
        strPaperSize = varargin(ismember(upper(varargin),cellstrPaperSizes));
        strPaperSize = strPaperSize{1};
        set(gcf, 'PaperType', strPaperSize);
    else
        % DEFAULT TO A4 PAPER SIZE
        set(gcf, 'PaperType', 'A4');   
    end           
    
    orient landscape

    drawnow
    

    % in case overwrite is set as an option
    strPrintName = fullfile(strRootPath,strFigureName);
    
    disp(sprintf('%s: storing %s',mfilename,strPrintName))
    print(gcf,'-dpdf',strPrintName);

    if strcmp(strOrigSetting,'off')
        set(0,'ShowHiddenHandles','off');
    end    
    
end


function boolIsDir=isdir(x)
    boolIsDir=fileattrib(char(x));
end


function B = detectEdge(A)
% edge detection for segmentation images.
%
% Usage:
%
%   B = detectEdge(A)
%
% by Berend Snijder
    if max(A(:)) > 0
        h = fspecial('laplacian',1);
        matColormap = colormap('Jet');
        B = rgb2gray(label2rgb(A,smoothcolormap(matColormap,max(A(:))),'k','shuffle'));
        B = imfilter(B,h);
        B = B>0;
    else
        B=A;
    end

end


function matColorMap = smoothcolormap(colormapSettings,intNumOfRows)

    if nargin<=1
        intNumOfRows = 255;
    end

    matColorMap = imresize(colormapSettings,[intNumOfRows,3],'lanczos2');
    matColorMap(matColorMap>1)=1;
    matColorMap(matColorMap<0)=0;

end



function c = BSsmooth(y)
    % moving average with 5% span over the data.
    y = y(:); % force column vector
    span = 5; % fix rolling average width
    n = length(y);
    span = min(span,n);
    width = span-1+mod(span,2); % force it to be odd

    % moving average calculation
    c = filter(ones(width,1)/width,1,y);
    cbegin = cumsum(y(1:width-2));
    cbegin = cbegin(1:2:end)./(1:2:(width-2))';
    cend = cumsum(y(n:-1:n-width+3));
    cend = cend(end:-2:1)./(width-2:-2:1)';
    c = [cbegin;c(width:end);cend];
end