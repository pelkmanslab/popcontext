function handles = PopulationContextCorrect(handles)

% Help for the Population Context Correct module:
% Category: Measurement
%
%
% SHORT DESCRIPTION:
%
% This module calculates the population context-corrected measurement as
% described in Snijder & Sacher, et al, in revision, 2011
%
% This module creates a bootstrapped multidimensional bin model of your
% selected measurement against selected population context parameters, and
% corrects your single-cell measurements with the trends described in the
% model. 
%
% The module will store 3 measurements: (1) "Raw", (2) "Predicted" and (3)
% "Corrected". "Raw" is the raw measurement, "Predicted" the average value
% observed by the bin-model for cells of that population context, and
% "Corrected" the value corrected with the selected method (by default a
% subtraction of raw - predicted).
%
%
% *************************************************************************
% Here is more detailed information on the different options:
%
%
% From which object does the measurement to correct come from?
%  --> Select the object from which the measurement was created that you
%  want to correct.
%
% Which category of measurements would you like to use?
%  --> Select the measurement category the measurement you would like to
%  correct comes from. See calculateMath module help for more description
%  of the logic with which measurement selection is done.
%
% Which feature do you want to use? (feature number or name)
%  --> Select the feature number or name of the measurement you want to
%  correct.
%
% Some features come from images, if so, from which?
%  --> If your measurement is derived from an image, like a TEXTURE or
%  INTENSITY measurement, select which image the measurement comes from.
%
% Some features use size scale (TEXTURE OR NEIGHBORS) or number of bins
% (RADIALDISTRIBUTION). If relevant, which did you use? 
%  --> Id.
%
% Do you want to log10-transform the measurement?
%  --> log10 transform of your selected measurement can be good for the
%  log-normalization of for instance intensity values.
% 
% Do you want to correct with Local Cell Density
%  --> Select Yes if you want to include the LCD as predictor in your bin
%  model. Generally, keep the number of predictors in your model to a
%  minimum. At least two predictors are required for the module to run.
%
% Do you want to correct with Object Size
%  --> Select Yes if you want to include the object size as predictor in
%  your bin model. 
%
% Do you want to correct with Edge
%  --> Select Yes if you want to include the cell islet edge measurement as
%  predictor in your bin model.
%
% Do you want to correct with Population Size
%  --> Select Yes if you want to include the population size measurement as
%  predictor in your bin model. Only switch this option on if you have
%  large fields of view and if you have many (>20) wells analyzed in the
%  CellProfiler run.
%
% Do you want to correct with Nearest Neighbor Distance
%  --> Select Yes if you want to include the nearest neighbor distance
%  measurement as predictor in your bin model.
%
% Do you want to correct with Distance From Islet Edge
%  --> Select Yes if you want to include the shortest distance of each cell
%  to the cell islet edge measurement as predictor in your bin model.
%
% How do you want to correct?
%  --> Allows you to select how you want to calculate the single-cell
%  corrected result. The default "subtract (measured - predicted)" is a
%  subtraction between the measured value and the bin-model described
%  value. Use this for instance for log10 intensity values. Other
%  possibilities include "log2 ratio (measured / predicted)" and "ratio
%  (measured / predicted)",
%
% If you have multiple MeasurePopulationContext modules, which module
% measurements (first, second, etc.) should we use? 
%  --> You can have multiple MeasurePopulationContext modules in the same
%  pipeline, and the results will depend on the settings value. This will
%  select the nth module of your pipeline for correction. Leave to value
%  "1" if you have only one such module in your pipeline. 
%
% Do you want save a PDF of the output figure?
%  --> If you select yes, a PDF will be created of the output figure and
%  stored with a timestamp in your default output directory. PDF creation
%  requires the figure to close. If the figure was open at the time of PDF
%  creation, the module will redraw the original figure.
%
%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% PopulationContextCorrect was developed & written by Berend Snijder,
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

%textVAR01 = From which object does the measurement to correct come from?
%choiceVAR01 = Nuclei
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
strObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = Which category of measurements would you like to use?
%inputtypeVAR02 = popupmenu category
Category = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Which feature do you want to use? (feature number or name)
%defaultVAR03 = 1
%inputtypeVAR03 = popupmenu measurement
FeatureNumber = handles.Settings.VariableValues{CurrentModuleNum,3};

%textVAR04 = Some features come from images, if so, from which?
%infotypeVAR04 = imagegroup
%inputtypeVAR04 = popupmenu
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Some features use size scale (TEXTURE OR NEIGHBORS) or number of bins (RADIALDISTRIBUTION). If relevant, which did you use?
%defaultVAR05 = 1
%inputtypeVAR05 = popupmenu scale
SizeScale = str2double(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Do you want to log10-transform the measurement?
%choiceVAR06 = No
%choiceVAR06 = Yes
%inputtypeVAR06 = popupmenu
strLog10TransformMeasurement = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Do you want to correct with Local Cell Density
%choiceVAR07 = Yes
%choiceVAR07 = No
%inputtypeVAR07 = popupmenu
strUseLCD = char(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = Do you want to correct with Object Size
%choiceVAR08 = Yes
%choiceVAR08 = No
%inputtypeVAR08 = popupmenu
strUseAREA = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = Do you want to correct with Edge
%choiceVAR09 = Yes
%choiceVAR09 = No
%inputtypeVAR09 = popupmenu
strUseEDGE = char(handles.Settings.VariableValues{CurrentModuleNum,9});

%textVAR10 = Do you want to correct with Population Size
%choiceVAR10 = No
%choiceVAR10 = Yes
%inputtypeVAR10 = popupmenu
strUsePOPSIZE = char(handles.Settings.VariableValues{CurrentModuleNum,10});

%textVAR11 = Do you want to correct with Nearest Neighbor Distance
%choiceVAR11 = No
%choiceVAR11 = Yes
%inputtypeVAR11 = popupmenu
strUseNEIGHBOR = char(handles.Settings.VariableValues{CurrentModuleNum,11});

%textVAR12 = Do you want to correct with Distance From Islet Edge
%choiceVAR12 = No
%choiceVAR12 = Yes
%inputtypeVAR12 = popupmenu
strUseDISTFROMEDGE = char(handles.Settings.VariableValues{CurrentModuleNum,12});

%textVAR13 = How do you want to correct?
%choiceVAR13 = subtract (measured - predicted)
%choiceVAR13 = log2 ratio (measured / predicted)
%choiceVAR13 = ratio (measured / predicted)
%inputtypeVAR13 = popupmenu
strCorrectionMethod = char(handles.Settings.VariableValues{CurrentModuleNum,13});

%textVAR14 = If you have multiple MeasurePopulationContext modules, which module measurements (first, second, etc.) should we use?
%defaultVAR14 = 1
%inputtypeVAR14 = popupmenu scale
intPopcornModule = str2double(handles.Settings.VariableValues{CurrentModuleNum,14});

%textVAR15 = Do you want save a PDF of the output figure?
%choiceVAR15 = No
%choiceVAR15 = Yes
%inputtypeVAR15 = popupmenu
strCreatePDF = char(handles.Settings.VariableValues{CurrentModuleNum,15});


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
    CPfigure(handles,'Text',ThisModuleFigureNumber);
    colormap('jet')
end


% on the first run, display a message on the output figure
if (SetBeingAnalyzed == handles.Current.StartingImageSet)

    % draw/activate figure if user did not close the window
    if any(findobj == ThisModuleFigureNumber)

        % leave a message
        currentfig = CPfigure(handles,'Text',ThisModuleFigureNumber);
        TextString = 'The PopulationContextCorrect module only runs at the last image cycle.';
        uicontrol(currentfig,'style','text','units','normalized','fontsize',handles.Preferences.FontSize,'HorizontalAlignment','center','string',TextString,'position',[.05 .45 .9 .1],'BackgroundColor',[.7 .7 .9],'tag','TextUIControl');

    end
    
    % do a quick check to see if a MeasurePopulationContext module is
    % present 
    matPopCornModulesIX = find(strcmp(handles.Settings.ModuleNames,'MeasurePopulationContext'));
    
    % if there is no population context measurement module, we should alert
    % the user and stop..
    if isempty(matPopCornModulesIX)
        CPmsgbox(sprintf('In order for the PopulationContextCorrect module to run we need at least a single MeasurePopulationContext module to be present before this module in the pipeline.'),'Missing MeasurePopulationContext','help')
        % should we error or just quit? just quit.
        error('In order for the PopulationContextCorrect module to run we need at least a single MeasurePopulationContext module to be present before this module in the pipeline.');
    elseif matPopCornModulesIX(intPopcornModule) > CurrentModuleNum
        CPmsgbox(sprintf('In order for the PopulationContextCorrect module to run we need at least a single MeasurePopulationContext module to be present before this module in the pipeline.'),'Missing MeasurePopulationContext','help')
        % should we error or just quit? just quit.
        error('In order for the PopulationContextCorrect module to run we need at least a single MeasurePopulationContext module to be present before this module in the pipeline.');
    end        
    
    

elseif SetBeingAnalyzed == handles.Current.NumberOfImageSets

    % do correction here... Few lines below adapted from the CalculateMath
    % module.

    try
        strFeatureName = CPgetfeaturenamesfromnumbers(handles, strObjectName, Category, FeatureNumber, ImageName, SizeScale);

        cellMeasurement = handles.Measurements.(strObjectName).(strFeatureName);
        strMeasurementName = sprintf('%s_%s',strObjectName,strFeatureName);
        
    catch objFoo
        error([objFoo.message '  Image processing was canceled in the ', ModuleName, ...
            ' module (#' num2str(CurrentModuleNum) ...
            ') because an error ocurred when retrieving the data.  '...
            'Likely the category of measurement you chose, ',...
            Category, ', was not available for ', ...
            strObjectName,' with feature number ' FeatureNumber ...
            ', possibly specific to image ''' ImageName ''' and/or ' ...
            'Texture Scale = ' num2str(SizeScale) '.']);
    end

    % get the object count per image
    matObjectCount = cat(1,handles.Measurements.Image.(['Count_',strObjectName]){:});
    
    % initialize the output as NaNs
    handles.Measurements.(strObjectName).(sprintf('PopContext_Raw_%s',strFeatureName)) = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)';
    handles.Measurements.(strObjectName).(sprintf('PopContext_Predicted_%s',strFeatureName)) = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)';
    handles.Measurements.(strObjectName).(sprintf('PopContext_Corrected_%s',strFeatureName)) = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)';
    
    % generate explicit meta data, i.e. image-index and object-index
    cellNucleiMetaData = arrayfun(@(x,y) cat(2,repmat(y,x,1), (1:x)'), matObjectCount,(1:numel(matObjectCount))', 'UniformOutput',false);

    % get all nuclei positions and meta data for current well
    matMeasurement = cat(1,cellMeasurement{:});
    matNucleiMetaData = cat(1,cellNucleiMetaData{:});
    
    % log10 transform if requested
    if strcmpi(strLog10TransformMeasurement,'yes')
        matMeasurement = log10(matMeasurement);
        if any(isinf(matMeasurement))
            CPmsgbox(sprintf('The log10 transform has led to %d Inf measurements... Setting them to NaN.',sum(isinf(matMeasurement))),'Setting Infs to NaN','help')
            matMeasurement(isinf(matMeasurement)) = nan;
        end
    end

    
    
    % get the population context measurements
    matPopCornModulesIX = find(strcmp(handles.Settings.ModuleNames,'MeasurePopulationContext'));
    
    % if there is no population context measurement module, we should alert
    % the user and stop..
    if isempty(matPopCornModulesIX)
        CPmsgbox(sprintf('In order for the PopulationContextCorrect module to run we need at least a single MeasurePopulationContext module to be present before this module in the pipeline.'),'Missing MeasurePopulationContext','help')
        % should we error or just quit? just quit.
        error('In order for the PopulationContextCorrect module to run we need at least a single MeasurePopulationContext module to be present before this module in the pipeline.');
    else   
        matPopCornModulesIX = matPopCornModulesIX(intPopcornModule);
    end
    
    % get settings from the population context correction module
    strPopcornObject = handles.Settings.VariableValues{matPopCornModulesIX,1};
    intCellDiameter = str2double(handles.Settings.VariableValues{matPopCornModulesIX,2});
    strWellImage = handles.Settings.VariableValues{matPopCornModulesIX,4};
    
    % measurements are stored with settings in them...
    strSettingDescription = sprintf('Per%s_Cell%d',strWellImage,intCellDiameter);
    
    % get the population desired context features and format the feature
    % list and final complete data...

    % first colom is the measurement that will be corrected, and add if
    % it's log10 transformed or not..
    matCompleteData = matMeasurement;
    if strcmpi(strLog10TransformMeasurement,'yes')
        cellstrFeatures{1} = ['log10 ',strMeasurementName];
    else
        cellstrFeatures{1} = strMeasurementName;
    end
    
    
    % add the features in order of preferred axes for drawing
    if strcmpi(strUseLCD,'Yes')
        matCompleteData = [matCompleteData, cat(1,handles.Measurements.(strPopcornObject).(['PopContext_',strSettingDescription,'_LocalCellDensity']){:})];
        cellstrFeatures{end+1} = 'Local Cell Density';
    end
    if strcmpi(strUseAREA,'Yes')    
        matCompleteData = [matCompleteData, cat(1,handles.Measurements.(strPopcornObject).(['PopContext_',strSettingDescription,'_Area']){:})];
        cellstrFeatures{end+1} = sprintf('%s Size',strPopcornObject);
    end
    if strcmpi(strUsePOPSIZE,'Yes')
        matCompleteData = [matCompleteData, cat(1,handles.Measurements.(strPopcornObject).(['PopContext_',strSettingDescription,'_PopulationSize']){:})];
        cellstrFeatures{end+1} = 'Population Size';
    end    
    if strcmpi(strUseEDGE,'Yes')
        matCompleteData = [matCompleteData, cat(1,handles.Measurements.(strPopcornObject).(['PopContext_',strSettingDescription,'_Edge']){:})];
        cellstrFeatures{end+1} = 'Cell Islet Edge';
    end
    if strcmpi(strUseNEIGHBOR,'Yes')
        matCompleteData = [matCompleteData, cat(1,handles.Measurements.(strPopcornObject).(['PopContext_',strSettingDescription,'_NearestNeighbor']){:})];
        cellstrFeatures{end+1} = 'Nearest Neighbor Distance';
    end
    if strcmpi(strUseDISTFROMEDGE,'Yes')
        matCompleteData = [matCompleteData, cat(1,handles.Measurements.(strPopcornObject).(['PopContext_',strSettingDescription,'_DistanceToEdge']){:})];
        cellstrFeatures{end+1} = 'Distance To Islet Edge';
    end

    
    % see if there is anything to build a model on.
    if numel(cellstrFeatures)==1
        CPmsgbox(sprintf('You must at least select one population context measurement to build a model with.'),'Incorrect settings','help')
        error('You must at least select one population context measurement to build a model with.');
    end

    % check for NANs and remove those...
    matNanIX = any([isnan(matMeasurement),isnan(matNucleiMetaData),isnan(matCompleteData)],2);
    if any(matNanIX)
        CPmsgbox(sprintf('Skipping %d cells with invalid measurements.',sum(matNanIX)),'Removing NaNs','help')
        matMeasurement(matNanIX,:) = [];
        matNucleiMetaData(matNanIX,:) = [];
        matCompleteData(matNanIX,:) = [];
    end
    
    % do bin correction, first column is readout, others are used for binning.
    % bootstrap over a bunch of different bin-numbers per dimension... and
    % keep cell number statistics.
    matCompleteDataBinReadout = NaN(size(matCompleteData,1),10);
    mat5thPercentileBinSize = NaN(1,10);
    intMinBins = 3;
    for i = 1:13
        [matCompleteDataBinReadout(:,i),~,~,~,matTensorTCN, matCompleteDataBinIndex] = doBinCorrection(matCompleteData, cellstrFeatures, intMinBins+(i-1), @nanmean);
        mat5thPercentileBinSize(i) = quantile(matTensorTCN(matCompleteDataBinIndex),0.05);
    end
    
    % see where the number of cells per bin (for the lower 5% of bins)
    % drops below 100. If it's always bad, revert to minimal number of bins
    intFinalBinNumber = find(mat5thPercentileBinSize>100,1,'last')+(intMinBins-1);
    if isempty(intFinalBinNumber);intFinalBinNumber=intMinBins;end

    % let's average models out from slightly undersampled (too few bins) to slightly
    % oversampled (too many bins), with more udnersampling. I.e., if
    % optimal is x, then let's average/smoothen out over x-2:x+1;
    intFinalBinNumberIX = intFinalBinNumber - (intMinBins-1);
    intBinRange = (max(intFinalBinNumberIX-2,1):min(intFinalBinNumberIX+1,size(matCompleteDataBinReadout,2)));
    matCompleteDataBinReadout = nanmean(matCompleteDataBinReadout(:,intBinRange),2);
    
    % recaulcalate model properties (bin edges, cell bin indices, etc.) for
    % final number of bins 
    [~,matBinEdges,matBinDimensions,~,~, matCompleteDataBinIndex, ~, matTCNPerBin] = doBinCorrection(matCompleteData, cellstrFeatures, intFinalBinNumber, @nanmean);

    % calculate correction according to user selected method
    switch lower(strCorrectionMethod)
        case 'log2 ratio (measured / predicted)'
            matCorrectedData = log2( matMeasurement ./ matCompleteDataBinReadout );
        case 'ratio (measured / predicted)'
            matCorrectedData = matMeasurement ./ matCompleteDataBinReadout;
        case 'subtract (measured - predicted)'
            matCorrectedData = matMeasurement - matCompleteDataBinReadout;
    end
    
    % now store the model prediction back in CP format at the respective
    % original image-index and object-index.
    for iImage = unique(matNucleiMetaData(:,1))'
        matImageIX = ismember(matNucleiMetaData(:,1),iImage);
        matObjectIX = matNucleiMetaData(matImageIX,2);
        handles.Measurements.(strObjectName).(sprintf('PopContext_Raw_%s',strFeatureName)){iImage}(matObjectIX) = matMeasurement(matImageIX);
        handles.Measurements.(strObjectName).(sprintf('PopContext_Predicted_%s',strFeatureName)){iImage}(matObjectIX) = matCompleteDataBinReadout(matImageIX);
        handles.Measurements.(strObjectName).(sprintf('PopContext_Corrected_%s',strFeatureName)){iImage}(matObjectIX) = matCorrectedData(matImageIX);
    end    
    
%     %%%%%%%%%%%%%%%%%%%
%     %%% DRAW FIGURE %%%
%     %%%%%%%%%%%%%%%%%%%

    % draw figure in any case, as we only draw it in the last run...
    if any(findobj == ThisModuleFigureNumber) || strcmpi(strCreatePDF,'yes')
        ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
        boolFigureWasClosed = CPwasfigureclosed(handles,ThisModuleFigureNumber);
        [currentfig] = drawFigure(handles,matCompleteDataBinIndex,matBinDimensions,matCompleteDataBinReadout,matBinEdges,intFinalBinNumber,cellstrFeatures,ThisModuleFigureNumber,strCorrectionMethod,matTCNPerBin,strPopcornObject,matCompleteData,matCorrectedData);
        if strcmpi(strCreatePDF,'yes')
            strPrintName = sprintf('%s_%s_ModelOverview.pdf',datestr(now,30),ModuleName);
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
        end
        if ~boolFigureWasClosed
            [currentfig] = drawFigure(handles,matCompleteDataBinIndex,matBinDimensions,matCompleteDataBinReadout,matBinEdges,intFinalBinNumber,cellstrFeatures,ThisModuleFigureNumber,strCorrectionMethod,matTCNPerBin,strPopcornObject,matCompleteData,matCorrectedData);
        end
    end
            
        
%     end    

    
    
    
end

% restore default colormap (we overwrite it in the beginning)
handles.Preferences.IntensityColorMap = strOldCMap;    


end % of function










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

function [currentfig] = drawFigure(handles,matCompleteDataBinIndex,matBinDimensions,matCompleteDataBinReadout,matBinEdges,intFinalBinNumber,cellstrFeatures,ThisModuleFigureNumber,strCorrectionMethod,matTCNPerBin,strPopcornObject,matCompleteData,matCorrectedData)


        
        % activate figure
        currentfig = CPfigure(handles,'Image',ThisModuleFigureNumber);
        % clear previous text
        clf(currentfig);

        % add fancy smooth plots :)
        foo = ind2sub2(matBinDimensions,matCompleteDataBinIndex);
        
        [matMeanSurf,matMeanSurfTCN] = wellfun(@nanmean,matCompleteDataBinReadout,foo(:,[1,2]));
        matMeanSurf(matMeanSurfTCN<4) = NaN;
        matSurfX = repmat(matBinEdges(:,2),[1,matBinDimensions(2)]);
        matSurfY = repmat(matBinEdges(1:matBinDimensions(2),3),[1,matBinDimensions(1)])';

        
        % calculate correlations between readout and predictors (before and
        % after correction)
        matAllCorrs = corr(matCompleteData);% correlations with raw readout
        matAllCorCorrs = corr([matCorrectedData, matCompleteData(:,2:end)]); % correlations with corrected readout
        
        % place some text (but not with uicontrol, as PDF creation can not
        % make this vectorized... (silly!)):
        % 0) (log10) measurement to correct 
        % 1) number of bins
        % 2) model parameters
        hAx = axes('position',[.05  .53  .42  .42]);
        axis off
        
        iLine = 0;
        iLineHeight = 0.095;
        TextString = sprintf('Population context model');
        text(0, 1-(iLine), TextString,'FontWeight','bold','FontSize',8);
        iLine = iLine + iLineHeight;

        TextString = sprintf('%d cells in %d bins (max %d bins/dim)', sum(matTCNPerBin), numel(matTCNPerBin),intFinalBinNumber);
        text(0, 1-(iLine), TextString,'FontSize',8);
        iLine = iLine + iLineHeight;        
        
        for i = 1:numel(cellstrFeatures)
            if i==1
                TextString = sprintf('Modeled: %s',cellstrFeatures{i});
            else
                TextString = sprintf('Predictor %d: %s (%d bins) (correlation raw=%.2f, corrected=%.2f)',i-1,cellstrFeatures{i},matBinDimensions(i-1),matAllCorrs(1,i),matAllCorCorrs(1,i));
            end
            text(0.05, 1-(iLine), TextString,'FontSize',8);
            iLine = iLine + iLineHeight;
        end
        TextString = sprintf('Corrected using %s',strCorrectionMethod);
        text(0.05, 1-(iLine), TextString,'FontSize',8);
        drawnow        
        
        
        if ~all(isnan(matMeanSurf(:)))
%             %%% A subplot of the figure window is set to display the original image.
%             hAx=subplot(5,5,[6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24],'Parent',ThisModuleFigureNumber);
            hAx=subplot(2,3,4,'Parent',ThisModuleFigureNumber);

            hold on
            surf(hAx,matSurfX(1:end-1,1:end-1),matSurfY(1:end-1,1:end-1),matMeanSurf(1:end-1,1:end-1),'linewidth',1,'FaceColor','interp','EdgeColor','none')% , 'FaceLighting','phong'
            % determine view depending on the correlation
            matDefaultView = [40, 50];
            if corr(nanmean(matMeanSurf,2),matSurfX(:,1),'rows','pairwise')
                matDefaultView(1) = -130;
            end
            view(hAx,matDefaultView);
            axis tight
            grid on
            set(hAx,'FontSize',7)
            xlabel(hAx,cellstrFeatures{2},'Interpreter','none')
            ylabel(hAx,cellstrFeatures{3},'Interpreter','none')
            zlabel(hAx,cellstrFeatures{1},'Interpreter','none')
            hold off
            drawnow

            % add a title/description below the surface plot
            title('Model projection on the first two model parameters','fontsize',8)
            %uicontrol(currentfig,'style','text','units','normalized','fontsize',8,'HorizontalAlignment','center','string',TextString,'position',[.03 .7 .9 .1],'tag','TextUIControl','FontWeight','bold','BackgroundColor','none');
            drawnow
        else
            % add a title/description below the surface plot
            TextString = 'No model created because there were too few cells';
            uicontrol(currentfig,'style','text','units','normalized','fontsize',8,'HorizontalAlignment','center','string',TextString,'position',[.03 .35 .9 .1],'tag','TextUIControl','FontWeight','bold');
            drawnow
            
        end

        % it's good to check the distribution of number of cells per bin, so let's
        % plot the (y=) number of cells per bin greater than x, where (x=) all
        % unique bin sizes.
        matUniqueTCN = unique(matTCNPerBin)';
        matCumSumTCN = NaN(size(matUniqueTCN));
        for i = 1:length(matUniqueTCN)
            matCumSumTCN(i) = sum(matTCNPerBin(matTCNPerBin>=matUniqueTCN(i)));
        end
        matCumSumTCN = 100*(matCumSumTCN/sum(matTCNPerBin));
        drawnow
        
        % plot bin size against cumulative #-cells for all bins bigger than x.
%         hAx=subplot(5,5,5,'Parent',ThisModuleFigureNumber);
        hAx=subplot(2,3,3,'Parent',ThisModuleFigureNumber);
        plot(hAx,matUniqueTCN,matCumSumTCN,'-ok')
        set(hAx,'fontsize',7)
        ylabel(hAx,sprintf('cumulative %% of cells\nin bins bigger than x'))
        xlabel(hAx,'bin size')
        title(hAx,sprintf('95%% of all cells is present in a bin bigger than %d cells.\n%.0f%% of all cells is present in a bin bigger than 20 cells.',matUniqueTCN(find(matCumSumTCN>=95,1,'last')),matCumSumTCN(find(matUniqueTCN>=20,1,'first'))))

        % plot total cell number surface plot
%         hAx=subplot(5,5,1,'Parent',ThisModuleFigureNumber);
        hAx=subplot(2,3,5,'Parent',ThisModuleFigureNumber);
        hold on
        surf(hAx,matSurfX(1:end-1,1:end-1),matSurfY(1:end-1,1:end-1),matMeanSurfTCN(1:end-1,1:end-1),'linewidth',1,'FaceColor','interp','EdgeColor','none')% , 'FaceLighting','phong'
        % determine view depending on the correlation
        matDefaultView = [40, 50];
        if corr(nanmean(matMeanSurf,2),matSurfX(:,1),'rows','pairwise')
            matDefaultView(1) = -130;
        end
        view(hAx,matDefaultView);
        axis tight
        grid on
        set(hAx,'FontSize',7)
        xlabel(hAx,cellstrFeatures{2},'Interpreter','none')
        ylabel(hAx,cellstrFeatures{3},'Interpreter','none')
        zlabel(hAx,sprintf('%s count',strPopcornObject),'Interpreter','none')
        title(sprintf('Number of %s',strPopcornObject),'fontsize',8)
        hold off
        drawnow

        % Pauli suggests to show a before and after correction plot. I
        % guess this would make sense against the strongest correlated
        % parameter? ...
        % actually a plotquant/plotmean type of figure would not be bad, as
        % there might be too many datapoints for a scatterplot.

        % find the strongest absolute correlated measure? (mutual
        % information or something might also be good, as non-linear
        % effects are also corrected for in most cases)
        % matAllCorrs = corr(matCompleteData);
        % slightly less like binary parameters for visualization...
        matCorrelationCorrectionFactor = 1-(0.1*(matBinDimensions==2));
        [~,iMax]=max(abs(matCorrelationCorrectionFactor .* matAllCorrs(1,2:end)));
        
        matYraw = matCompleteData(:,1);
        matYcor = matCorrectedData(:);
        matXData = matCompleteData(:,iMax+1);
        
        % discard any nans and infs datapoints
        matBadIX = any(isnan(matXData) | isnan(matYcor) | isnan(matYraw) | ...
                        isinf(matXData) | isinf(matYcor) | isinf(matYraw),2);
        matYraw(matBadIX) = [];
        matYcor(matBadIX) = [];
        matXData(matBadIX) = [];
            
        % show quantile bin mean +- stdev values
        iBins = 10;
        % calculate bin edges
        matBinEdges = quantile(matXData,linspace(0,1,iBins+1));
        matBinEdges(end) = matBinEdges(end)+0.01;
        % get bin indices
        [~,matXbin] = histc(matXData,matBinEdges);
        % calculate statistics
        matYbinRaw = NaN(iBins,2);
        matYbinCor = NaN(iBins,2);
        matYbinCount = NaN(iBins,1);
        for iBin = unique(matXbin)'
            matCurBinIX = matXbin==iBin;
            matYbinRaw(iBin,:) = [mean(matYraw(matCurBinIX)),std(matYraw(matCurBinIX))];
            matYbinCor(iBin,:) = [mean(matYcor(matCurBinIX)),std(matYcor(matCurBinIX))];
            matYbinCount(iBin) = sum(matCurBinIX);
        end

        % draw red and blue line and surface plots
        
        % skip datapoints with less than two cells or with NaNs
        matBadIX = matYbinCount<3 | any(isnan(matYbinRaw),2) | any(isnan(matYbinCor),2);

        figAx1=subplot(2,3,6,'Parent',ThisModuleFigureNumber);
        hold on
        % draw raw: red
        %fill([matBinEdges(~matBadIX),fliplr(matBinEdges(~matBadIX))],[(matYbinRaw(~matBadIX,1)-matYbinRaw(~matBadIX,2))',fliplr((matYbinRaw(~matBadIX,1)+matYbinRaw(~matBadIX,2))')],[1 0.5 0.5],'linestyle','none','FaceAlpha',0.75);
        plot(figAx1,matBinEdges(~matBadIX),matYbinRaw(~matBadIX,1)-0.5*matYbinRaw(~matBadIX,2),'--r');
        plot(figAx1,matBinEdges(~matBadIX),matYbinRaw(~matBadIX,1)+0.5*matYbinRaw(~matBadIX,2),'--r');
        plot(figAx1,matBinEdges(~matBadIX),matYbinRaw(~matBadIX,1),'-ok','MarkerFaceColor','r','MarkerSize',5);
        hold off
        
        % create second axis for overlay
        figAx2=axes();
        % draw corrected: green
        hold on
        % unfortunately the alpha face value does not work on PDF...
        %fill([matBinEdges(~matBadIX),fliplr(matBinEdges(~matBadIX))],[(matYbinCor(~matBadIX,1)-matYbinCor(~matBadIX,2))',fliplr((matYbinCor(~matBadIX,1)+matYbinCor(~matBadIX,2))')],[0.3 1 0.3],'linestyle','none','FaceAlpha',0.3);
        plot(figAx2,matBinEdges(~matBadIX),matYbinCor(~matBadIX,1)-0.5*matYbinCor(~matBadIX,2),'--g');
        plot(figAx2,matBinEdges(~matBadIX),matYbinCor(~matBadIX,1)+0.5*matYbinCor(~matBadIX,2),'--g');
        plot(figAx2,matBinEdges(~matBadIX),matYbinCor(~matBadIX,1),'-ok','MarkerFaceColor','g','MarkerSize',5);
        hold off
        matPositionFirstAxis = get(figAx1,'Position');
        set(figAx1,'Position',matPositionFirstAxis,'YAxisLocation','left','FontSize',7,'YColor',[0.5 0 0]);
        ylabel(figAx1,cellstrFeatures{1},'Interpreter','none');
        xlabel(figAx1,cellstrFeatures{iMax+1},'Interpreter','none');
        set(figAx2,'Position',matPositionFirstAxis,'YAxisLocation','right','Color','none','TickLength',[0 0],'YColor',[0 0.5 0],'XTick',[],'XTickLabel',[],'fontsize',7)
        ylabel(figAx2,sprintf('corrected %s',cellstrFeatures{1}),'Interpreter','none');
        title(sprintf('comparison raw (red) and corrected (green)\ncorrelation with ''%s'' raw: %.2f and corrected: %.2f',cellstrFeatures{iMax+1},corr(matXData,matYraw),corr(matXData,matYcor)),'Interpreter','none')
        hold off
        drawnow
        

%         % create new axis to overlay cell number data on...
%         plot(matBinEdges,matYbinCount,':ok','MarkerFaceColor',[1 1 1],'Color',[.3 .3 .3],'MarkerSize',5);    
%         ylabel(sprintf('%s count',strPopcornObject))
% 
        
        
        
end


function [matCompleteDataBinReadout,matBinEdges,matBinDimensions,matTensorII,matTensorTCN, matCompleteDataBinIndex, matIIPerBin, matTCNPerBin] = doBinCorrection(matCompleteData, cellstrFeatures, varargin)
%
% help for doBinCorrection
%
% Usage:
%
% [matCompleteDataBinReadout,matBinEdges,matBinDimensions,matTensorII,matTensorTCN, matCompleteDataBinIndex, matIIPerBin, matTCNPerBin] = doBinCorrection(matCompleteData, cellstrFeatures, varargin)
%
% [BS] Does the shizzle ultra fast. Mi nizzle extreme! That's what 3 years
% of PhD are for...
%
% Takes first column of matCompleteMetaData as readout of interest, and
% makes a multidimensional matrix from the remaining columns in
% matCompleteMetaData(:,2:end), using quantile binning, maximally 14 bins
% per dimension (less if bins have less than 14 unique discrete values)
%
% To get some shnizzle figures feedback in the mix!
%
% doBinCorrection(..., matBinEdges) 
%
% if matBinEdges has the right dimensions (columns equals the number of
% data columns, and rows is less than or equal to the max number of bins),
% it will overwrite automaticcally generated binning. 
%
% doBinCorrection(..., @function_handle) 
%
% will run a different function handle per bin, where @mean is default.
%
% doBinCorrection(..., boolColumnVectorMatrix) 
%
% where boolColumnVectorMatrix has as many rows as matCompleteData, only 1
% column, and is logical, true values mean that cell is left out of the
% tensor, but will get a predicted value.


if nargin==0
    strSettingsFile = npc('\\nas-biol-imsb-1\share-3-$\Data\Users\50K_final_reanalysis\ProbModel_Settings.txt');
    strRootPath = npc('\\nas-biol-imsb-1\share-3-$\Data\Users\50K_final_reanalysis\VSV_CNX');
    [matCompleteData, cellstrFeatures] = getRawProbModelData2(strRootPath,strSettingsFile);
end

% Parameter for the maximum number of bins per features, can be user
% supplied.
boolUserSuppliedMaxBins = any(cellfun(@(x) isnumeric(x) & numel(x)==1, varargin));
if boolUserSuppliedMaxBins
    % get the max number of bins from the most likely candidate in varargin
    intMaxBins = varargin{cellfun(@(x) isnumeric(x) & numel(x)==1, varargin)};
%     fprintf('%s: Max number of edges determined by user: %d.\n',mfilename,intMaxBins)
else
    % set max number of bins to default 12
    intMaxBins = 12;
end


% if varargin contains something that seriously looks like matBinEdges
% (same number of columns as matCompleteData, and very few rows, less than
% 24), user this instead of calculating your own
boolUserSuppliedBinEdges = any(cellfun(@(x) isnumeric(x) & size(x,2)==size(matCompleteData,2) & size(x,1)<=intMaxBins+1, varargin));

if boolUserSuppliedBinEdges
    % get the bin edges from the most likely candidata in varargin
%     fprintf('%s: Binning edges defined by user.\n',mfilename)
    matBinEdges = varargin{cellfun(@(x) size(x,2)==size(matCompleteData,2)  & size(x,1)<=intMaxBins+1, varargin)};
else
    % get quantile bin edges (i.e. quantile binning!)
    matBinEdges = quantile(matCompleteData,linspace(0,1,intMaxBins));
end

% check if the user supplied a different function handle to create the
% tensor with
boolUserSuppliedFunction = any(ismember(cellfun(@class,varargin,'UniformOutput',false),'function_handle'));
if boolUserSuppliedFunction
    strFunctionHandles = varargin{ismember(cellfun(@class,varargin,'UniformOutput',false),'function_handle')};
%     fprintf('%s: Running function ''%s'' per bin, as defined by user.\n',mfilename,char(strFunctionHandles))
else
    strFunctionHandles = @nanmean;
end

% check if the user supplied a different function handle to create the
% tensor with
boolUserSuppliedTensorExclusionData = any(cellfun(@(x) islogical(x) & size(x,1)==size(matCompleteData,1) & size(x,2)==1, varargin));
if boolUserSuppliedTensorExclusionData
    matCompleteDataToExclude = varargin{cellfun(@(x) islogical(x) & size(x,1)==size(matCompleteData,1) & size(x,2)==1, varargin)};
%     fprintf('%s: User supplied exclusion criteria for keeping %d%% of cells out of the tensor.\n',mfilename,round(100*nanmean(matCompleteDataToExclude)))
end


% do binning, assuming first column is readout!
numFeatures = size(matCompleteData,2)-1;
matBinDimensions = NaN(1,numFeatures);
matCompleteDataBinIX = NaN(size(matCompleteData,1),numFeatures);
for iFeature = 1:numFeatures
    
    % do binning. taking unique bin edges solves binary/non-binary readout
    % issue as well as bins with discrete [1,2,3...] numbers that are fewer
    % than the number of bins suggested :)
    matBinEdgesPerColumn = unique(matBinEdges(:,iFeature+1));
    [foo,matCompleteDataBinIX(:,iFeature)] = histc(matCompleteData(:,iFeature+1), matBinEdgesPerColumn);
    clear foo

    % Store the multidimensional matrix dimensions for current column
    matBinDimensions(iFeature) = length(matBinEdgesPerColumn);
    
end



% since unique-rows is slower than getting sub-indices for each bin and
% unqiue on those, let's do this.
if size(matCompleteDataBinIX,2)==1
    % note that if there is only one explaining param, the index equals the
    % sub-index.
    matSubIndices = matCompleteDataBinIX;
else
    matSubIndices = sub2ind2(matBinDimensions,matCompleteDataBinIX);
end

% sort completedata, completemetadata, binindices, and subindices so that
% subindices are sorted from lowest to highest values
[matSubIndices,matSortIX] = sort(matSubIndices);
% matCompleteDataBinIX=matCompleteDataBinIX(matSortIX,:);
% matCompleteMetaData=matCompleteMetaData(matSortIX,:);
matCompleteData=matCompleteData(matSortIX,:);
if boolUserSuppliedTensorExclusionData
    matCompleteDataToExclude=matCompleteDataToExclude(matSortIX,:);
end

% now get the unique values of the sorted subindices, and return the 'last'
% of each occuring unique subindex. this means we can then stepwise loop
% over all unique subindices and batch-wise calculate the readout values
% without ever having to do a 'find' command! sweet fastness.
[matUniqueBinSubIndices, matBinIX1, matBinIX2] = unique(matSubIndices,'last');

% loop over the each bin-value, and batch-wise calculate the corresponding
% readout of all cells belonging to that bin (i.e. mean infection index,
% and total cell number per bin)

% first check if the user already (perhaps) supplied a tensor for
% correction...
boolUserSuppliedTensor = any(cellfun(@(x) isnumeric(x) & isequal(size(x),size(matUniqueBinSubIndices)),varargin));
if boolUserSuppliedTensor
    matIIPerBin = varargin{cellfun(@(x) isnumeric(x) & isequal(size(x),size(matUniqueBinSubIndices)),varargin)};
%     fprintf('%s: Correcting using user supplied tensor.\n',mfilename,char(strFunctionHandles))
else
    % otherwise, initialize
    matIIPerBin = NaN(size(matUniqueBinSubIndices));
end

matTCNPerBin = NaN(size(matUniqueBinSubIndices));
matPreviousIX = 1;
for iBin = 1:size(matUniqueBinSubIndices,1)
    matCurrentIX = (matPreviousIX:matBinIX1(iBin));
    
    % if user supplied criteria to keep data out of the tensor, kick it out
    % here...
    if boolUserSuppliedTensorExclusionData
        matCurrentIX(matCompleteDataToExclude(matCurrentIX)) = [];
    end
    
    % if user requested a different function to be run...
    if ~boolUserSuppliedTensor
        if boolUserSuppliedFunction
            matIIPerBin(iBin) = strFunctionHandles(matCompleteData(matCurrentIX,1));
        else
            matIIPerBin(iBin) = mean(matCompleteData(matCurrentIX,1));
        end
    end
    matTCNPerBin(iBin) = length(matCurrentIX);
    matPreviousIX = matBinIX1(iBin)+1;
end

matCompleteDataBinReadout = matIIPerBin(matBinIX2);


% shuffle matCompleteData and back to original sorting
matOrigSortingIX = 1:size(matCompleteData,1);
matOrigSortingIX = matOrigSortingIX(matSortIX);
[foo,matOrigSortingIX] = sort(matOrigSortingIX);

clear foo
% (is there another way to sort back from original sorting?)
matCompleteDataBinReadout = matCompleteDataBinReadout(matOrigSortingIX,:);
matCompleteDataBinIndex = matUniqueBinSubIndices(matBinIX2(matOrigSortingIX,:));

% if it's requested, also return actual tensors with IIs and TCNs
if nargout>3
    % to actually make this a multidimensional matrix (tensor), do this:
    matTensorII = NaN(matBinDimensions);
    matTensorTCN = NaN(matBinDimensions);
    matTensorII(matUniqueBinSubIndices) = matIIPerBin;
    matTensorTCN(matUniqueBinSubIndices) = matTCNPerBin;
end


end



function matSubIndices = sub2ind2(matBinDimensions,matCompleteDataBinIX) %#ok<STOUT,INUSL>

    str2exec = 'matSubIndices = sub2ind(matBinDimensions';
    for i = 1:size(matCompleteDataBinIX,2)
        str2exec = [str2exec,sprintf(',matCompleteDataBinIX(:,%d)',i)]; %#ok<AGROW>
    end
    str2exec = [str2exec,');'];
    eval(str2exec);
    
end


function matCompleteDataBinIX = ind2sub2(matBinDimensions,matSubIndices) %#ok<STOUT,INUSD>

    str2exec = '[I1';
    for i = 2:size(matBinDimensions,2)
        str2exec = [str2exec,sprintf(',I%d',i)]; %#ok<AGROW>
    end
    str2exec = [str2exec,']'];
    
    eval(sprintf('%s = ind2sub(matBinDimensions,matSubIndices);',str2exec));
    eval(sprintf('matCompleteDataBinIX = %s;',str2exec));
    
end


function [dataOutput, matOutputTCN] = wellfun(fhandle, matCompleteData, matCompleteMetaData, isOutputCell)
% WELLFUN parses getRawProbModel2 data structures to evaluate any function
% handle on all cells per well.
%
%
% [dataOutput, matOutputTCN] = wellfun(fhandle, matCompleteData, matCompleteMetaData, isOutputCell)
%
% Evaluates the function handle of "fhandle" on all the datapoints in
% all columns of matCompleteData, per each well as indicated by
% matCompleteMetaData. The user must specify beforehand wether the
% dataOutput should be a cell array (if output per well contains more than
% one datapoint) or wether dataOutput can be a matrix (default if not
% specified).
%
% dataOutput will contain a matrix with per well and plate -index the
% output from fhandle run on all datapoints per well
%
% matOutputTCN will contain the number of datapoints per well
%
% If matCompleteMetaData only has two columns, the first is the
% well-row-number and the second is the well-column-number.
%
% If matCompleteMetaData contains 3 or more columns, the third is the plate
% number. In this case, all output will be 3D, of shape [row,column,plate].
%
% Example:
%
%    % settings file for getRawProbModelData2
%    strSettingsFile = npc('\\nas-biol-imsb-1\share-3-$\Data\Users\50K_final_reanalysis\ProbModel_Settings.txt');
%    % data path for getRawProbModelData2
%    strRootPath = '\\nas-biol-imsb-1\share-3-$\Data\Users\50K_final_reanalysis\VSV_CNX';
%    % get data from path with settingsfile
%    [matCompleteData, strFinalFieldName, matCompleteMetaData] = getRawProbModelData2(strRootPath,strSettingsFile);
%    % function handle
%    fhandle = @(x) median(sum(x));
%    % run wellfun
%    [dataOutput, matOutputTCN] = wellfun(fhandle, matCompleteData, matCompleteMetaData)
%
%
%   See also getRawProbModelData2 cellfun arrayfun
%
% Copyright: Berend Snijder, 2010. Diggety.

% default output is matrix, due to potential memory saving
if (nargin<4)
    isOutputCell = false;
end

% if no plate numbers are passed, process only well indices
intMetaDataColSize = size(matCompleteMetaData,2);
intMetaDataCols = 1:intMetaDataColSize;

% find out how many well rows & columns are present
matRowColPlateNums = max(matCompleteMetaData(:,intMetaDataCols));

% reformat well indices as plate indices, giving unique indices per well.
matWellSubIndices = sub2ind2(matRowColPlateNums,matCompleteMetaData(:,intMetaDataCols));

% sort data according to well indices
[matWellSubIndices, matSortIX] = sort(matWellSubIndices);

% sort input data
matCompleteData = matCompleteData(matSortIX,:);

% now get the unique values of the sorted subindices, and return the 'last'
% of each occuring unique subindex. this means we can then stepwise loop
% over all unique subindices and batch-wise calculate the readout values
% without ever having to do a 'find' command! sweet fastness.
[matUniqueWellSubIndices, matBinIX1] = unique(matWellSubIndices,'last');

% loop over the each well-index value, and batch-wise evaluate the function
% handle "fhandle" for all cells for that well index.

% output can be either a cell or a matrix, but must be a-priori indicated
% by user. defaults to matrix.
if isOutputCell
    dataOutput = cell(matRowColPlateNums);
else
    dataOutput = NaN(matRowColPlateNums);
end

% initialize datapoint count per well
matOutputTCN = NaN(matRowColPlateNums);

matPreviousIX = 1;
for iBin = 1:size(matUniqueWellSubIndices,1)
    % get data sub-indices
    matCurrentIX = (matPreviousIX:matBinIX1(iBin));
    
    % evaluate input function
    matWellOutput = feval(fhandle, matCompleteData(matCurrentIX,:));
    
    % do a check. if output of first well contains more than one datapoint,
    % set isOutputCell to true and notify user. Only do this check at the
    % first point.
    if iBin==1
        if ~isOutputCell && numel(matWellOutput)>1
            % fprintf('%s: for first well output evaluated in more than a single datapoint. setting output to cell array. please see ''help wellfun'' for more info\n',mfilename)
            isOutputCell  = true;
            dataOutput = cell(matRowColPlateNums);
        end
    end
    
    % put output in dataOutput
    if isOutputCell
        dataOutput{matUniqueWellSubIndices(iBin)} = matWellOutput;
    else
        dataOutput(matUniqueWellSubIndices(iBin)) = matWellOutput;
    end
    
    % count number of datapoints
    matOutputTCN(matUniqueWellSubIndices(iBin)) = numel(matCurrentIX);
    
    % increment data sub-indices to go to start of next data point
    matPreviousIX = matBinIX1(iBin)+1;
end

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