%% Seaweed ANN (MLP)
%% Works best for cases {a12,f12,all} using [h1,h2] = [25,25], a=0.0001, L=1/30
clear all
close all
clc

%load seaweed spectroscopy data 
Data = load("Negin_TimeTempDepencdency_TPC_WithAbsorbance.txt");

Data_FucusVes_WETOH = Data(1:9,:);
Data_FucusVes_W = Data(10:18,:);
Data_AscoNod_WETOH = Data(19:27,:);
Data_AscoNod_W = Data(28:end,:);
Data_all = Data(:,:);


% Which data to use? 
    % Ascophyllum Nodosum with HC+W = a1
    % Ascophyllum Nodosum with HC+W+ETOH = a2
    % Fucus Vesicolusus with HC+W = f1
    % Fucus Vesicolusus with HC+W+ETOH = f2
Seaweed_Solvent_Combo = 'all';         %   a1, a2, f1, f2, or all.

% hidden layer size (h)
h = [25,25] ;
%Training algorithm
trainAlg = 'trainbr';
% Learning rate a
a = 0.00001;
%number of epochs
epochs = 3000;
% Input delay:
d1 = 1;
% Feedback delay:
d2 = 1;
% Interpolation ratio:
ratio = 1/(30);
L=ratio;
% Training:Testing Division
cut = 0.6;
% Interpolation method and colours: 
InpMethod = 'Spl';    %Replace this with whichever function abbreviation (Mak, Lin, Cub, Spl, or MultSpl)
% plot colours
colour = 'b';

switch d2
    case 2
        aa=2;
    case 1
        aa=0;
end
switch Seaweed_Solvent_Combo
    case 'f1'
        Data = Data_FucusVes_W;
        Title = 'TPC for Fucus Vesicolusus using HC+W Solvent';
        subTitle1 = 'All Datapoints';
        subTitle2 = 'Predicted Datapoints';
        colour = 'b';
    case 'f2' 
        Data = Data_FucusVes_WETOH;
        Title = 'TPC for Fucus Vesicolusus using HC+W+ETOH Solvent';
        subTitle1 = 'All Datapoints';
        subTitle2 = 'Predicted Datapoints';
        colour = 'b';
    case 'a1'
        Data = Data_AscoNod_W;
        Title = 'TPC for Ascophyllum Nodosus using HC+W Solvent';
        subTitle1 = 'All Datapoints';
        subTitle2 = 'Predicted Datapoints';   
        colour = 'r';
    case 'a2'
        Data = Data_AscoNod_WETOH;
        Title = 'TPC for Ascophyllum Nodosus using HC+W+ETOH Solvent';
        subTitle1 = 'All Datapoints';
        subTitle2 = 'Predicted Datapoints';
        colour = 'r';
    case 'all'
        Data = Data_all;
        Title = 'TPC for Ascophyllum Nodosus and Fucus Vesicolusus using HC+W+ETOH or HC+W Solvent';
        subTitle1 = 'All Datapoints';
        subTitle2 = 'Predicted Datapoints';
        colour = 'm';
    case 'a12'
        Data = vertcat(Data_AscoNod_W,Data_AscoNod_WETOH);
        Title = 'TPC for Ascophyllum Nodosus using data from Both Solvents';
        subTitle1 = 'All Datapoints';
        subTitle2 = 'Predicted Datapoints';
        colour = 'r';
    case 'f12'
        Data = vertcat(Data_FucusVes_W,Data_FucusVes_WETOH);
        Title = 'TPC for Fucus Vesicolusus using both Solvents';
        subTitle1 = 'All Datapoints';
        subTitle2 = 'Predicted Datapoints';
        colour = 'b';
end

Data_preInterp = Data;

[entries,attributes] = size(Data);
entries_breakpoint1 = round(entries*cut); %this is cutting out % of entries
NumNonInterp = length(Data_preInterp(round(entries_breakpoint1)+1:end,1));


%%%%%%%%%%%%%%%Interpolation after data division%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 0.01;
x = 1:length(Data(1:round(entries_breakpoint1)));
switch InpMethod
    case 'Mak'
        method = 'makima'
        colour2 = '[0.82,0.00,1.00]'
    case 'Lin'
        method = 'linear'
        colour2 = '[0.87,0.88,0.00]'
    case 'Cub'
        method = 'cubic'
        colour2 = 'b'
    case 'Spl'
        method = 'spline'
        colour2 = '[0.95,0.52,0.00]'
    case 'MultSpl'
        colour2 = '[0.90,0.49,0.00]'
end
%%
switch InpMethod
    case 'MultSpl'
        % Define the query points (for interpolation)
        xq = linspace(1, 24, 24*L); % Interpolation points (e.g., 100 points between 1 and 24)
        
        % Initialize the interpolated results
        interpolated_Data = zeros(length(xq), size(Data, 2));
        
        % Perform cubic spline interpolation for each variable
        for col = 1:size(Data, 2)
            % Interpolate using csaps (Cubic Smoothing Splines)
            % The second argument (s) is the smoothing parameter (0 = least smooth (exact Interpolation), 1 = interpolation)
            s = 1;
            spline_fn = csaps(x, Data(:, col), s); 
            
            % Evaluate the spline at the query points
            interpolated_data(:, col) = fnval(spline_fn, xq);
        end
        % % Plot the results for visualization
        % figure;
        % for col = 1:size(Data, 2)
        %     subplot(2, 2, col);
        %     plot(x, Data(:, col), 'o', 'MarkerFaceColor', 'r', 'DisplayName', 'Original Data');
        %     hold on;
        %     plot(xq, interpolated_Data(:, col), 'b-', 'DisplayName', 'Interpolated Data');
        %     title(['Variable ', num2str(col)]);
        %     legend;
        %     grid on;
        % end
    otherwise
        x = [1:round(entries_breakpoint1)]';
        F = griddedInterpolant(x,Data(1:round(entries_breakpoint1),:), method);
        qx = 1:ratio:length(Data(1:entries_breakpoint1,1));        
        Vq = F(qx);
        Data_Int = Vq;        
        Data = vertcat(Data_Int,Data(round(entries_breakpoint1):end,:));
end
% Assign input and output data
[entries2,attributes] = size(Data);
entries_breakpoint = round(entries2*cut); %this is cutting out % of entries
Data_inputs = Data(:,1:attributes-1); 
Data_output = Data(:, attributes); 
interpCut = length(Data_output)-NumNonInterp;

%%
% Create a NAR neural network
net = narxnet(1:d1, 1:d2, h,'open', trainAlg);  % Modify the architecture as needed
                                               %narxnet(1:numInputDelays, 1:numFeedbackDelays, hiddenlayersize(50), feedbackMode, trainFunction)
                                               %  See also TRAINGDM, TRAINGDA, TRAINGDX, TRAINLM, TRAINRP,
                                               %                     TRAINCGF, TRAINCGB, TRAINBFG, TRAINCGP, TRAINOSS.

% % For segregated data division

% Based on cut percentage
% net.divideFcn = 'divideind';  % You can use other division functions
% net.divideParam.trainInd = 1:round(0.5*entries_breakpoint); %for training 
% net.divideParam.valInd = round(0.5*entries_breakpoint):entries_breakpoint;  % for validating
% net.divideParam.testInd = entries_breakpoint:length(Data);  %for testing

% Based on number of non interpolated data points
% net.divideFcn = 'divideind';  % You can use other division functions
% net.divideParam.trainInd = 1:round(length(Data_output(1:end))-NumNonInterp-aa); %for training 
% net.divideParam.valInd = round(length(Data_output(1:end))-NumNonInterp-aa):(length(Data_output(1:end))-round(0.5*NumNonInterp)-aa);  % for validating
% net.divideParam.testInd = (length(Data_output)-round(0.5*NumNonInterp)-aa):length(Data_output);  %for testing

net.divideFcn = 'divideind';  % You can use other division functions
net.divideParam.trainInd = 1:0.5*round(length(Data_output(1:end))-NumNonInterp-aa); %for training 
net.divideParam.valInd = 0.5*round(length(Data_output(1:end))-NumNonInterp-aa):(length(Data_output(1:end))-NumNonInterp-aa);  % for validating
net.divideParam.testInd = (length(Data_output)-round(NumNonInterp)-aa):length(Data_output);  %for testing

% Configure the network
net.layers{1}.transferFcn = 'tansig';  % You can choose other activation functions

net.layers{2}.transferFcn = 'tansig';  % Output layer activation function




 %%%%%%%%%%%%%%%%%%%%%%%%% Training Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set training parameters
net.trainParam.epochs = epochs;        % Maximum number of training epochs
% net.trainParam.goal = 10^(-99);         % Performance goal
net.trainParam.max_fail = epochs;        % Maximum number of validation failures
net.trainParam.min_grad = 10^(-9999);     % Minimum gradient for convergence
% net.trainParam.mu = 0.00000001;           % Initial mu (if using 'trainlm' or traingdx)
net.trainParam.lr = a;           % Learning rate (if using 'trainlm' or traingdx)
% net.trainParam.show = 10;           % Epochs between displays
net.trainParam.showCommandLine = false;   % Display training progress in command line
net.trainParam.showWindow = true; % Display training progress in a separate window
net.performParam.normalization;




% Prepare input time series
inputTimeSeries = tonndata(Data_inputs, false, false);
targetTimeSeries = tonndata(Data_output, false, false);

% Prepare input and target time series with delays
[Data_inputs, inputStates, layerStates, Data_output] = preparets(net, inputTimeSeries, {}, targetTimeSeries);


% Train the NARX neural network
[net, tr] = train(net, Data_inputs, Data_output, inputStates, layerStates);

%%%%%%%%% Make predictions on the training set
predTrainData_outputs = net(Data_inputs, inputStates, layerStates);



predictedOutput = cell2mat(predTrainData_outputs);
expectedOutput = cell2mat(Data_output);
error = (-predictedOutput + expectedOutput)./expectedOutput;



%%
D_inputs = cell2mat(Data_inputs)';
Temp_axis = D_inputs(:,2);

switch Seaweed_Solvent_Combo
    case {'all','a12','f12'}
         % Plot the actual vs. predicted values
        figure('Name', 'model and error');
        fontname('Times New Roman');
        % subplot(2,1,1)
        plot(Temp_axis((round(interpCut)):end), cell2mat(Data_output((round(interpCut)):end)'), 'kd',LineWidth=1);
        hold on;
        plot(Temp_axis((round(interpCut)):end), cell2mat(predTrainData_outputs((round(interpCut)):end)'), Color=colour, Marker='o', MarkerSize=8, LineWidth=lw);
        fontname('Times New Roman');
        legend('Actual', 'Predicted', 'Position',[0.717023807548341 0.697148487959284 0.17392857340404 0.0754761920202346]);
        xlabel('Temperature (\circC)','FontSize',12);
        ylabel('Total Phenolic Content (mg/L)','FontSize',12);
        title(Title, subTitle2,'FontSize',10);
        % subplot(2,1,2)
        %     plot(error((round(interpCut+0.5*NumNonInterp)):end)', 'kd')
        %     % ylim(subplot(2,1,2),[-1 0.5])
        %     fontname('Times New Roman');
        % ylabel("Relative Error",'FontSize',12)
    otherwise
       % Plot the actual vs. predicted values
        figure('Name', 'model and error');
        fontname('Times New Roman');
        % subplot(2,1,1)
        plot(D_inputs((round(interpCut+0.5*NumNonInterp)):end,2),cell2mat(Data_output(round(interpCut+0.5*NumNonInterp):end)'), 'k:x',LineWidth=1);
        hold on;
        plot(D_inputs((round(interpCut+0.5*NumNonInterp)):end,2), cell2mat(predTrainData_outputs((round(interpCut+0.5*NumNonInterp)):end)'), Color=colour, Marker='o', MarkerSize=8, LineWidth=lw);
        fontname('Times New Roman');
        legend('Actual', 'Predicted', 'Position',[0.717023807548341 0.697148487959284 0.17392857340404 0.0754761920202346]);
        xlabel('Temp (C)','FontSize',12);
        ylabel('Total Phenolic Content (mg/L)','FontSize',12);
        title(Title, subTitle2,'FontSize',10);
        % subplot(2,1,2)
        %     plot(error(entries_breakpoint -aa:end)', 'k-')
        %     % ylim(subplot(2,1,2),[-1 0.5])
        %     fontname('Times New Roman');
        % ylabel("Relative Error",'FontSize',12)
end

switch Seaweed_Solvent_Combo
    case {'all','a12','f12'}
         % Plot the actual vs. predicted values
        figure('Name', 'model and error');
        fontname('Times New Roman');
        % subplot(2,1,1)
        plot(cell2mat(Data_output(:)'), 'k:x',LineWidth=1);
        hold on;
        plot(cell2mat(predTrainData_outputs(:)'), Color=colour, Marker='o', MarkerSize=8, LineWidth=lw);
        fontname('Times New Roman');
        legend('Actual', 'Predicted', 'Position',[0.717023807548341 0.697148487959284 0.17392857340404 0.0754761920202346]);
        xlabel('Data Point','FontSize',12);
        ylabel('Total Phenolic Content (mg/L)','FontSize',12);
        title(Title, subTitle1,'FontSize',10);
        % subplot(2,1,2)
        %     plot(error(:)', 'k-')
        %     % ylim(subplot(2,1,2),[-1 0.5])
        % ylabel("Relative Error",'FontSize',12)
    otherwise
       % Plot the actual vs. predicted values
        figure('Name', 'model and error');
        fontname('Times New Roman');
        % subplot(2,1,1)
        plot(D_inputs(:,1),cell2mat(Data_output(:)'), 'k:',LineWidth=2);
        hold on;
        plot(D_inputs(:,1), cell2mat(predTrainData_outputs(:)'), Color=colour, Marker='o', MarkerSize=8, LineWidth=lw);
        fontname('Times New Roman');
        legend('Actual', 'Predicted', 'Position',[0.717023807548341 0.697148487959284 0.17392857340404 0.0754761920202346]);
        xlabel('Time (min)','FontSize',12);
        ylabel('Total Phenolic Content (mg/L)','FontSize',12);
        title(Title, subTitle1,'FontSize',10);
        % subplot(2,1,2)
        %     plot(error(:)', 'k-')
        %     % ylim(subplot(2,1,2),[-1 0.5])
        % ylabel("Relative Error",'FontSize',12)
end


%% 

NARX = cell2mat(predTrainData_outputs((entries_breakpoint -aa):end))';
errNARX = error((round(interpCut)):end)';


PredOutput_test = cell2mat(predTrainData_outputs((round(interpCut)):end)');
expectedOutput_test = expectedOutput((round(interpCut)):end)';

% PredOutput_test = cell2mat(predTrainData_outputs((entries_breakpoint -aa:end)));
% expectedOutput_test = expectedOutput((entries_breakpoint -aa:end));

mdl2 = fitlm(PredOutput_test,expectedOutput_test);
disp(mdl2)
disp(['Standard deviation or error: ', num2str(std(errNARX))])
disp(['Mean Absolute Error: ', num2str(mean(abs(errNARX)))])




% NARX = cell2mat(predTrainData_outputs((entries_breakpoint -aa):end))';
% errNARX = error(entries_breakpoint -aa:end)';
% 
% PredOutput_test = cell2mat(predTrainData_outputs(end-round(cut*NumNonInterp):end)');
% expectedOutput_test = expectedOutput(end-round(cut*NumNonInterp):end)';
% 
% % PredOutput_test = cell2mat(predTrainData_outputs((entries_breakpoint -aa:end)));
% % expectedOutput_test = expectedOutput((entries_breakpoint -aa:end));
% 
% mdl2 = fitlm(PredOutput_test,expectedOutput_test);
% disp(mdl2)
% disp(['Standard deviation or error: ', num2str(std(errNARX))])
% disp(['Mean Absolute Error: ', num2str(mean(abs(errNARX)))])


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All possible adjustable training parameters:%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 'trainFcn': Specifies the training function.
% Typical values: 'trainlm' (Levenberg-Marquardt), 'trainbfg' (BFGS
% quasi-Newton), 'trainrp' (Rprop), etc. 
% 
% 'trainParam.epochs': Maximum
% number of training epochs.
% % Typical values: 100, 200, 500, etc. 
% 
% 'trainParam.goal': Performance goal.
% Typical values: 1e-6, 1e-5, etc. 
% 
% 'trainParam.max_fail': Maximum number of validation failures. 
% Typical values: 6, 10, etc. 
% 
% 'trainParam.min_grad': Minimum gradient for
% convergence. 
% Typical values: 1e-6, 1e-5, etc. 
% 
% 'trainParam.lr': Learning rate for the
% Levenberg-Marquardt algorithm. 
% Typical values: 0.01, 0.001, etc. (if using 'trainlm') 
% 
% 'trainParam.mu':
% Initial mu (Levenberg-Marquardt parameter). 
% Typical values: 0.01, 0.001, etc. 
% 
% (if using 'trainlm') 'trainParam.show':
% Epochs between displays. 
% Typical values: 10, 25, etc. 
% 
% 'trainParam.showCommandLine': Display
% training progress in command line. 
% Typical values: true or false. 
% 
% 'trainParam.showWindow': Display training
% progress in a separate window.
% Typical values: true or false.


%List of training functions:
%   trainlm         %   
%   trainbfg        %
%   trainrp         %%%
%   trainscg        %
%   traincgb        %%%
%   traincgf        %%%
%   traincgp        %
%   trainposs       %
%   traingdx        %



% Possible transfer functions:
% 
% compet - Competitive transfer function.
% elliotsig - Elliot sigmoid transfer function.
% hardlim - Positive hard limit transfer function.
% hardlims - Symmetric hard limit transfer function.
% logsig - Logarithmic sigmoid transfer function.
% netinv - Inverse transfer function.
% poslin - Positive linear transfer function1.
% purelin - Linear transfer function1.
% radbas - Radial basis transfer function1.
% radbasn - Radial basis normalized transfer function1.
% satlin - Positive saturating linear transfer function1.
% satlins - Symmetric saturating linear transfer function1.
% softmax - Soft max transfer function1.
% tansig - Symmetric sigmoid transfer function1.
% tribas - Triangular basis transfer function1.



