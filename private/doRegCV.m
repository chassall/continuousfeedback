function [allErrors,bestBeta] = doRegCV(data,X,regtype,condIs,breakpoints,lambdas,k)
%DOREGCV Regularization with cross-validation
%   Detailed explanation goes here (TODO)

% data: channels X samples
% X: design matrix, samples X regressors
% regtype: 'ident', 'onediff', 'twodiff', 'threediff' (see pinv_reg.m)
% condIs: points to include, e.g. {1:size(data,2)} for everything, or can
% split up for speed, e.g. {1:n, n+2:size(data,2)}
% breakpoints: last index of each waveform in X, e.g. [1000 2000] for
% signals that go from 1-1000, 1001-2000
% lambdas: reg. parameters to try
% k: number of folds for the cross-validation, e.g. 10
% Returns
% allErrors: sum of squared error, averaged across each sensor, for each
% lambda (you should plot this and it should be U-shaped)
% bestBeta: beta estimates using the best lambda (regressors X sensors)
% Example usage:
% [allErrors,bestBeta]  = doRegCV(data,X,'onediff',{1:size(X,2)},[1000 2000],[100 1000 10000 100000 1000000],10);

allErrors = nan(1,length(lambdas));
tempBetas = [];

% Loop through each lambda
disp('looping through lambdas');
for li = 1:length(lambdas)

    lambda = lambdas(li);
    disp(lambda);
    cv = cvpartition(size(data,2),'KFold',k); % Make CV train/test partitions
    theseErrors = nan(1,k);
    
    % Cross-validation loop (loop through folds)
    for ki = 1:k
        disp(ki);
        
        % TRAIN
        thisTrainingI = training(cv,ki);
        thisTestingI = test(cv,ki);

        % Loop through user-specified sections of DM
        for c = 1:length(condIs)

            % Compute pseudoinverse of design matrix with regularization
            thisPDM = pinv_reg(X(thisTrainingI,condIs{c}),lambda,regtype,breakpoints{c});
            
            % TODO: some issue when lambda = 0. Could be issue with pinv?

            % Compute beta
            thisBeta(condIs{c},:) = thisPDM * data(:,thisTrainingI)';
        end

        tempBetas(li,:,:) = thisBeta;
        
        % TEST
        thisResidual = data(:,thisTestingI)' - full(X(thisTestingI,:)) * thisBeta; % compute residual
        thisError = sum(thisResidual.*thisResidual); % error for each sensor
        theseErrors(ki) = mean(thisError); % mean across sensors
    end
    allErrors(li) = mean(theseErrors);
end

[~,bestI] = min(allErrors);
bestBeta = squeeze(tempBetas(bestI,:,:));

end