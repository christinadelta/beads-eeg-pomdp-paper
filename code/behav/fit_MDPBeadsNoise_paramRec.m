function fMDP_output = fit_MDPBeadsNoise_paramRec(Xfit,R,sequences,choicedata)

% run model whith Xfit = estimated cost-sample

% unpack the R structure
% unpack the R structure
R.sample    = R.cs;
beta        = Xfit(1);
maxDraws    = 10;
Ntrials     = size(sequences,1);
Qsad        = NaN*ones(Ntrials, maxDraws, 3);
% ll          = 0;
urns        = R.urntype;

% loop over sequences and run backward induction
for trial = 1:Ntrials

    nDrawsActual    = maxDraws; % now we need the entire draw length
    trialDraws      = sequences(trial,1:nDrawsActual);  % Include bead at guess (e.g., 3)
    nDraws          = length(trialDraws);  % 3
    trial_urn       = urns(trial);  % Sets q = 0.8 or 0.6

    % I need to work on the coding of the beads in the sequence.. wil come
    % back to this
    % Recode 2s to 0s for backward induction 
    if trial_urn == 0  % Majority colour is green

        % Original: 1 = green (majority), 2 = blue (minority)
        trialDraws(trialDraws == 2)     = 0;  % 2 (blue) becomes 0, 1 (green) stays 1
        
    else  % trial_urn == 1, majority colour is blue

        % Original: 1 = blue (majority), 2 = green (minority)
        trialDraws                      = (trialDraws == 2);   % 2 (green) becomes 1, 1 (blue) becomes 0
    end
    
    % Recode 2s to 0s for backward induction 
    % trialDraws(trialDraws  == 2) = 0;

    % Run backward utility
    for draw = 1 : nDraws

        Qsad(trial, draw, 1:3)  = backWardUtility(trialDraws, draw, maxDraws, R)';
        vVec                    = Qsad(trial, draw, 1:3);
        cprob(trial, draw, :)   = exp(beta*vVec)./sum(exp(beta*vVec));
        
    end % end of draws
end % end of trials loop


% use the choice probabilities to compute choices (in a probabilistic sampling decision model)
N                               = 1000;
[cprob_samples,model_urnchoice] = MDPProbSampling(cprob,maxDraws,Ntrials,N);

% map blue and green urns to majority urns and compute performance 
for dri = 1:Ntrials
    trial_urn = urns(dri);  % true urn: 1 = blue majority, 0 = green majority
    
    % Model’s choice
    if trial_urn == 1  % blue is majority
        if model_urnchoice(dri) == 2  % model picks majority (blue)
            correctResp(dri)    = 1;
        else  % model picks minority (green)
            correctResp(dri)    = 0;
        end
    
    elseif trial_urn == 0  % trial_urn == 0, green is majority
        
        if model_urnchoice(dri) == 1  % model picks majority (green)
            correctResp(dri)    = 1;
        else  % model picks minority (blue)
            correctResp(dri)    = 0;
        end
    end
end

% store model output structure
fMDP_output.performance      = mean(correctResp);
fMDP_output.avsamples        = mean(cprob_samples)-1;
fMDP_output.samples          = cprob_samples-1;
fMDP_output.choice           = correctResp;
fMDP_output.actionVals       = Qsad;
fMDP_output.modelurns        = model_urnchoice;

end % fit mdp function

%%%%%%%%%%%
function Qsa = backWardUtility(drawSequence, draw, maxDraws, R)

utility = zeros(maxDraws, maxDraws+1);

ng = sum(drawSequence(1:draw));

for drawi = maxDraws : -1 : (draw + 1)
        
    [utility] = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, R);
    
end
    
Qsa = actionValueBeads(utility, R, draw, ng, draw, maxDraws);
end

%%%%%%%%%%%%%%%%
function utility_t = stateUtilityBeads(utility, drawi, draw, maxDraws, ng, R)

utility_t = zeros(maxDraws, maxDraws+1);

futureDraws = drawi - draw;

ndf = drawi;

for greenDraws = 0 : futureDraws
    
    ngf = ng + greenDraws;

    Qsa = actionValueBeads(utility, R, ndf, ngf, drawi, maxDraws);

    utility_t(ndf, ngf+1) = max(Qsa);        
    
end
end

%%%%%%%%%%%%%%%%%%
function Qsa = actionValueBeads(utility, R, nd, ng, drawi, maxDraws)

pg = PG(R.thisq, nd, ng);

pb = 1 - pg;

QG = R.correct*pg + R.error*pb;
QB = R.correct*pb + R.error*pg;

if drawi < maxDraws

    QD = R.sample + pb*((1-R.thisq)*utility(nd+1, ng+1+1) +   (R.thisq)*(utility(nd+1, ng+1))) + ...
                    pg*(  (R.thisq)*utility(nd+1, ng+1+1) + (1-R.thisq)*(utility(nd+1, ng+1)));
                
else
    
    QD = 0;
    
end

Qsa = [QG; QB; QD];

end 

%%%%%%%%%%%%%%
function p = PG(q, nd, ng)

p = 1/(1 + (q/(1-q))^(nd-2*ng));

end

%%%%%%%%%%%%%%%%%
function [model_avsamples, model_urnchoice] = MDPProbSampling(cprob, maxDraws, Ntrials, N)
    
    % Initialise outputs
    model_avsamples     = zeros(Ntrials, 1);  % Average draws per trial
    model_urnchoice     = zeros(Ntrials, 1);  % Most frequent urn choice per trial

    for sequence = 1:Ntrials
        
        choiceProbs = squeeze(cprob(sequence, :, :));  % 10 x 3 array of probabilities
        
        % Compute cumulative probabilities for the three actions
        cum_probs   = cumsum(choiceProbs, 2);  % Cumulative probabilities across actions
        
        % Generate random samples from a uniform distribution
        test        = rand(N, maxDraws);
        
        % Temporary array to store urn choices for this sequence
        urn_choices_this_test   = zeros(1, N);
        samples_this_test       = zeros(1, N);  % Pre-allocate to avoid dynamic resizing

        for iteration = 1:size(test,1)
            
            % Find first sequence position where either Urn1 or Urn2 is chosen
            stop_idx            = find(test(iteration,:) < cum_probs(:,1)' | test(iteration,:) < cum_probs(:,2)', 1, 'first');
            
            if ~isempty(stop_idx)
                
                % A choice was made before maxDraws
                samples_this_test(iteration) = stop_idx;
                stop_draw = stop_idx;
            else
                % No choice made within maxDraws, force a choice
                samples_this_test(iteration)    = maxDraws;
                stop_draw                       = maxDraws;
                
                % Choose urn with highest probability at maxDraws
                if choiceProbs(maxDraws, 1) > choiceProbs(maxDraws, 2)
                    urn_choices_this_test(iteration) = 1;  % Urn1
                else
                    urn_choices_this_test(iteration) = 2;  % Urn2
                end

                continue;  % Skip to next iteration since no further check needed
            end

            % Determine which urn was chosen at the stopping draw
            if test(iteration, stop_draw) < choiceProbs(stop_draw, 1)
                urn_choices_this_test(iteration) = 1;  % Urn1
            else 
                urn_choices_this_test(iteration) = 2;  % Urn2
            end
        end

        model_avsamples(sequence, 1) = mean(samples_this_test); % Average sample over iterations
        model_urnchoice(sequence, 1) = mode(urn_choices_this_test); % Most frequent urn choice for this sequence
    end
end