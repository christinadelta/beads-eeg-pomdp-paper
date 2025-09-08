function ll = mdpBeadsCost_paramRec(params,R,sequences,choicedata)

% unpack the R structure
R.sample    = params(1);
beta        = R.beta;
maxDraws    = 10;
Ntrials     = size(sequences,1);
Qsad        = NaN*ones(Ntrials, maxDraws, 3);
ll          = 0;
urns        = R.urntype;

% loop over sequences and run backward induction
for trial = 1:Ntrials

    trialChoices    = choicedata{1, trial};
    nDrawsActual    = sum(trialChoices(:, 3));  % Number of draws (e.g., 2)
    trialDraws      = sequences(trial,1:nDrawsActual + 1);  % Include bead at guess (e.g., 3)
    nDraws          = length(trialDraws);  % 3
    trial_urn       = urns(trial);  % Sets q = 0.8 or 0.6

    % I need to work on the coding of the beads in the sequence.. wil come
    % back to this
    % Recode 2s to 0s for backward induction 
    if trial_urn == 0  % Majority colour is green
        
        % Original: 1 = green (majority), 2 = blue (minority)
        trialDraws(trialDraws == 2) = 0;  % 2 (blue) becomes 0, 1 (green) stays 1
    
    else  % trial_urn == 1, majority colour is blue
        
        % Original: 1 = blue (majority), 2 = green (minority)
        trialDraws = (trialDraws == 2);   % 2 (green) becomes 1, 1 (blue) becomes 0
    end

    % Recode 2s to 0s for backward induction 
    % trialDraws(trialDraws  == 2) = 0;

    % Run backward utility
    for draw = 1 : nDraws

        Qsad(trial, draw, 1:3)  = backWardUtility(trialDraws, draw, maxDraws, R)'; 
        vVec                    = Qsad(trial, draw, 1:3);
        cprob(trial, draw, :)   = exp(beta*vVec)./sum(exp(beta*vVec));
    
    end

    % I need to make sure that the urn choice is coded correctly for the
    % computations of ll... if subject chose the blue urn = 2, if subject
    % chose the green urn = 1;
    lastRow = trialChoices(end, :);  % Get last row: [1 0 0] or [0 1 0]
    if lastRow(1) == 1
        seqChoice = 2;  % Blue
    elseif lastRow(2) == 1
        seqChoice = 1;  % Green
    else
        error('Invalid final choice');  % Shouldn’t happen with your data
    end

    % seqChoice = (choiceData{trial}(end) == 2) + 1; 

    if nDraws-1 > 0 & nDraws < maxDraws
        ll = ll - sum(log(squeeze(cprob(trial, nDraws-1, 3)))) - log(squeeze(cprob(trial, nDraws, seqChoice)));
    elseif nDraws < maxDraws
        ll = ll - log(squeeze(cprob(trial, nDraws, seqChoice)));
    end  

end % end of trials loop

end % end of MDP function 

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