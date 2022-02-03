function [X, acc, p_rand] = check_classifier( num_trials, alpha, num_classes, varargin )
% CHECK_CLASSIFIER better than random?
%
% Calculates the minimum number of correct classifications or minimum
% classification accuracy for NC trials under the confidence limit alpha.
% If no output arguments are provided the results are printed in text form.
%
% function [X, acc] = check_classifier( num_trials, alpha, num_classes )
%
% Inputs:
%    num_trials : number of trials per class. A scalar corresponds to an
%                 equal number of trials for each class. In the 2 class case,
%                 num_trials can be a row vector with C elements that specifies
%                 the number of trials for each class seperately.
%                 Classifier output probabilities are assumed to be the
%                 same as trial probabilities.
%                 If num_trials is a matrix of size [2xC], the first row is
%                 the number of trials for each class, and the second row
%                 are the classifier output probabilities for each class.
%         alpha : confidence limit (default: 0.05)
%   num_classes : number of classes (default: 2 or number of elements in num_trials)
%
% Optional:
%        'plot'   : create histogram plot
%        'method' : how to obtain the confidence interval
%                   'nomial' : use the inverse of the binomial distribution.
%                   'cp'     : Clopper-Pearson interval
%                   'wald'   : Wald interval
%                   'awald'  : (default) Adjusted Wald interval
%                   'wilson' : Wilson score interval
%
% Outputs:
%     X : minimum number of correct classifications
%   acc : minimum accuracy
%
% ========================================================================
%  A single classification result can either be correct or incorrect. When
%  class c appears with a probabilty of p_c, the random classifier selects
%  this class with the same probabilty p_c (assumption). Thus the probabilty
%  of correctly classifying class c by chance is p_c*p_c. The total
%  probability of correct classification (i.e. the classification accuracy
%  due to chance, p0) is the sum of p_c*p_c for all c. Thus the random
%  classifier produces correct results that are binomially distributed with
%  probability p.
%  From the CDF of that binomial distribution, the treshold for non-random
%  classification can be determined at a given alpha-error.

if nargin < 1
	help check_classifier
	error( 'Not enough input arguments.' )
end

doplot = false;
method = 'awald';

for i = 1 : length(varargin)
    if strcmpi( varargin{i}, 'plot' )
        doplot = true;
    end
    if strcmpi( varargin{i}, 'method' )
        method = varargin{i+1};
        i = i+1;
    end
end

if nargin < 2
	alpha = 0.05;
end

if nargin < 3
	num_classes = max(2,length(num_trials));
end

if num_classes < 2
    error( 'Less than 2 classes does not make sense.' )
end

if isvector(num_trials)
    if length( num_trials ) == 1
        num_trials = repmat( num_trials, num_classes, 1 );
    end
    rate_outputs = num_trials;
else
    rate_outputs = num_trials(2,:);
    num_trials = num_trials(1,:);
end

num_trials = num_trials(:);
rate_outputs = rate_outputs(:) / sum( rate_outputs );

if num_classes ~= length( num_trials )
    error( 'num_trials must be scalar or a vector of length num_classes.' )
end

N = sum( num_trials );              % total number of trials
p_rand = sum( num_trials .* rate_outputs / N  );  % TP rate of random classifier
%p_rand = sum( (num_trials/N).^2 );  % TP rate of random classifier
%p_rand = sum( (num_trials/N).*ones(length(num_trials),1)/length(num_trials) );  % TP rate of random classifier
%p_rand = max( num_trials/N );

if strcmpi( method, 'nomial' )
    X_critical = binoinv( 1-alpha, N, p_rand ) + 1;
    acc_critical = X_critical / N;
elseif strcmpi( method, 'cp' )
    k = N * p_rand;
    f = finv( 1-alpha, 2*k+2, 2*n - 2*k );
    acc_critical = (k+1) * f / ( N - k + (k+1) * f );
    X_critical = ceil(N * acc_critical);
elseif strcmpi( method, 'wald' )
    z = norminv( 1-alpha, 0, 1 );
    acc_critical = p_rand + z * sqrt( p_rand * (1-p_rand) / N );
    X_critical = ceil(N * acc_critical);
elseif strcmpi( method, 'awald' )
    z = norminv( 1-alpha, 0, 1 );
    p = (( N * p_rand ) + 2) / (N + 4);
    acc_critical = p + z * sqrt( p * (1-p) / (N+4) );
    X_critical = ceil(N * acc_critical);
elseif strcmpi( method, 'wilson' )
    z = norminv( 1-alpha, 0, 1 );
    acc_critical = ( p_rand + z^2/(2*N) + z*sqrt( p_rand*(1-p_rand)/N + z^2/(4*N^2) ) ) / ( 1 + z^2/N );
    X_critical = ceil(N * acc_critical);
else
    error( ['Unknown confidence interval method: ' method] )
end


if doplot
    if( X_critical <= N )
        bar( 0:X_critical-1, pdf(1:X_critical), 'b' )
        hold on
        bar( X_critical:N, pdf(X_critical+1:N+1), 'r' )
        hold off
    else
        bar( 0:N, pdf, 'b' )
    end
    
    range = find(max(pdf)./pdf < 200);
    range(2:end-1) = [];
    axis( [range(1),range(2),0,max(pdf)] );
    
    xlabel( 'Number of correctly classified trials' )
    ylabel( 'Probability' )
    title( sprintf('%d-class Random Classification Limit (%d trials total, \\alpha=%f)', num_classes, N, alpha) )
end

if( X_critical > N )
    disp( sprintf('With this number of trials, no significant classification\nbetter than random is possible at alpha=%f.',alpha) )
    acc = -1;
    X = -1;
    return
end

if( nargout == 0 )
    disp( ['Number of classes               : ' num2str(num_classes)] )
    disp( ['Class propabilities             : [', num2str(num_trials'/N) ']'] )
    disp( ['Alpha                           : ' num2str(alpha)] )
    disp( ['Minimum correct classifications : ' num2str(X_critical)] )
    disp( ['Minimum classification accuracy : ' num2str(acc_critical)] )
	%disp( sprintf( 'A %d-class classifier must have at least %d correct\nclassifications out of %d trials to be better than\nrandom within a confidence limit of %f.', num_classes, X_critical, N, alpha ) )
	%disp( sprintf( 'This corresponds to a minimum accuracy of %f.', acc_critical ) )
else
	acc = acc_critical;
	X = X_critical;
end
