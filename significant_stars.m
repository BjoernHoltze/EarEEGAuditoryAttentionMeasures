function stars = significant_stars(p_value)
%% returns the number of stars according to the p-value
% input:    p_value     [double] p-values (range: from 0 to 1]
% output:   stars           [string] string containing the right amount of stars
%
% author: Björn Holtze
% date: 16.10.2020

    if p_value <= 0.001
        stars = '***';
    elseif p_value <= 0.01
        stars = '**';
    elseif p_value <= 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end
end

