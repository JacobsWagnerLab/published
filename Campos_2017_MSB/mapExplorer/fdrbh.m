function [h, fdr] = fdrbh(p, q, method)
%% -------------------------------------------------------------------------
% function [h, fdr] = fdrbh(pvals, q, method)
% @auther: Manuel Campos
% 
% @date: August 14 2016
% 
% @copyright 2015-2016 Yale University
%=========================================================================
% ********************** input **********************
%pvals:     p-values to correct
%q:         desired FDR threshold. Default is set at 0.05
%method:	{'indep' | 'dep' } Implement Benjamini and Hockberg FDR
%           correction under the assumption of independence or dependence
%           of teh tests leading to the different p-values.
%
% ********************** output **********************
%h:         logical vector of FDR corrected pvalues below the threshold q
%fdr:       Corrected p-values
%
%=========================================================================
% This script executes the Benjamini & Hochberg (1995) (method 'indep')
% or the Benjamini & Yekutieli (2001) procedure for controlling the false
% discovery rate (FDR) of a family of hypothesis tests. FDR is the expected
% proportion of rejected hypotheses that are mistakenly rejected 
% (i.e., the null hypothesis is actually true for those tests). 
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
if min(p(:))<0 || max(p(:))>1
    error('Invalid p-values');
end

% Set the default level of significance at 5%
if nargin<2
    q = 0.05;
end
% Set the default procedure to the Benjamini & Hochberg assuming
% inpendences between tests.
if nargin<3
    method = 'indep';
end

switch method
    case 'indep'
        m = numel(p);
        [p_ord, idx] = sort(p);
        fdr_ord =  p_ord(:) .* (m./(1:m))';
        % Running min in reverse order (in-place)
        bioinfoprivate.cumminmaxmex(fdr_ord,'min','reverse');
        fdr(idx) = fdr_ord;
    % If the dependencies between the tests is takien into consideration by
    % setting the method to 'dep', apply the Benjamini & Yekutieli
    % procedure.
    case 'dep'
        m = numel(p);
        corm = m*sum(1./(1:m));
        [p_ord, idx] = sort(p);
        fdr_ord =  p_ord(:) .* (corm./(1:m))';
        % Running min in reverse order (in-place)
        bioinfoprivate.cumminmaxmex(fdr_ord,'min','reverse');
        fdr(idx) = fdr_ord;
    otherwise
        error('Improper method definition');
end
h = fdr<=q;
end

