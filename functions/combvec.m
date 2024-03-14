%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMBVEC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = combvec(input)
%Changed slightly from the source to accept cell array
%
%
%  Example
%
%    a{1} = [1 2];
%    a{2} = [3 4; 3 4];
%    a3 = CombVec(a)
%    a3 =
%        1     2     1     2
%        3     3     4     4
%        3     3     4     4

if isempty(input)
    out = [];
else
    %input = num2cell(input,2);
    out = input{1};
    for i=2:length(input)
        cur = input{i};
        out = [copyb(out,size(cur,2)); copyi(cur,size(out,2))];
    end
end

%=========================================================
    function b = copyb(mat,s)
        
        [~,mc] = size(mat);
        inds    = 1:mc;
        inds    = inds(ones(s,1),:).';
        b       = mat(:,inds(:));
    end
%=========================================================
    function b = copyi(mat,s)
        
        [~,mc] = size(mat);
        inds    = 1:mc;
        inds    = inds(ones(s,1),:);
        b       = mat(:,inds(:));
    end
end%function combvec