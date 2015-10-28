%
% Calculate the discrete Frechet distance (the coupling measure)
% between two paths P and Q.
%
% Input: P and Q
%        Arrays where the rows are the waypoints, and columns are 
%        x,y,z values of each waypoint.
%
% Output: cm The coupling measure
%         cm_seq The coupling sequence
%
%
% David Butterworth, 2015
%
% Algorithm by:
% T. Eiter and H. Mannila, "Computing Discrete Frechet Distance",
% Technical Report CD-TR 94/64, Christian Doppler Laboratory for Expert
% Systems, Vienna University of Technology, 1994
%
% Back-tracking algorithm by:
% Zachary Danziger, 2011
%

function [cm, cm_seq] = DiscreteFrechetDistance(P, Q)

p = size(P,1); % Number of points in P
q = size(Q,1); % Number of points in Q

% ca: array [1..p, 1..q]
% (all values must be -1 at start)
ca = -1 .* ones(p,q);

% function d(i,j)
% Calculate Euclidean distance between two points u_i, v_i
% on the paths P and Q.
%(also called L2 or Pythagorean distance)
d = @(u,v) (sqrt(sum((u-v).^2)));

% function c(i,j)
% Calculate the coupling measure between two paths P and Q
% by iterating over the waypoints i = 1 to p on path P, j = 1 to q on path Q.
function ret = c(i,j)
    if (ca(i,j) > -1)
        ret = ca(i,j);
    elseif (i == 1) && (j == 1)
        ca(i,j) = d(P(1,:),Q(1,:));
        ret = ca(i,j);
    elseif (i > 1) && (j == 1)
        ca(i,j) = max( c(i-1,1), d(P(i,:),Q(1,:)) );
        ret = ca(i,j);
    elseif (i == 1) && (j > 1)
        ca(i,j) = max( c(1,j-1), d(P(1,:),Q(j,:)) );
        ret = ca(i,j);
    elseif (i > 1) && (j > 1)
        ca(i,j) = max( min([c(i-1,j), c(i-1,j-1), c(i,j-1)]), d(P(i,:),Q(j,:)) );
        ret = ca(i,j);
    else
        ca(i,j) = inf;
    end
end

% Return the coupling measure
cm = c(p, q);

% If DiscreteFrechetDist() is called with a second output argument,
% calculate the coupling sequence,
if (nargout == 2)
    cm_seq = getCouplingSequence(p, q, ca);
end

end

% Calculate the coupling sequence, which is the sequence of steps
% along each path that result in the coupling measure.
function cm_seq = getCouplingSequence(p, q, ca)
    cm_seq = zeros(p+q+1, 2);
    
    % Create matrix where first row and first column are inf,
    % with matrix ca contained within.
    padded_ca = [inf*ones(1,q+1); ...
                 [inf*ones(p,1), ca]];
    Pi = p + 1;
    Qi = q + 1;
    count = 1;
    
    % Iterate from p+1 down to 3, and q+1 down to 3
    while ((Pi ~= 2) || (Qi ~= 2))

        % Step down the gradient of matrix padded_ca
        [min_value, min_idx] = min([padded_ca(Pi-1,Qi), padded_ca(Pi-1,Qi-1), padded_ca(Pi,Qi-1)]);
        if (min_idx == 1)
            cm_seq(count,:) = [Pi-1, Qi];
            Pi = Pi - 1;
        elseif (min_idx == 2)
            cm_seq(count,:) = [Pi-1, Qi-1];
            Pi = Pi - 1;
            Qi = Qi - 1;
        elseif (min_idx == 3)
            cm_seq(count,:) = [Pi, Qi-1];
            Qi = Qi - 1;
        end
        count = count + 1;
    end
    
    % Subtract 1 (the padding value) from the indices in cm_seq
    for i = 1:numel(cm_seq)
        if (cm_seq(i) > 0)
            cm_seq(i) = cm_seq(i) - 1;
        end
    end
    
    % Find the index of the last non-zero value in cm_seq
    last_value_idx = find(cm_seq(:,1)==0,1,'first') - 1;
    
    % Get the non-zero values
    cm_seq = [ cm_seq(1:last_value_idx,:) ];

    % Flip order of rows from bottom-to-top
    cm_seq = flipud(cm_seq);
    
    % Add the last point of P and Q
    cm_seq = [cm_seq; ...
              p, q];
end
