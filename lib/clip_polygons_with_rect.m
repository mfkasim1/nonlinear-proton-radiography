%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to clip 2 polygons with clipping polygon as rectangle.
% The algorithm uses Sutherland-Hodgman algorithm.
% Conditions:
%   * The clipping polygon must be rectangle and the subject polygon can be non-convex.
%   * All coordinates must be ordered in CW direction
%   * The last point is not the first point
% Input:
%   * xc: array that specifies the x-coordinate of the clipping polygon (px1) (must be a rectangle)
%   * yc: array that specifies the y-coordinate of the clipping polygon (px1) (must be a rectangle)
%   * xs: array that specifies the x-coordinate of the subject polygon (px1)
%   * ys: array that specifies the y-coordinate of the subject polygon (px1)
% Output:
%   * xp: array that specifies the x-coordinate of the resulted polygon (px1)
%   * yp: array that specifies the y-coordinate of the resulted polygon (px1)
% Extended:
%   * can include jacobian vector, so xc, yc, xs, ys, and xp, yp will have dimension of (px3)
%   * the rows are: x : [xA, dxA/dxi, dxA/dyi; xB, dxB/dxi, dxB/dyi; ...]
%   * the rows are: y : [yA, dyA/dxi, dyA/dyi; yB, dyB/dxi, dyB/dyi, ...]

function [xp, yp] = clip_polygons_with_rect(xc, yc, xs, ys)
    
    % make the coordinates (x,y) for each polygon
    clippingPolygon = [xc, yc];
    subjectPolygon = [xs, ys];
    
    % get the boundary points
    xcmin = min(xc(:,1)); xcmax = max(xc(:,1));
    ycmin = min(yc(:,1)); ycmax = max(yc(:,1));
    
    outputList = subjectPolygon;
    outputListIdx = size(outputList, 1) + 1;
    
    for (i = [1:size(clippingPolygon,1)])
        % break if there are no point left
        if (length(outputList) == 0) break; end;
        
        % copy the output list and clear it
        inputList = outputList(1:outputListIdx-1,:);
        outputListIdx = 1; % index in outputList array to be filled
        
        S = inputList(end,:);
        for (iE = [1:size(inputList,1)])
            E = inputList(iE,:);
            % SEedge = [S; S-E];
            
            % check if E is inside the clipEdge
            if (isInside(xcmin, xcmax, ycmin, ycmax, E, i))
                
                % check if S is not inside the clipEdge
                if (~isInside(xcmin, xcmax, ycmin, ycmax, S, i))
                    
                    % add the intersection from S to E with the clipEdge
                    outputList(outputListIdx,:) = getIntersection(xcmin, xcmax, ycmin, ycmax, S, S-E, i);
                    outputListIdx = outputListIdx + 1;
                end
                outputList(outputListIdx,:) = E;
                outputListIdx = outputListIdx + 1;
            
            % check if S is inside the clipEdge
            elseif (isInside(xcmin, xcmax, ycmin, ycmax, S, i))
                
                % add the intersection from S to E with the clipEdge
                outputList(outputListIdx,:) = getIntersection(xcmin, xcmax, ycmin, ycmax, S, S-E, i);
                outputListIdx = outputListIdx + 1;
            end
            
            S = E;
        end
    end
    
    if (outputListIdx == 1)
        xp = [];
        yp = [];
    else
        xp = outputList(1:outputListIdx-1,1);
        yp = outputList(1:outputListIdx-1,2:end);
    end
end

function ret = isInside(xmin, xmax, ymin, ymax, point, i)
    if (i == 1)
        ret = point(1) >= xmin;
    elseif (i == 2)
        ret = point(2) <= ymax;
    elseif (i == 3)
        ret = point(1) <= xmax;
    else
        ret = point(2) >= ymin;
    end
end

function intersection = getIntersection(xmin, xmax, ymin, ymax, p0, grad, i)
    if (i == 1)
        intersection = [xmin, p0(2) + grad(2)/grad(1) * (xmin - p0(1))]; % y0+dy/dx*(x-x0);
    elseif (i == 2)
        intersection = [p0(1) + grad(1)/grad(2) * (ymax - p0(2)), ymax]; % x0+dx/dy*(y-y0);
    elseif (i == 3)
        intersection = [xmax, p0(2) + grad(2)/grad(1) * (xmax - p0(1))]; % y0+dy/dx*(x-x0);
    else
        intersection = [p0(1) + grad(1)/grad(2) * (ymin - p0(2)), ymin]; % x0+dx/dy*(y-y0);
    end
end

