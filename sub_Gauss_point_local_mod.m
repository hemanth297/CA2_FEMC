
function [ GI ]  =  sub_Gauss_point_local_mod(n)
GI = struct;
switch(n)
    case 1
       xi_1d = [ 0 ; 0 ];
wt_1d = [ 1 ; 1 ];

xi_2d = [   [xi_1d(1); xi_1d(2); xi_1d(2); xi_1d(1)]  , ...
            [xi_1d(1); xi_1d(1); xi_1d(2); xi_1d(2)]   ];
wt_2d = [ 1 ; 1 ; 1 ; 1];

    case 2
xi_1d = [ -1 ; 1 ] / sqrt(3);
wt_1d = [ 1 ; 1 ];
% For 2D,
% --------
% | 4  3 |
% | 1  2 |
% --------
xi_2d = [   [xi_1d(1); xi_1d(2); xi_1d(2); xi_1d(1)]  , ...
            [xi_1d(1); xi_1d(1); xi_1d(2); xi_1d(2)]   ];
wt_2d = [ 1 ; 1 ; 1 ; 1];
end

GI.xi_1d = xi_1d;
GI.weight_1d = wt_1d;
GI.xi_2d = xi_2d;
GI.weight_2d = wt_2d;
end




% function [w2D,pt2D]=gaussValues2DQuad(n)    
%     switch (n)
%         case 1
%             w=2; pt=0;
%         case 2
%             w=[1,1]; pt=[-1/sqrt(3), 1/sqrt(3)];
%         case 3
%             w=[5/9, 8/9, 5/9]; pt=[-sqrt(3/5), 0, sqrt(3/5)];
%         case 4
%             w=[(18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36, (18-sqrt(30))/36]; 
%             pt=[-sqrt(3/7-2/7*(sqrt(6/5))), sqrt(3/7-2/7*(sqrt(6/5))),...
%                 -sqrt(3/7+2/7*(sqrt(6/5))), sqrt(3/7+2/7*(sqrt(6/5)))];
%         case 5
%             w=[(322+13*sqrt(70))/900, (322+13*sqrt(70))/900, 128/225,...
%                 (322-13*sqrt(70))/900,(322-13*sqrt(70))/900]; 
%             pt=[-1/3*sqrt(5-2*sqrt(10/7)), 1/3*sqrt(5-2*sqrt(10/7)),0,...
%                 -1/3*sqrt(5+2*sqrt(10/7)), 1/3*sqrt(5+2*sqrt(10/7))];
%         otherwise
%             error('No data are defined for this value');
%     end
%     pt2D=[];
%     w2D=[];
%     for i=1:n
%         for j=1:n
%             pt2D=[pt2D; [pt(i),pt(j)]];
%             w2D=[w2D,w(i)*w(j)];
%         end
%     end    
% end