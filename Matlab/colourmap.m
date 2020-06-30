function cmap=colourmap(map)
%PLOTTING three colour maps; redblue shaded, vegetation and both 
% This function creates a colourmap according to input 'map':
% map=1, shades of blue going to red colourmap
% map=2, vegetation colourmap
% map=3, both (for overlaying plots)
%
% Inputs:
% - map (number indicates which map is needed)
%
% Outputs:
% - cmap ([r g b] values of specified colourmap (120 by 3))

% Vera Somers, March 2020

%set-up blue red map 1

n = 122*0.5;
r = (0:n-1)'/max(n-1,1);
g = r;
r = [r; ones(n,1)];
g = [g; flipud(g)];
b = flipud(r);
map1 = [r g b];
map1(all(map1==1,2),:)=[];


%set-up map 2 landscape

if map==2
    x=1;
else
    x=20;
end
     

map2=[repmat([0 0 0],x,1)
    repmat([0 0.5 0],x,1)
    repmat([0 1 0],x,1)
    repmat([1 1 0],x,1)
    repmat([0.5 0.5 0.5],x,1)
    repmat([0 0 1],x,1)];

%create correct output
if map==1
   cmap=map1;
elseif map==2
   cmap=map2;
else
   cmap=[map1;map2];
end


end
