function Plotting(Sveg,Slike,p0,KK,DD)
%PLOTTING risk map and resource allocation map 
% This function creates two figures; a resource allocation map and risk map
% based on the given spreading process and optimization problem
%
% Inputs:
% - Sveg (vegetation state matrix, eucalyptic forest=1, grassland=2,
% Desert=3, city=4, water/unburnable=5)
% - Slike (likelihood state matrix)
% - p0 (node impact (or priority) vector before resource allocation)
% - KK (resource allocation matrix)
%
% Outputs:
% - figure (1) resource allocation map
% - figure (2) risk map
% - (optional) figure (3) node resource allocation map

% Vera Somers, March 2020


%set up to invoke latex compatibility
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontSize',12)

%initializing parameters
[rows,cols]=size(Slike);
n=rows*cols;
rX=repmat(linspace(0,cols,cols),[rows 1]);
rX=rX(:)';
rY=repmat(linspace(rows,0,rows),[1 cols]);

%call colormap function
cmap=colourmap(3);


%set up resource allocation matrix for plotting
NewK=(KK+KK')-eye(size(KK,1)).*diag(KK);


for k1=1:n
    for k2=1:n
        if abs(NewK(k1,k2))<1E-6
           NewK(k1,k2)=0;
        elseif NewK(k1,k2)>1
           NewK(k1,k2)=1;
        end
    end
end


%resource allocation map
figure(1)
hold on
Sveg2=Sveg/max(max(Sveg));
t=imshow(Sveg2,'InitialMagnification','fit','Colormap',cmap);
alpha(0.2)
G=graph(NewK);
for i=1:size(G.Edges.Weight)
    if G.Edges.Weight(i)>=1
       G.Edges.Weight(i)=0.98;
    end
end
r=plot(G);
r.LineWidth = 2;  
r.EdgeCData=G.Edges.Weight; 
c=colorbar;
c.Label.String = 'Resource allocation';
c.Label.Interpreter='latex';
c.TickLabelInterpreter='latex';
c.Label.FontSize=12;
caxis([0 1])
r.MarkerSize=0.1;   
r.NodeCData =ones(1,1000);
r.XData=rX;
r.YData=fliplr(rY); 
t.XData=rX;
t.YData=fliplr(rY);
axis equal
set(gca,'XTick',[]);
set(gca,'YTick',[]);
%add nodelabels
r.NodeLabel=[];
ax = gca;
ax.Visible = 'off';
C1=t.CData+1;
C2=G.Edges.Weight;
set(t,'CData',C1);
r.EdgeCData=C2;
caxis([0 2])
c.Limits=[0 1];

%paper layout for printing figures
ax.XLim=[0 40];
ax.YLim=[-0.5 25.5];
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left+0.025 bottom ax_width-0.08 ax_height];

fig = gcf;
fig.PaperPositionMode = 'auto';
fig.PaperPosition=[3.0784   10.8151   14.8431    8.0698];
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];




%risk map

cmap2=colourmap(1);
p02 = Slike.*reshape(p0,rows,cols);   
pp = p02/max(max(p02));
figure(2);
imshow(pp,'InitialMagnification','fit','Colormap',cmap2)
c=colorbar;
c.Label.String = 'Risk';
c.Label.Interpreter='latex';
c.TickLabelInterpreter='latex';
c.Label.FontSize=12;

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left-0.041 bottom ax_width ax_height];

fig = gcf;
fig.PaperPositionMode = 'auto';
fig.PaperPosition=[3.0784   10.8151   14.8431    8.0698];
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

%extra plot for delta allocation
if exist('DD','var')==1
    vijS = reshape(DD,rows,cols);
    figure(3);
    imshow(vijS,'InitialMagnification','fit','Colormap',cmap2)
    c=colorbar;
    c.Label.String = 'Resource Allocation';
    c.Label.Interpreter='latex';
    c.TickLabelInterpreter='latex';
    c.Label.FontSize=12;
    
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left-0.041 bottom ax_width ax_height];
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig.PaperPosition=[3.0784   10.8151   14.8431    8.0698];
    %-0.0000    1.5478   14.8167    8.0169
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];

end

end
