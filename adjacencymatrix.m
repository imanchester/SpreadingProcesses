function AM=adjacencymatrix(rows,cols)
%ADJACENCYMATRIX 8-node neighbour Adjancency matrix 
% This function creates an adjacency matrix for an landscape of size
% rows by cols, where a node is connected to its 8 direct neighbours (horizontal, vertical and diagonal) 
%   
% Inputs:
% - rows (number of rows)
% - cols (number of columns)
%
% Outputs:
% - AM (adjacency matrix

% Vera Somers, March 2020

%parameter set-up
n=rows*cols;
AM=zeros(n);

%create adjacency matrix

for ii=0:(rows-1)
    for jj=0:(cols-1)
        i=jj*rows + ii;
        if ii > 0
           AM(i,i+1)=1;
           AM(i+1,i)=1;
        end
        if jj > 0
           AM(i+1-rows,i+1)=1;
           AM (i+1,i+1-rows)=1;   
        end
        if ii>0 && jj<(cols-1)
            AM(i+1,i+rows)=1;
           AM(i+rows,i+1)=1;
        end
        if ii>0 && jj>0
           AM(i+1,i-rows)=1;
           AM(i-rows,i+1)=1;
        end
        
    end
end
end