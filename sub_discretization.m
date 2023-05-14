function [Mesh] = sub_discretization(Model)

Mesh = struct;
% Nodes are uniformlyl populated.% The nodal coordinates are stored in x_node.% x_node : num_nodes - by - 2
x_min = Model.domain.minmax(1,1);
x_max = Model.domain.minmax(2,1);
y_min = Model.domain.minmax(1,2);
y_max = Model.domain.minmax(2,2);

num_elem_x  =  Model.mesh.interval(1);
num_elem_y  =  Model.mesh.interval(2);

x_line     = linspace ( x_min , x_max , num_elem_x+1 )';
y_line     = linspace ( y_min , y_max , num_elem_y+1 )';
one_vector = ones ( size(x_line) );

x_node = [];
for j = 1 : num_elem_y+1
    x_node = [ x_node ; [x_line, y_line(j)*one_vector] ];
end

Mesh.x_node = x_node ;

%% Connectivity
num_elem_x  =  Model.mesh.interval(1) ;
num_elem_y  =  Model.mesh.interval(2) ;
num_element =  num_elem_x * num_elem_y;  % the number of finite elements in the domain

NEN = 4;  % num of nodes per element
connectivity = zeros ( NEN , num_element );

local_to_global_block = [ 0 ; 1 ; 1+(num_elem_x+1) ; 0+(num_elem_x+1) ] ;

x_max = Model.domain.minmax(2,1);
y_max = Model.domain.minmax(2,2);

element_index = 0;
for node_index = 1 : length(x_node)  % loop over nodes
    x = x_node(node_index,:);
    if x(1) < x_max-1e-8  &&  x(2) < y_max-1e-8  % skip the nodes on the right and upper boundaries.
        element_index = element_index + 1;
        % The (node-index)-th node is the lower-left-corner node of the
        % currunt element.
        connectivity( : , element_index )  =  node_index + local_to_global_block;
    end
end

Mesh.connectivity = connectivity ;
%plot_mesh ( Mesh.x_node , Mesh.connectivity ); title('Undeformed Mesh'); % plot the mesh to check if meshing was correctly done.

end






