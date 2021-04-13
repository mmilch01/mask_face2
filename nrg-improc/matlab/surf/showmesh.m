% name = '~/matlab/meshes/nefertiti.off';
% if(~exist(vertex) || ~exist(faces))
%     [vertex,faces]=read_mesh(name);
%     nbr_landmarks=400;
%     nverts=max(size(vertex));
%     options.edge_color=1;
%     landmark = farthest_point_sampling_mesh( vertex,faces, [], nbr_landmarks, options);
%     [D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmark);
%     [vertex_voronoi, faces_voronoi]=compute_voronoi_triangulation_mesh(Q,vertex,faces);    
%     [normal, normalf]=compute_normal(vertex_voronoi,faces_voronoi);
% end;

function showmesh(faces,vertex,options)
%vertex=fv.vertices;
%faces=fv.faces;
face_color=options.face_color;
if(isempty(face_color)) face_color=[.7 .7 .7]; end;
edge_color=options.edge_color;
if(isempty(edge_color)) edge_color=[.7 .7 .7]; end;

alpha=options.alpha;
if(isempty(alpha)) alpha=1; end;

line_style=options.line_style;
if(isempty(line_style)) line_style='-';end;


patch('vertices',vertex,'faces',faces,...
    'facecolor',face_color,...
    'edgecolor',edge_color,...
    'CData',[NaN NaN NaN],...
    'LineStyle',line_style,...
    'FaceAlpha',alpha,...
    'FaceLighting','none');
%colormap rgb;%gray(256);
%lighting flat;
%camproj('perspective');
axis square;
axis tight;
axis off;
daspect([1 1 1]);
view([-1 0.7 0.5]);
%camlight;
%lighting gouraud;
