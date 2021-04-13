function[]=test_rect_mesh(faces, vertices, faces_b, vertices_b, ...
    faces_t,vertices_t)
options.face_color=[.7 .7 .7];
options.edge_color=[1 1 1];
options.alpha=1;
options.line_style='-';
showmesh(faces',vertices',options);

hold on
%show top mesh
options.edge_color=[0 0 1];
options.face_coror=[.3 .3 1];
options.alpha=0;
options.line_style='--';
showmesh(faces_t',vertices_t',options);

%show mod mesh
% options.edge_color=[.7,.7,0];
% options.face_color=[.2,.2,0];
% options.alpha=0;
% options.line_style=':';
% showmesh(faces_mod',vertices_mod',options);

%show bottom mesh
options.edge_color=[0.3 0.7 0.2];
options.alpha=0;
options.line_style='--';
showmesh(faces_b',vertices_b',options);

hold off