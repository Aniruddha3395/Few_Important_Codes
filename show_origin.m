function show_origin()

hold on;
scatter3d(part_ptcloud,'.b');
hold on;
scatter3(m0(1),m0(2),m0(3),100,'d','filled','b');
hold on;
quiver3(m0(1),m0(2),m0(3),1,0,0,100,'r');hold on;quiver3(m0(1),m0(2),m0(3),0,1,0,100,'g');hold on;quiver3(m0(1),m0(2),m0(3),0,0,1,100,'b');
hold on;
scatter3d(scan_ptcloud,'.r');
hold on;
end