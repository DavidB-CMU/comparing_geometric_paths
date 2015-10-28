%
% Example 3:
% Mirrored paths
%

P = [0.0, 1.0, 1.0;
     1.0, 1.0, 1.0;
     2.0, 1.0, 1.0;
     3.0, 1.0, 1.0;
     4.0, 1.0, 1.0;
     5.0, 0.0, 1.0;
     6.0, 1.0, 1.0;
     7.0, 1.0, 1.0;
     8.0, 1.0, 1.0;
     9.0, 1.0, 1.0];

Q = [0.0, 2.0, 1.0;
     1.0, 2.0, 1.0;
     2.0, 2.0, 1.0;
     3.0, 2.0, 1.0;
     4.0, 2.0, 1.0;
     5.0, 3.0, 1.0;
     6.0, 2.0, 1.0;
     7.0, 2.0, 1.0;
     8.0, 2.0, 1.0;
     9.0, 2.0, 1.0];

% Calculate the discrete Frechet distance
[cm, cm_seq] = DiscreteFrechetDistance(P, Q);

% Print the coupling sequence and the distance between points
fprintf('Coupling sequence: \n');
fprintf('  Path P    Path Q    Distance \n');
for i = 1:length(cm_seq)
    u = P(cm_seq(i,1),:);
    v = Q(cm_seq(i,2),:);
    dist = sqrt(sum((u-v).^2));
    fprintf('    %2i        %2i       %f  \n', cm_seq(i,1), cm_seq(i,2), dist);
end

frechet_dist_str = ['Frechet distance = ' num2str(cm)];
fprintf('\n%s \n\n', frechet_dist_str);

% Plot the two paths
figure(1);
hold on;
plot3(Q(:,1),Q(:,2),Q(:,3), '*-r', 'linewidth',2, 'markerfacecolor','r');
plot3(P(:,1),P(:,2),P(:,3), '*-b', 'linewidth',2, 'markerfacecolor','b');
legend('Q','P','location','best');

% Plot the coupling sequence
for i = 1:length(cm_seq)
  line_x_values = [P(cm_seq(i,1),1), Q(cm_seq(i,2),1)];
  line_y_values = [P(cm_seq(i,1),2), Q(cm_seq(i,2),2)];
  line_z_values = [P(cm_seq(i,1),3), Q(cm_seq(i,2),3)];
  line(line_x_values, line_y_values, line_z_values, 'color',[0.5,0.5,0.5]);
end

xlabel('x');
ylabel('y');
zlabel('z', 'rot',0);
axis([min(P(:,1))-1, max(P(:,1))+1, -5, 5, -1, 3])
view([-45, 25])

text(0.0,0.0, frechet_dist_str, 'FontSize',10)
