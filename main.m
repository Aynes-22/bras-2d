%  Espace de travail d'un bras 2R avec obstacles

clear; close all; clc;

%% Parametres robot et obstacles
L1 = 1;
L2 = 1;

% obstacle 1 : segment vertical
A = [3*L1/2, L1];
seg1 = [A; A + [0, L1/3]];

% obstacle 2 : cercle
C2 = [L1, 0];
R2 = L1/4;

%% Discretisation espace articulaire
N = 400;
th = linspace(0, 2*pi, N);
[Th1, Th2] = meshgrid(th, th);

%% Cinematique directe
P1x = L1 * cos(Th1);
P1y = L1 * sin(Th1);

P2x = P1x + L2 * cos(Th1 + Th2);
P2y = P1y + L2 * sin(Th1 + Th2);

%% Jacobien analytique  det(J) = L1*L2*sin(theta2)
detJ = L1 * L2 * sin(Th2);

%% Detection de collision
% l'origine est fixe (0,0) mais on a besoin de matrices N×N
% pour que les operations soient compatibles avec P1x, P1y
O_x = zeros(N, N);
O_y = zeros(N, N);

% coordonnees des extremites de l'obstacle 1
Pob1x = seg1(1,1);  Pob1y = seg1(1,2);
Qob1x = seg1(2,1);  Qob1y = seg1(2,2);

% segment Origine->P1 contre les deux obstacles
c1 = segCroise(O_x, O_y, P1x, P1y, Pob1x, Pob1y, Qob1x, Qob1y);
c2 = distPointSeg(C2(1), C2(2), O_x, O_y, P1x, P1y) <= R2;

% segment P1->P2 contre les deux obstacles
c3 = segCroise(P1x, P1y, P2x, P2y, Pob1x, Pob1y, Qob1x, Qob1y);
c4 = distPointSeg(C2(1), C2(2), P1x, P1y, P2x, P2y) <= R2;

collision = c1 | c2 | c3 | c4;

%% C-space et workspace
CspaceObs = collision;
ptsTotal = [P2x(:), P2y(:)];
ptsLibre = ptsTotal(~collision(:), :);

%% Resultats
ratioLibre = 100 * size(ptsLibre,1) / size(ptsTotal,1);
disp(['Workspace libre : ' num2str(ratioLibre, '%.1f') ' %']);

%% Figure 1 : C-space
figure('Name','C-space','NumberTitle','off','Color','w');
imagesc(th, th, double(CspaceObs));
colormap([1 1 1; 0.15 0.15 0.15]);
axis xy square; hold on;
[~, hS] = contour(th, th, abs(detJ) < 1e-3, 1, 'r-', 'LineWidth', 1.8);
xlabel('\theta_1 [rad]', 'FontSize', 12);
ylabel('\theta_2 [rad]', 'FontSize', 12);
title('C-space : zones interdites et singularites', 'FontSize', 13);
legend(hS, 'Singularites (det J = 0)', 'Location', 'northeast', 'FontSize', 10);
xticks(0:pi/2:2*pi);
xticklabels({'0','pi/2','pi','3pi/2','2pi'});
yticks(0:pi/2:2*pi);
yticklabels({'0','pi/2','pi','3pi/2','2pi'});
text(0.02, 0.97, sprintf('Zones interdites : %.1f %%', 100-ratioLibre), ...
     'Units','normalized','Color','w','FontSize',10, ...
     'FontWeight','bold','VerticalAlignment','top');

%% Figure 2 : Workspace
figure('Name','Workspace','NumberTitle','off','Color','w');
scatter(ptsTotal(:,1), ptsTotal(:,2), 1, [0.82 0.82 0.82], 'filled'); hold on;
scatter(ptsLibre(:,1), ptsLibre(:,2), 1, [0.1 0.35 0.75], 'filled');
plot(seg1(:,1), seg1(:,2), 'k-', 'LineWidth', 3);
viscircles(C2, R2, 'Color','k', 'LineStyle','--', 'LineWidth', 1.5);
plot(0, 0, 'ks', 'MarkerFaceColor','k', 'MarkerSize', 7);
axis equal; grid on;
xlabel('x [m]', 'FontSize', 12);
ylabel('y [m]', 'FontSize', 12);
title('Workspace total (gris) et libre (bleu)', 'FontSize', 13);
legend({'Total','Libre','Obstacle 1','Obstacle 2','Base'}, ...
       'Location','northwest', 'FontSize', 10);
text(0.02, 0.05, sprintf('Workspace libre : %.1f %%', ratioLibre), ...
     'Units','normalized','FontSize',10,'FontWeight','bold','Color',[0.1 0.35 0.75]);

%% Fonctions

function c = segCroise(Q1x, Q1y, Q2x, Q2y, Px1, Py1, Px2, Py2)
% teste si le segment mobile Q1->Q2 croise le segment fixe P1->P2
% utilise le test d'orientation (produit vectoriel 2D)
    o1 = (Px2-Px1).*(Q1y-Py1) - (Py2-Py1).*(Q1x-Px1);
    o2 = (Px2-Px1).*(Q2y-Py1) - (Py2-Py1).*(Q2x-Px1);
    o3 = (Q2x-Q1x).*(Py1-Q1y) - (Q2y-Q1y).*(Px1-Q1x);
    o4 = (Q2x-Q1x).*(Py2-Q1y) - (Q2y-Q1y).*(Px2-Q1x);
    c  = (o1.*o2 < 0) & (o3.*o4 < 0);
end

function d = distPointSeg(Cx, Cy, Ax, Ay, Bx, By)
% distance du point fixe (Cx,Cy) au segment A->B
% projection orthogonale clampee sur [0,1]
    vx = Bx-Ax;  vy = By-Ay;
    t  = ((Cx-Ax).*vx + (Cy-Ay).*vy) ./ max(vx.^2+vy.^2, 1e-10);
    t  = max(0, min(1, t));
    d  = sqrt((Cx - Ax - t.*vx).^2 + (Cy - Ay - t.*vy).^2);
end
