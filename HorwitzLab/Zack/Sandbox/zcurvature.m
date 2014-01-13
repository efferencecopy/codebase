function [Kout, K, V] = zcurvature(F)
[Fx Fy Fz] = gradient(F);
% dude = abs([Fx(:) Fy(:) Fz(:)])';
% woah = bsxfun(@le, dude, eps*ones(size(dude)));
% [~,idx,~] = unique(woah','rows');
% idx = idx(all(woah(:,idx)));
% Fx(idx) = 1e-12; Fy(idx) = 1e-12; Fz(idx) = 1e-12;

[Fxx Fxy Fxz] = gradient(Fx);
[Fyx Fyy Fyz] = gradient(Fy);
[Fzx Fzy Fzz] = gradient(Fz);

Fx = Fx(:); Fy = Fy(:); Fz = Fz(:);
Fxx = Fxx(:); Fxy = Fxy(:); Fxz = Fxz(:);
Fyx = Fyx(:); Fyy = Fyy(:); Fyz = Fyz(:);
Fzx = Fzx(:); Fzy = Fzy(:); Fzz = Fzz(:);

K = zeros(length(Fx),2); % principal curvatures
V = zeros(length(Fx)*2,3); % corresponding eigenvectors

% Ktest = zeros(size(Fx));
% for i = 1:length(Fx)
%     delF = [Fx(i) Fy(i) Fz(i)];
%     H = [Fxx(i) Fxy(i) Fxz(i);Fyx(i) Fyy(i) Fyz(i);Fzx(i) Fzy(i) Fzz(i)];
%     Ktest(i) = -det([H delF'; delF 0]) / norm(delF)^4;
% end
% Ktest = reshape(Ktest,size(F));

% for i = 1:length(Fx)
%     s = 1; % Positive curvature on positive side of F, otherwise set s = -1
%     N = [Fx(i) Fy(i) Fz(i)]; % The gradient of f at (x0,y0,z0)
%     M = null(N); % Vectors in M and the gradient N are mutually orthogonal
%     H = [Fxx(i) Fxy(i) Fxz(i);Fyx(i) Fyy(i) Fyz(i);Fzx(i) Fzy(i) Fzz(i)]; % The Hessian of f at (x0,y0,z0)
%     [U,D] = eig(M.'*H*M); % These eigenvalues are proportional to curvatures
%     K(i,:) = -s*diag(D).'/norm(N); % Convert to the true curvatures
%     idx = 2*(i-1)+1;
%     V(idx:idx+1,:) = (M*U)'; % Obtain the corresponding eigenvectors in curvature directions
% end

Kout = ((Fz.*(Fxx.*Fz - 2*Fx.*Fxz) + Fzz.*Fx.^2).*(Fz.*(Fyy.*Fz - 2*Fy.*Fyz) + Fzz.*Fy.^2) - (Fz.*(-Fx.*Fyz + Fxy.*Fz - Fxz.*Fy) + Fx.*Fy.*Fzz).^2) ./ (Fz.^2.*(Fx.^2 + Fy.^2 + Fz.^2).^2);
Kout = reshape(Kout,size(F));