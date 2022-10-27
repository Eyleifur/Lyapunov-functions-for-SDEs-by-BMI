% Í fylki S geymum við stuðlana við breyturnar Q_iQ_j.
% Þetta fall gefur vísi á rétta röð Q-ana
% T.d fyrir n = 2 þá er röðin á breytunum í S þessi:
% [Q_1Q_1, Q_1Q_2, Q_1Q_3, Q_2Q_2, Q_2Q_3, Q_3Q_3]
% Þegar i = 1:m og j = i:m þá gefur QQVariableIndex(n,i,j) vísinn sem vísar
% á breytu Q_iQ_j í fylkinu S.
function p = QQVariableIndex(m, i, j)
    p = (i - 1)*m + j - (i - 1)*i / 2;
end


% Fyrir kerfi af stærð n þá höfum við m = n(n+1)/2 mismunandi Q breytur. Síðan
% höfum við samsetningar af Q-unum, samtals m(m+1)/2.

% (i - 1)*m + j segir hve langt við erum komin inn í mxm fylki (1-index).
% Viljum upptalningu á efra hornalínufylkinu svo við drögum frá fjölda staka í neðra
% hornlínufylkinu. Þannig þegar við erum í línu i þá er fjöldi staka í
% neðra hornlínufylki í ixi fylki jafn i(i-1)/2.