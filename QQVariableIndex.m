% � fylki S geymum vi� stu�lana vi� breyturnar Q_iQ_j.
% �etta fall gefur v�si � r�tta r�� Q-ana
% T.d fyrir n = 2 �� er r��in � breytunum � S �essi:
% [Q_1Q_1, Q_1Q_2, Q_1Q_3, Q_2Q_2, Q_2Q_3, Q_3Q_3]
% �egar i = 1:m og j = i:m �� gefur QQVariableIndex(n,i,j) v�sinn sem v�sar
% � breytu Q_iQ_j � fylkinu S.
function p = QQVariableIndex(m, i, j)
    p = (i - 1)*m + j - (i - 1)*i / 2;
end


% Fyrir kerfi af st�r� n �� h�fum vi� m = n(n+1)/2 mismunandi Q breytur. S��an
% h�fum vi� samsetningar af Q-unum, samtals m(m+1)/2.

% (i - 1)*m + j segir hve langt vi� erum komin inn � mxm fylki (1-index).
% Viljum upptalningu � efra hornal�nufylkinu svo vi� dr�gum fr� fj�lda staka � ne�ra
% hornl�nufylkinu. �annig �egar vi� erum � l�nu i �� er fj�ldi staka �
% ne�ra hornl�nufylki � ixi fylki jafn i(i-1)/2.