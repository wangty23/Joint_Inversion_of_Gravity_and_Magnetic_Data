function showboxf(A)
[nqy,~]=size(A);
hold on
for i=1:nqy
    plot([A(i,1),A(i,1)],[A(i,3),A(i,4)],'k','LineWidth',1);
    plot([A(i,2),A(i,2)],[A(i,3),A(i,4)],'k','LineWidth',1);
    plot([A(i,1),A(i,2)],[A(i,3),A(i,3)],'k','LineWidth',1);
    plot([A(i,1),A(i,2)],[A(i,4),A(i,4)],'k','LineWidth',1);
end
end