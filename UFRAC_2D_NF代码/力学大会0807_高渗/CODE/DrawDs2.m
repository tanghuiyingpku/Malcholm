function DrawDs2(Ds,nAllEle,AllEle)
global isMechActive_global;
figure;
Ds = Ds*5000;
for i = 1 : nAllEle
    plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'b','Linewidth',0.5);
    hold on;
%     plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'r.','Markersize',4);
   
    if (abs(AllEle(i,10) - 3) < 1e-6) ||  abs(AllEle(i,10) - 1) < 1e-6 || abs(AllEle(i,10) - 6) < 1e-6
        text(AllEle(i,8), AllEle(i,9),num2str(i));
        if isMechActive_global(i) == 0
            plot(AllEle(i,8) ,AllEle(i,9),'m.','Markersize',11);
        end
        if isMechActive_global(i) == -2
            plot(AllEle(i,8) ,AllEle(i,9),'k.','Markersize',11);
        end
        if isMechActive_global(i) == -3
            plot(AllEle(i,8) ,AllEle(i,9),'R.','Markersize',11);
        end
        if isMechActive_global(i) == 1
            plot(AllEle(i,8) ,AllEle(i,9),'b.','Markersize',11);
        end
        
        if Ds(i) < -1e-10
          %  plot([AllEle(i,8)+abs(Ds(i))/2*AllEle(i,5),AllEle(i,8)-abs(Ds(i))/2*AllEle(i,5)],[AllEle(i,9)-abs(Ds(i))/2*AllEle(i,6),AllEle(i,9)+abs(Ds(i))/2*AllEle(i,6)],'b','LineWidth',0.5);
            
            plot([AllEle(i,8)+abs(Ds(i))/2*AllEle(i,5),AllEle(i,8)-abs(Ds(i))/2*AllEle(i,5)],[AllEle(i,9)-abs(Ds(i))/2*AllEle(i,6),AllEle(i,9)+abs(Ds(i))/2*AllEle(i,6)],'b','Markersize',4,'LineWidth',1.0);
        else
           % plot([AllEle(i,8)+abs(Ds(i))/2*AllEle(i,5),AllEle(i,8)-abs(Ds(i))/2*AllEle(i,5)],[AllEle(i,9)-abs(Ds(i))/2*AllEle(i,6),AllEle(i,9)+abs(Ds(i))/2*AllEle(i,6)],'r','LineWidth',0.5);
            
            plot([AllEle(i,8)+abs(Ds(i))/2*AllEle(i,5),AllEle(i,8)-abs(Ds(i))/2*AllEle(i,5)],[AllEle(i,9)-abs(Ds(i))/2*AllEle(i,6),AllEle(i,9)+abs(Ds(i))/2*AllEle(i,6)],'r','Markersize',4,'LineWidth',1.0);
        end
        % plot([AllEle(i,1) AllEle(i,3)],[AllEle(i,2) AllEle(i,4)],'k^');
    end
end

axis equal;
% axis([-0.6 0.6 2.5 4]);

end