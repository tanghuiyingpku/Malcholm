function Xf2 = UpdataXf(height,hbanking,hbanking0,densfv,densf0v,ConnList,Schindex,e,e0,nAct,nfluid,Vfl,Xf,Cp,dt,Lv)
global well nwell;
global  Index Mat;
Vinit = 0;
Vupdate = 0;
Xf2 = Xf;
MatXf = zeros(nAct,nAct);
RhsXf = zeros(nAct,1);
for j = 1 : nfluid
    MatXf = MatXf *0;
    RhsXf = RhsXf *0;
    for i = 1 : nAct
        sumcp = sum(Cp(i,:));
        flow_height = height - sum(hbanking(i,:));
        flow_height0 = height - sum(hbanking0(i,:));
        nConn = ConnList(i,2);
       % if Cp2(i,j) > 1e-8
            for kk = 1 : nConn
                if ConnList(i,kk+2) > -0.1
                    Conj =  ConnList(i,kk+2);
                else
                    continue;
                end
                num = Index(Conj);
                if num < 1e-6
                    continue;
                end
                ei = 0.5*(e(i) + e(num));
                if Vfl(i,kk) > -1e-6 %&& Cp(num,j) > 1e-8
                    % From outside to the current grid
                    MatXf(i,num) = MatXf(i,num) - densfv(num)*Vfl(i,kk) *ei*flow_height*(1-sum(Cp(num,:)));
                else
                    MatXf(i,i) = MatXf(i,i) - densfv(i)*Vfl(i,kk) *ei*flow_height*(1-sumcp);
                end
            end
      %  end
        Vinit = Vinit + Xf(i,j)*Lv(i)*flow_height0*e0(i)*(1-sumcp)*densf0v(i);
        MatXf(i,i) = MatXf(i,i)+1/dt*e(i)*flow_height*Lv(i)*(1-sumcp)*densfv(i);
        RhsXf(i) = RhsXf(i)+1/dt*e0(i)*flow_height0*Lv(i)*Xf(i,j)*(1-sumcp)*densf0v(i);
    end
    %Source term from well
    inj = 0;
    dens0 = CalcFLDens(j,Mat.Pp);
    for ii = 1 : nwell
        for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
            iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
            num = Index(well{ii}.Perfindex(iele,:));
            if min(num) > 1e-6
                Cp1 = well{ii}.Sch(Schindex).PfCp(jj,:);
                Cp2 = well{ii}.Sch(Schindex).PfCp(jj,:);
                inj = inj + well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1)*well{ii}.Sch(Schindex).PfXf(jj,j);
                inj = inj + well{ii}.Sch(Schindex(ii)).Pf_Q(jj,2)*well{ii}.Sch(Schindex).PfXf(jj,j);
                RhsXf(num(1)) = RhsXf(num(1)) + well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1)*well{ii}.Sch(Schindex).PfXf(jj,j)*dens0*(1-sum(Cp1));
                RhsXf(num(2)) = RhsXf(num(2)) + well{ii}.Sch(Schindex(ii)).Pf_Q(jj,2)*well{ii}.Sch(Schindex).PfXf(jj,j)*dens0*(1-sum(Cp2));
            end
        end
    end
   
%     if det(MatXf) < 1e-16
%         Cp2(:,j) = Cp2(:,j)*0;
%     else
%         Cp2(:,j) = MatXf\RhsXf;
%     end
    Xf2(:,j) = MatXf\RhsXf;
    for i  =  1 : nAct
        sumcp = sum(Cp(i,:));
        flow_height = height - sum(hbanking(i,:));
        Vupdate= Vupdate + Xf2(i,j)*Lv(i)*flow_height*e(i)*(1-sumcp)*densfv(i);
    end
    error = (Vupdate - Vinit - inj*dt);  
end

end