function Cp2 = UpdateCp_x(Cp,nAct,nprop,height,hbanking,hbanking0,ConnList,Vprop_x,Schindex,e,e0,dt,Lv)
global well nwell;
global  Index;
global Fluid
diam = zeros(nprop,1);
for ii = 1 : nprop
    diam(ii) = Fluid.prop{ii}.diam;
end
Vinit = 0;
Vupdate = 0;
Cp2 = Cp;
MatCp = zeros(nAct,nAct);
RhsCp =zeros(nAct,1);
for j = 1 : nprop
    MatCp = MatCp *0;
    RhsCp = RhsCp *0;
    for i = 1 : nAct
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
            if Conj < 1e-6
                continue;
            end
            num = Index(Conj);
            if num < 1e-6
                continue;
            end
            ei = 0.5*(e(i) + e(num));
            if Vprop_x(i+(j-1)*nAct,kk) > -1e-6 %&& Cp(num,j) > 1e-8
                %           From outside to the current grid
                k=1;
                
                if Cp(num,j) > 1e-5
                    if diam(j)*2 > e(i)
                        k=0;
                    else
                        if diam(j)*4 > e(i)
                            k = (e(i)-diam(j)*2)/2/diam(j);
                        else
                            k=1;
                        end
                    end
                    MatCp(i,num) = MatCp(i,num) - k*Vprop_x(i+(j-1)*nAct,kk) *ei*flow_height;
                end
            else
                k=1;
                
                if Cp(i,j) > 1e-5
                    if diam(j)*2 >  e(num)
                        k=0;
                    else
                        if diam(j)*4 > e(num)
                            k = (e(num)-diam(j)*2)/2/diam(j);
                        else
                            k=1;
                        end
                    end
                    MatCp(i,i) = MatCp(i,i) - k*Vprop_x(i+(j-1)*nAct,kk) *ei*flow_height;
                end
            end
        end
        %  end
        Vinit = Vinit + Cp(i,j)*Lv(i)*flow_height0*e0(i);
        MatCp(i,i) = MatCp(i,i)+1/dt*e(i)*flow_height*Lv(i);
        RhsCp(i) = RhsCp(i)+1/dt*e0(i)*flow_height0*Lv(i)*Cp(i,j);
    end
    %Source term from well
    inj = 0;
    for ii = 1 : nwell
        for jj = 1 :well{ii}.Sch(Schindex(ii)).nPf
            iele =  well{ii}.Sch(Schindex(ii)).Pf(jj);
            num = Index(well{ii}.Perfindex(iele,:));
            if min(num) > 1e-6
                inj = inj + well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1)*well{ii}.Sch(Schindex(ii)).PfCp(jj,j);
                inj = inj + well{ii}.Sch(Schindex(ii)).Pf_Q(jj,2)*well{ii}.Sch(Schindex(ii)).PfCp(jj,j);
                RhsCp(num(1)) = RhsCp(num(1)) + well{ii}.Sch(Schindex(ii)).Pf_Q(jj,1)*well{ii}.Sch(Schindex(ii)).PfCp(jj,j);
                RhsCp(num(2)) = RhsCp(num(2)) + well{ii}.Sch(Schindex(ii)).Pf_Q(jj,2)*well{ii}.Sch(Schindex(ii)).PfCp(jj,j);
            end
        end
    end
    
    %     if det(MatCp) < 1e-16
    %         Cp2(:,j) = Cp2(:,j)*0;
    %     else
    %         Cp2(:,j) = MatCp\RhsCp;
    %     end
    Cp2(:,j) = MatCp\RhsCp;
    for i  =  1 : nAct
        flow_height = height - sum(hbanking(i,:));
        Vupdate= Vupdate + Cp2(i,j)*Lv(i)*flow_height*e(i);
    end
    error = (Vupdate - Vinit - inj*dt);
end
end