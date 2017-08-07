   function [hasFind,c]= FindOpenTheta(a,b,KI1,KI2)
        [c,fc] = BiSearch(a,b,KI1,KI2);
        g = Fnc2ndOrder(c,KI1,KI2);
        if g > 0 && abs(fc) < 1e-8
            hasFind = 1;
        else
            hasFind = 0;
            c  = -1;
        end
    end
    function f = DeltFnc(theta,KI1v,KI2v)
        f = cosd(theta/2)*(KI1v*sind(theta) + KI2v*(3*cosd(theta)-1))/1e6;
    end
    function g = Fnc2ndOrder(theta,KI1v,KI2v)
        g = 2*KI1v*cosd(theta/2)*(3*cosd(theta)-1) - KI2v*sind(theta/2)*(9*cosd(theta)+5);
        g = g/1e6;
    end
    function [c,fc] = BiSearch(a,b,KI1,KI2)
        dlt = 1e-6;
        dlt2 = 1e-10;
        while abs(a-b) > dlt
            c = (a+b)/2;
            fb = DeltFnc(b,KI1,KI2);
            fc = DeltFnc(c,KI1,KI2);
            if abs(fc) <dlt2
                break;
            else
                if fc*fb < 0
                    a = c;
                else
                    b = c;
                end
            end
        end
    end