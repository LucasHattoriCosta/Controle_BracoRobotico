cont = 1;
lastcont = 1;
Kd = 0.1:0.05:3;%a
Kp = 20:0.2;70%b
Ki = 100:2:1000;%c
for a = 1:length(Kd);
    for b = 1:length(Kp);
        for c = 1:length(Ki);
            Ncc = [Kd(a),Kp(b),Ki(c)];
            Dcc = [1, 0];
            Gcc = tf(Ncc,Dcc);
            GccGH = Gcc*GcGH;
            T1 = feedback(GccGH,1);
            T = 20;
            delT = 0.01;
            t = 0:delT:T;
            u = ones(1,T/delT+1);
            y = lsim(T1,u,t);
            ITAE = 0;
            for i = 1:length(t);
                inc = t(i)*abs(1-y(i));
                ITAE = ITAE + inc;
            end
            Results(cont,:) = [Kd(a),Kp(b),Ki(c),ITAE];
            cont = cont +1;
            lastcont = lastcont + 1;
            porcent = 100*cont/(length(Kd)*length(Kp)*length(Ki));
            if lastcont > (length(Kd)*length(Kp)*length(Ki))/100
                disp(vpa(porcent,3));
                lastcont = 0;
            end
        end
    end
end
ResultsOrd = sortrows(Results,4);
                
          
    