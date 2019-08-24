% INITIAL CONDITIONS FOR STABILITY CHECK OF DROPLET ASSEMBLIES
% And Random configurations

% These configurations are symmetric versions of the shapes reported in
% Shen et al Adv. Sc. 2016


switch intl
    case 2.1
        x0=[0;2.2];
        y0=[0;0];
    case 3.1
        x0=[0;2.2;4.4];
        y0=[0;0;0];
    case 3.2
        x0=[0;2.2;1.1];
        y0=[0;0;1.9053];
    case 4.1
        x0=[0;2.2;4.4;6.6];
        y0=[0;0;0;0];
    case 4.2
        x0=[0;-2.2;-2.2-round(sqrt(3),3)*1.1;-2.2-round(sqrt(3),3)*1.1];
        y0=[0;0;1.1;-1.1];
    case 4.3
        x0=[0;sqrt(3)*1.1;sqrt(3)*1.1;2*sqrt(3)*1.1];
        y0=[0;1.1;-1.1;0];
    case 5.1
        x0=[0;2.2;4.4;6.6;8.8];
        y0=[0;0;0;0;0];
    case 5.2
        x0=[0;-2.2;-4.4;-4.4-sqrt(3)*1.1;-4.4-sqrt(3)*1.1];
        y0=[0;0;0;1.1;-1.1];
    case 5.3
        x0=[0;-2.2;-2.2-sqrt(3)*1.1;-2.2-sqrt(3)*1.1;-2.2-2*sqrt(3)*1.1];
        y0=[0;0;1.1;-1.1;0];
    case 5.4
        x0=[0;2.2;4.4;1.1;3.3];
        y0=[0;0;0;sqrt(3)*1.1;sqrt(3)*1.1];
    case 5.5
        x0=[0;2.2;4.4;1.1;1.1];
        y0=[0;0;0;sqrt(3)*1.1;-sqrt(3)*1.1];
    case 6.1
        x0=[0;2.2;4.4;6.6;8.8;11];
        y0=[0;0;0;0;0;0];
    case 6.2
        x0=[4.4-sqrt(3)*1.1;4.4-sqrt(3)*1.1;4.4;6.6;8.8;11];
        y0=[1.1;-1.1;0;0;0;0];
    case 6.3
        x0=[4.4-sqrt(3)*1.1;4.4-sqrt(3)*1.1;4.4;6.6;6.6+sqrt(3)*1.1;6.6+sqrt(3)*1.1];
        y0=[1.1;-1.1;0;0;1.1;-1.1];
    case 6.4
        x0=[2.2;0;-2.2;-2.2-sqrt(3)*1.1;-2.2-sqrt(3)*1.1;-2.2-2*sqrt(3)*1.1];
        y0=[0;0;0;1.1;-1.1;0];
    case 6.5
        x0=[0;2.2;2.2;2.2;2.2+sqrt(3)*1.1;2.2+sqrt(3)*1.1];
        y0=[0;0;-2.2;2.2;1.1;-1.1];
    case 6.6
        x0=[0;2.2;4.4;1.1;3.3;2.2];
        y0=[0;0;0;sqrt(3)*1.1;sqrt(3)*1.1;2*sqrt(3)*1.1];
    case 6.7
        x0=[0;2.2;4.4;1.1;3.3;5.5];
        y0=[0;0;0;sqrt(3)*1.1;sqrt(3)*1.1;sqrt(3)*1.1];
    case 6.8
        x0=[0;1.1;1.1;2.2;3.3;3.3];
        y0=[0;sqrt(3)*1.1;-sqrt(3)*1.1;0;-sqrt(3)*1.1;sqrt(3)*1.1];
    case 'rand1'  % RANDOM INITIAL CONFIGURATION
        dd=-1;
        ddrop=1e3;
        p=0;
        while dd<=0
            x0=0:2.2:2.2*(n-1);
            y0=1.8*rand(n,1);
            for jj=1:n
                for kk=1:n
                    if kk~=jj
                        ddrop=min(ddrop,norm([x0(jj)-x0(kk) y0(jj)-y0(kk)])-2);
                    end
                end
            end
            dd=max(dd,ddrop);
            p=p+1;
        end
        disp(p)
end