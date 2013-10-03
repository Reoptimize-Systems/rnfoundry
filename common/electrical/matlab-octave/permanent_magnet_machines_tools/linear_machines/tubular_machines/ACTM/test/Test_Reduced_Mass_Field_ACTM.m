% Test_Reduced_Mass_Field_ACTM
WmVWp = 0.5;
WpVRm =  0.5;
RoVRm = 1.2;
RshoVRm = 0.2;
Rm = 0.15;

Rs2VRmRshodiff = 0.9; 
Rs1VRmRshodiff = 0.95; 
Ws2VWs = 0.5*0.5;
Ws1VWs = 0.75*0.5;

%RunFEMMSimNew_ACTM(WmVWp, WpVRm, RshoVRm, Rm, mode, Rs2VRmRshodiff, Rs1VRmRshodiff, Ws2VWs, Ws1VWs);

FEMMFL = [0 0 0; 0 0 0];

for mode = 0:1
    
    j = mode + 1;
    
    RunFEMMSimNew_ACTM(WmVWp, WpVRm, RshoVRm, Rm, mode, Rs2VRmRshodiff, Rs1VRmRshodiff, Ws2VWs, Ws1VWs);

    cwVWp = 1/3;

    Ro = RoVRm * Rm;
    Wp = WpVRm * Rm;

    pos = [(0.21667)*Wp,Wp*0.5,Wp-(Wp*0.21667)];

    femmScale = 'm';
    control = 1;

    g = Rm / 33;

    [Ntot(j), dc] = CoilTurns((Ro-g-Rm)*(Wp/3), 0.55, 2.5/1000);

    i = 1;
    %     y = [(pos(i)+0.01*pos(i))-(Wp/6), (pos(i)+0.01*pos(i))-(Wp/6), (pos(i)+0.01*pos(i))+(Wp/6), (pos(i)+0.01*pos(i))+(Wp/6)] * 100;
    %     x = [Rm+g, Ro-0.001, Rm+g, Ro-0.001] * 100;

    y = [pos(i)-(Wp/6), pos(i)-(Wp/6), pos(i)+(Wp/6), pos(i)+(Wp/6)];
    x = [Rm+g, Ro-(0.005*Ro), Rm+g, Ro-(0.005*Ro)];

    mi_drawrectangle(x(1),y(1),x(4),y(4));

    mi_selectsegment((Rm+g),pos(i));
    mi_selectsegment((Ro-(0.005*Ro)),pos(i));
    mi_selectsegment(((Rm+g)+((Ro-(Rm+g))/2)),(pos(i)-(Wp/6)));
    mi_selectsegment(((Rm+g)+((Ro-(Rm+g))/2)),(pos(i)+(Wp/6)));

    mi_setsegmentprop('', 0, 1, 0, 5);

    mi_clearselected;

    nodex = (Rm+g) + ((Ro-(0.005*Ro)) - (Rm+g))/2;
    nodey = pos(i);
    mi_addcircprop('Circuit', 0, 1);
    mi_addblocklabel(nodex, nodey);
    mi_selectlabel(nodex,nodey);
    mi_setblockprop('2.5mm Magnet Wire', 0, (Ro-Rm-g)/40, 'Circuit', 0, 5, Ntot(j));

    mi_clearselected;
    mi_analyse;
    mi_loadsolution;
    vals = mo_getcircuitproperties('Circuit');

    FEMMFL(j,i) = vals(3);

    for i = 2:size(pos,2)

        mi_selectgroup(5);
        mi_movetranslate2(0,(pos(i)-pos(i-1)), 4);
        mi_clearselected;
        mi_analyse;
        mi_loadsolution;

        vals = mo_getcircuitproperties('Circuit');

        FEMMFL(j,i) = vals(3);

        mo_close;

    end

    mi_close;

end
