function h = contourplot(indepvar1, indepvar2, depvar)

    [XI,YI] = meshgrid(min(indepvar1):max(indepvar1)/500:max(indepvar1),min(indepvar2):max(indepvar2)/500:max(indepvar2));

    ZI = griddata(indepvar1,indepvar2,depvar,XI,YI);

    [C,h] = contour(XI,YI,ZI);

    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2, 'LineWidth', 2);

    colormap winter %hot

    clabel(C,h,'FontSize',24);
    
end