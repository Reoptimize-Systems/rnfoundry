function reportstrs = designreport_RADIAL_SLOTTED(design, simoptions, reportstrs)
% produces a slotted radial flux electrical machine design report in LaTeX
% format, and optionally prduces the pdf using pdflatex
%
% designreport_RADIAL_SLOTTED can also be used to extend a report produced
% by a higher level function by passing in the existing report strings.
%
% Syntax
%
% designreport_RADIAL_SLOTTED(design, simoptions)
% designreport_RADIAL_SLOTTED(design, simoptions, reportstrs)
% reportstrs = designreport_RADIAL_SLOTTED(...)
%
% 
% Input
%
%
% Output
%
%   reportstrs - cell array of strings containing the report
%

    if nargin < 3
        reportstrs = {};
    end
    % make up report on things specific to radial slotted
    
    % generate the dimension table body
    tabledata = { ...
        '$R_{yi}$', 'Inner armature yoke radius. [mm]', design.Ryi * 1000;
        '$R_{yo}$', 'Outer armature yoke radius. [mm]', design.Ryo * 1000;
        '$R_{ym}$', 'Mean armature yoke radius. [mm]', design.Rym * 1000;
        '$R_{ci}$', 'Inner coil radius. [mm]', design.Rci * 1000;
        '$R_{co}$', 'Outer coil radius. [mm]', design.Rco * 1000;
        '$R_{cm}$', 'Mean coil radius. [mm]', design.Rcm * 1000;
        '$R_{tsb}$', 'Tooth shoe base radius [mm]', design.Rtsb * 1000;
    };

    % add the tooth shoe radius at the gap if present
    if isfield (design, 'Rtsg')
        tabledata = [ tabledata; {'$R_{tsg}$', 'Tooth shoe at gap radius. [mm]', design.Rtsg * 1000 } ];
    end

    % tooth surface location
    if strncmpi (design.ArmatureType, 'e', 1)

        tabledata = [ tabledata; { '$R_{ai}$', 'Tooth surface radius. [mm]', design.Rai * 1000 } ];

    elseif strncmpi (design.ArmatureType, 'i', 1)

        tabledata = [ tabledata; {'$R_{ao}$', 'Tooth surface radius. [mm]', design.Rao * 1000 } ];

    end
    
    % add the coil base curve radius if it's present
    if isfield (design, 'Rcb')
        tabledata = [ tabledata; {'$R_{cb}$', 'Coil base curve beginning radius. [mm]', design.Rcb * 1000 } ];
    end
    
    tabledata = [ tabledata; { ...
        '$y_d$', 'Coil pitch in slots. (No. Of Slots)', design.yd;
        '$\theta_{c}$', 'Coil pitch. [rad]', design.yd*design.thetas;
        '$\tau_{c}$', 'Coil pitch at mean coil radius. [mm]', design.yd*design.thetas*design.Rcm * 1000;
        '$\theta_{s}$', 'Slot centers pitch angle. [rad]', design.thetas;
        '$\theta_{sg}$', 'Slot opening angle. [rad]', design.thetasg;
        '$\tau_{sg}$', 'Slot opening pitch. [mm]', design.thetasg * design.Rtsg * 1000;
        '$t_c$', 'Coil slot height in the radial direction. [mm]', design.tc(1) * 1000;
        '$t_{cb}$', 'Coil slot base height in the radial direction. [mm]', design.tc(2) * 1000;
        '$t_y$', 'Coil yoke thickness in the radial direction. [mm]', design.ty * 1000;
        '$t_{sb}$', 'Coil shoe thickness where shoe meets tooth. [mm]', design.tsb * 1000;
        '$t_{sg}$', 'Coil shoe thickness at coil slot opening. [mm]', design.tsg * 1000;
        '$\theta_{cg}$', 'Internal slot pitch angle at gap end. [rad]', design.thetacg;
        '$\tau_{cg}$', 'Internal slot pitch at gap end. [mm]', design.thetacg.*design.Rtsb * 1000;
        '$\theta_{cy}$', 'Internal slot pitch angle at yoke end (coil base). [rad]', design.thetacy;
        '$\tau_{cg}$', 'Internal slot pitch at gap end. [mm]', design.thetacg.*design.Rtsb * 1000;
    } ];

    % radial flux armature dimensions
    if strncmpi (design.ArmatureType, 'e', 1)

        tabledata = [ tabledata; { '$\tau_{cy}$', 'Internal slot pitch at yoke end (coil base). [mm]', design.thetacy.*design.Ryi * 1000 } ];

    elseif strncmpi (design.ArmatureType, 'i', 1)

        tabledata = [ tabledata; {'$\tau_{cy}$', 'Internal slot pitch at yoke end (coil base). [mm]', design.thetacy.*design.Ryo * 1000 } ];

    end
    
    % generate the LaTex table of the outputs
    colheadings = {};
    rowheadings = {};
    colsep = ' & '; 
    rowending = ' \\';
    fms = '.3f';

    fname = tempname;
    
    fid = fopen(fname, 'w');
    
    if fid == -1
        error('Temporary file could not be opened.');
    end
    
    displaytable(tabledata, colheadings, [30,60,20], fms, rowheadings, fid, colsep, rowending);
    
    % close the file
    fclose(fid);
    
    radarmdimstablestrs = txtfile2cell(fname);
    
    h = figure;
    plot (design.MagFEASimPositions, design.ArmatureToothFluxDensity);
    set (gcf, 'color', 'w');
    xlabel ('Position (Fraction of Pole)');
    ylabel ('Flux Density Magnitude (T)');
    matlab2tikz ([fname, '.tikz'], 'width', '8cm', 'showInfo', false);
    close (h)
    
    fluxdensitystrs = txtfile2cell([fname, '.tikz']);
    
    h = figure;
    plot (design.MagFEASimMaterials.ArmatureYoke.BHPoints(:,2), ...
          design.MagFEASimMaterials.ArmatureYoke.BHPoints(:,1));
	hold on;
    xlim = get (gca, 'XLim');
    plot ([design.MagFEASimMaterials.ArmatureYoke.BHPoints([1,end],2)], ...
          [design.ArmatureToothFluxDensityPeak, design.ArmatureToothFluxDensityPeak], ...
           'k:' );
    hold off;
%     ind = find ( design.MagFEASimMaterials.ArmatureYoke.BHPoints(:,1) > 2.0 );
%     if isempty (ind), ind = size (design.MagFEASimMaterials.ArmatureYoke.BHPoints, 1); end
    set (gca, 'XLim', [xlim(1), 16000]);
    set (gcf, 'color', 'w');
    xlabel ('H');
    ylabel ('B');
    matlab2tikz ([fname, '.tikz'], 'width', '8cm', 'mathmode', false, 'showInfo', false);
    close (h)
    
    bhcurvestrs = txtfile2cell([fname, '.tikz']);
    
    % now make the table of values with simulation data relevant to only
    % sloted radial flux machines
    
radslottedreportstrs = [{...
'% Beginning section of report generated by designreport_RADIAL_SLOTTED.m';
'\section{Armature Flux Information}';
'The flux linkage has been sampled in a line going down the centre of a tooth, ';
'and the absolute value of the magnitude of the flux density obtained. This ';
'has been repeated for every position used to obtain the other machine data.';
'The resulting plot of these values at each position is shown in Figure \ref{fig:peakarmflux}.';
sprintf('The peak of these values is %3.2f~T', design.ArmatureToothFluxDensityPeak);
'\begin{figure}[h]';
'\centering';
};
fluxdensitystrs;
{ ...
'\caption{Magnitude of flux density along tooth center line at multiple rotor positions.}';
'\label{fig:peakarmflux}';
'\end{figure}';
'';
'The B-H curve for the material used in the armature is shown in Figure \ref{fig:armbhcurve}.';
'\begin{figure}[h]';
'\centering';
};
bhcurvestrs;
{ ...
sprintf('\\caption{B-H used curve for armature iron material (%s). Dotted line shows peak magnutude of B in tooth.}', design.MagFEASimMaterials.ArmatureYoke.Name);
'\label{fig:armbhcurve}';
'\end{figure}';
'\FloatBarrier';
'\section{Dimensions of the design}'; 
'The dimensions of the slotted radial flux design.';
'';
'\begin{table}[htb]';
'\centerline{';
'\begin{tabular}{lll}'
'\toprule';
'Dimension & Description & Value \\';
'\midrule';
}; ...
radarmdimstablestrs;
{...
'\bottomrule';
'\end{tabular}';
'}';
'\caption{Armature dimensions.}';
'\end{table}';
'\begin{table}[htb]';
'\centerline{';
'\begin{tabular}{ll}'
'\toprule';
'Component & Material \\';
'\midrule';
['Armature Iron &', design.MagFEASimMaterials.ArmatureYoke.Name, '\\'];
'\bottomrule';
'\end{tabular}';
'}';
'\caption{Armature materials.}';
'\end{table}';
'% FloatBarrier requires the placeins package';
'\FloatBarrier{}';
'';
}; ...
];
    
    % append stuff common to all maradial machines
    radslottedreportstrs = designreport_RADIAL(design, simoptions, radslottedreportstrs);

        % append the report strings
    reportstrs = [ reportstrs;
                   radslottedreportstrs ];
               
end
