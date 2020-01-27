% (From January 21, 2020)
% Analysis of Mouse Behavior 

%% Loading data

mydir = 'C:\ActivePassiveAnalysis'

folders = dir(mydir)
files = arrayfun( @(x) sprintf('%s\\%s\\DisplacementInOdor.mat',x.folder,x.name), folders, 'UniformOutput', false );
files = files(3:end);
output = cellfun( @(x) load( x ), files, 'ErrorHandler', @(x,y) [0], 'UniformOutput', false );

displacements = cellfun( @(x) x.displacement, output , 'ErrorHandler', @(x,y) [0], 'UniformOutput', false );


%% Adding another animal (this is a template, the filename as well as the numbers in lines 32 and 33 may need to be changed for any new animals)

t0 = 60000;
tend = 200000;
newfile = 'C:\ActivePassiveAnalysis\PO372\DisplacementInOdor.mat';
output = load( newfile )
displacements = output.displacement;
newtable = table( {displacements}, {newfile}, [1], [2], ...
    'VariableNames', {'Data','File','CSplus','CSminus'} );

csplus_to_add =  cell2mat( rowfun( @(x,y) x{1}(:,[t0:10:tend],y), newtable(:,{'Data','CSplus'}) ,'OutputFormat', 'cell' ));
csminus_to_add = cell2mat( rowfun( @(x,y) x{1}(:,[t0:10:tend],y), newtable(:,{'Data','CSminus'}) ,'OutputFormat', 'cell' ));

csplus_displacements = [ csplus_displacements; csplus_to_add ];
csminus_displacements = [ csminus_displacements; csminus_to_add ];

include_csminus = [ include_csminus, [191:200] ]; % Including all csminus trials for the new animal
include_csplus = [ include_csplus, [192:200] ]; % Excluding first csplus trial for this animal
%% Creating a table of displacements and the condition (do not run if a new animal is added)

mytable = table( displacements, files, repmat(1,numel(files),1), [ repmat(2,10,1);repmat(3,9,1) ],...
    'VariableNames', {'Data','File','CSplus','CSminus'} );
%
mytable([9:10],:).CSminus = [3;3];
mytable([11:13],:).CSminus = [2;2;2];

%%

t0 = 60000;
tend = 200000;

csplus_displacements =  cell2mat( rowfun( @(x,y) x{1}(:,[t0:10:tend],y), mytable(:,{'Data','CSplus'}) ,'OutputFormat', 'cell' ));
csminus_displacements = cell2mat( rowfun( @(x,y) x{1}(:,[t0:10:tend],y), mytable(:,{'Data','CSminus'}) ,'OutputFormat', 'cell' ));

csplus_locs = arrayfun(@(x) changepts2( csplus_displacements(x,:) ), [1:size(csplus_displacements,1)], 'UniformOutput', false );
csminus_locs = arrayfun(@(x) changepts2( csminus_displacements(x,:) ), [1:size(csminus_displacements,1)], 'UniformOutput', false );

%% Obtaining the change points by manual inspection
% (CS plus)

keypressFxn1 = @(x) fprintf('%i', str2double(x) );

f = figure('WindowKeyPressFcn', @(handle,event) keypressFxn1(event.Key), 'Tag', 'Fig1', 'Position', [0, 350, 560, 420] );

for_figures = table ( [1:size(csplus_displacements,1)]', csplus_locs', csplus_displacements,...
    'VariableNames', {'Index','Locs','Displacement'});

for i=1:size(for_figures,1)
   
    plot( for_figures(i,:).Displacement );
    lines_ = for_figures(i,:).Locs{1};
    
    rowfun( @(x,y) text(y,25,sprintf('%i',x) ), table( [1:numel(lines_)]', lines_ ) );
    arrayfun( @(x) line([x,x],[-20,30],'color','k'), lines_ );
    mychangept(i) = inputdlg('Test');
     
end

% (CS minus)

keypressFxn1 = @(x) fprintf('%i', str2double(x) );

f = figure('WindowKeyPressFcn', @(handle,event) keypressFxn1(event.Key), 'Tag', 'Fig1', 'Position', [0, 350, 560, 420] );

for_figures = table ( [1:size(csminus_displacements,1)]', csminus_locs', csminus_displacements,...
    'VariableNames', {'Index','Locs','Displacement'});

for i=1:size(for_figures,1)
   
    plot( for_figures(i,:).Displacement );
    lines_ = for_figures(i,:).Locs{1};
    
    rowfun( @(x,y) text(y,25,sprintf('%i',x) ), table( [1:numel(lines_)]', lines_ ) );
    arrayfun( @(x) line([x,x],[-20,30],'color','k'), lines_ );
    mychangept_csminus(i) = inputdlg('Test');
     
end

%% Prepare data (only if you found the locations above)

mychangept = cellfun( @(x) str2double(x), mychangept );
mychangept_csminus = cellfun( @(x) str2double(x), mychangept_csminus );

t=table(csplus_locs', mychangept');
for i = 1:size(t,1) if and(t(i,:).Var2>0,~isnan(t(i,:).Var2)); ingress_csplus(i) = t(i,:).Var1{1}(t(i,:).Var2); end; end;

t=table(csminus_locs', mychangept_csminus');
for i = 1:size(t,1) if and(t(i,:).Var2>0,~isnan(t(i,:).Var2)); ingress_csminus(i) = t(i,:).Var1{1}(t(i,:).Var2); end; end;

%% START FROM HERE Once loaded %%

%% Fourier (with boxplot)
%% 

% The variables ingress_csplus and ingress_csminus contain the indices of
% the ingress (obtained from hand selecting within Pia_Jan23_2020.m)

N = 20;
trials = 10;
t0 = 2000;
tend = 3000;

animals_lbl = cell2mat(arrayfun(@(x) repmat(x,trials,1), [1:N], 'UniformOutput', false) ); animals_lbl = animals_lbl(:);
trials_lbl = repmat( [1:trials], 1, N )';

csminus_ = csminus_displacements( include_csminus, : );
csplus_ = csplus_displacements( include_csplus, : );

csplustable = table( csplus_, repmat(t0,size(csplus_,1),1), repmat(tend,size(csplus_,1),1), animals_lbl(include_csplus), trials_lbl(include_csplus) );
csplus_coeffs = rowfun(@(x,y,z) real(myfourier(x,y,z)), csplustable(:,[1:3]) );

csminustable = table( csminus_, repmat(t0,size(csminus_,1),1), repmat(tend,size(csminus_,1),1), animals_lbl(include_csminus), trials_lbl(include_csminus) );
csminus_coeffs = rowfun(@(x,y,z) real(myfourier(x,y,z)), csminustable(:,[1:3]) );

csminus_results = table( csminus_coeffs.Var1, animals_lbl(include_csminus), trials_lbl(include_csminus) );
csplus_results = table( csplus_coeffs.Var1, animals_lbl(include_csplus), trials_lbl(include_csplus) );

figure;
boxplot( csminus_results.Var1, categorical( csminus_results.Var2 ) )
camroll(-90);
ylim([0,0.01])
ylabel('CS minus');

figure;
boxplot( csplus_results.Var1, categorical( csplus_results.Var2 ) )
camroll(-90);
ylim([0,0.01])
ylabel('CS plus');

%% Fourier on all the stuff

t0 = 2000;
tend = 3000;

csplustable = table( csplus_displacements, repmat(t0,size(csplus_displacements,1),1), repmat(tend,size(csplus_displacements,1),1), animals_lbl, trials_lbl );
csminustable = table( csminus_displacements, repmat(t0,size(csminus_displacements,1),1), repmat(tend,size(csminus_displacements,1),1), animals_lbl, trials_lbl );

csplus_coeffs = table2array( rowfun(@(x,y,z) real(myfourier(x,y,z)), csplustable(:,[1:3]) ));
csminus_coeffs = table2array( rowfun(@(x,y,z) real(myfourier(x,y,z)), csminustable(:,[1:3]) ));

csplus_max = max(csplus_displacements,[],2);
csminus_max = max(csminus_displacements,[],2);

csplus_table = table( [1:N*trials]', cell2mat( arrayfun(@(x) repmat(x,trials,1)', [1:N] , 'UniformOutput', false) )',...
    repmat([1:trials]',N,1), [csplus_max], [csplus_coeffs] );
csminus_table = table( [1:N*trials]', cell2mat( arrayfun(@(x) repmat(x,trials,1)', [1:N] , 'UniformOutput', false) )',...
    repmat([1:trials]',N,1), [csminus_max], [csminus_coeffs] );

csplus_table = csplus_table( include_csplus, : );
csminus_table = csminus_table( include_csminus, : );

csplus_table.Properties.VariableNames = {'Index','Mouse','Trial','csplus_max','csplus_coeffs'};
csminus_table.Properties.VariableNames = {'Index','Mouse','Trial','csminus_max','csminus_coeffs'};

%%

ingress_cutoff = 10;
fourier_cutoff = 0.002;

csplus_table_discrete = table( [1:N*trials]', cell2mat( arrayfun(@(x) repmat(x,trials,1)', [1:N] , 'UniformOutput', false) )',...
    repmat([1:trials]',N,1), and(csplus_max>ingress_cutoff,not(csplus_coeffs>fourier_cutoff)), and(not(csplus_max>ingress_cutoff),(csplus_coeffs>fourier_cutoff)), and([csplus_max>ingress_cutoff], [csplus_coeffs>fourier_cutoff]), and(not([csplus_max>ingress_cutoff]),not([csplus_coeffs>fourier_cutoff])) );

csminus_table_discrete = table( [1:N*trials]', cell2mat( arrayfun(@(x) repmat(x,trials,1)', [1:N] , 'UniformOutput', false) )',...
    repmat([1:trials]',N,1), and(csminus_max>ingress_cutoff,not(csminus_coeffs>fourier_cutoff)), and(not(csminus_max>ingress_cutoff),(csminus_coeffs>fourier_cutoff)), and([csminus_max>ingress_cutoff], [csminus_coeffs>fourier_cutoff]), and(not([csminus_max>ingress_cutoff]),not([csminus_coeffs>fourier_cutoff])) );

csplus_table_discrete.Properties.VariableNames = {'Index','Mouse','Trial','Ingress','Tremble','Both','Neither'};
csminus_table_discrete.Properties.VariableNames = {'Index','Mouse','Trial','Ingress','Tremble','Both','Neither'};

summarytbl_plus = arrayfun(@(mouse) sum(table2array(csplus_table_discrete(csplus_table_discrete.Mouse==mouse,[4:7]))), [1:N], 'UniformOutput', false)
summarytbl_plus = reshape( cell2mat( summarytbl_plus )', 4, N )';


summarytbl_minus = arrayfun(@(mouse) sum(table2array(csminus_table_discrete(csminus_table_discrete.Mouse==mouse,[4:7]))), [1:N], 'UniformOutput', false)
summarytbl_minus = reshape( cell2mat( summarytbl_minus )', 4, N )';

figure; subplot(1,2,1); imagesc( summarytbl_plus ); subplot(1,2,2); imagesc( summarytbl_minus )
