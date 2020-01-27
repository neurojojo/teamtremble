% January 21, 2020 
% Analysis of Mouse Behavior 
%
% Fourier decomposition for pre-ingress period

mydir = 'C:\ActivePassiveAnalysis'

folders = dir(mydir)
files = arrayfun( @(x) sprintf('%s\\%s\\DisplacementInOdor.mat',x.folder,x.name), folders, 'UniformOutput', false );
files = files(3:end);
output = cellfun( @(x) load( x ), files, 'ErrorHandler', @(x,y) [0], 'UniformOutput', false );

displacements = cellfun( @(x) x.displacement, output , 'ErrorHandler', @(x,y) [0], 'UniformOutput', false );

%%

mytable = table( displacements, files, repmat(1,numel(files),1), [ repmat(2,10,1);repmat(3,9,1) ], 'VariableNames', {'Data','File','CSplus','CSminus'} );

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

%% CS plus

keypressFxn1 = @(x) fprintf('%i', str2double(x) );

f = figure('WindowKeyPressFcn', @(handle,event) keypressFxn1(event.Key), 'Tag', 'Fig1', 'Position', [0, 350, 560, 420] );

for_figures = table ( [1:size(csplus_displacements,1)]', csplus_locs', csplus_displacements,...
    'VariableNames', {'Index','Locs','Displacement'});

for i=1:190%41
   
    plot( for_figures(i,:).Displacement );
    lines_ = for_figures(i,:).Locs{1};
    
    rowfun( @(x,y) text(y,25,sprintf('%i',x) ), table( [1:numel(lines_)]', lines_ ) );
    arrayfun( @(x) line([x,x],[-20,30],'color','k'), lines_ );
    mychangept(i) = inputdlg('Test');
     
end

%% CS minus

keypressFxn1 = @(x) fprintf('%i', str2double(x) );

f = figure('WindowKeyPressFcn', @(handle,event) keypressFxn1(event.Key), 'Tag', 'Fig1', 'Position', [0, 350, 560, 420] );

for_figures = table ( [1:size(csminus_displacements,1)]', csminus_locs', csminus_displacements,...
    'VariableNames', {'Index','Locs','Displacement'});

for i=1:190%41
   
    plot( for_figures(i,:).Displacement );
    lines_ = for_figures(i,:).Locs{1};
    
    rowfun( @(x,y) text(y,25,sprintf('%i',x) ), table( [1:numel(lines_)]', lines_ ) );
    arrayfun( @(x) line([x,x],[-20,30],'color','k'), lines_ );
    mychangept_csminus(i) = inputdlg('Test');
     
end

%% 
clearvars displacement mytable displacements

%%
mychangept = cellfun( @(x) str2double(x), mychangept );
mychangept_csminus = cellfun( @(x) str2double(x), mychangept_csminus );

%%

t=table(csplus_locs', mychangept');
for i = 1:size(t,1) if and(t(i,:).Var2>0,~isnan(t(i,:).Var2)); ingress_csplus(i) = t(i,:).Var1{1}(t(i,:).Var2); end; end;

t=table(csminus_locs', mychangept_csminus');
for i = 1:size(t,1) if and(t(i,:).Var2>0,~isnan(t(i,:).Var2)); ingress_csminus(i) = t(i,:).Var1{1}(t(i,:).Var2); end; end;

%%

figure; subplot(1,2,1); hist(ingress_csplus); subplot(1,2,2); hist(ingress_csminus)

%% CHecking

x = 140
figure; plot( csminus_displacements(x,:) );
line( [ingress_csminus(x),ingress_csminus(x)], [-30,30] )


%%

include_csplus = find( or(ingress_csplus==0,ingress_csplus>3000));
figure; plot( csplus_displacements(include_csplus, : )')

%%

output = rowfun( @(trial,tend) myfourier( displacement_(trial,:), 8000, tend,'k' ), mytable(:,[1,2]) );

%%

