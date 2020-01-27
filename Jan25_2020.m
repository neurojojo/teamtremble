%% Fourier (with boxplot)
%% 

% The variables ingress_csplus and ingress_csminus contain the indices of
% the ingress (obtained from hand selecting within Pia_Jan23_2020.m)

N = 19;
trials = 10;
t0 = 2000;
tend = 3000;

animals_lbl = cell2mat(arrayfun(@(x) repmat(x,trials,1), [1:N], 'UniformOutput', false) );
animals_lbl = animals_lbl(:);

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

figure; plot( csplus_table.csplus_max, csplus_table.csplus_coeffs, 'r.' )
hold on; plot( csminus_table.csminus_max, csminus_table.csminus_coeffs, 'k.' )

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

%% Compare minus and plus

%csplus_table = csplus_table(include_csplus,:);
%csminus_table = csminus_table(include_csminus,:);

csplus_by_mouse = arrayfun(@(x) nanmean(csplus_table(csplus_table.Mouse==x,:).csplus_coeffs), [1:N] );
csminus_by_mouse = arrayfun(@(x) nanmean(csminus_table(csminus_table.Mouse==x,:).csminus_coeffs), [1:N] );

figure; arrayfun( @(x) line([1,2],[csminus_by_mouse(x),csplus_by_mouse(x)]), [1:N] )