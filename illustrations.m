clear all;
close all;


%% parameters
%%% methods
HYBRID = 1;
SIMULATED = 2;
COMBINATORIAL = 3;

%%% constante definition for an automatic domain definition
AUTO = 1;

%%% number of simulations of the simulated annealing method
n_sims = 10 ;

%%% number of sampling positions (used for grids and lines ...)
n_sampling_positions = 10 ;

%%% figure management
nFig = 1;

%% data loading
prec_dat = load('data/data_irsn');
imp_dat = load('data/data_irsn_imprecise');


%% sampling positions definition 
%  if vector nx2: n sampling points of a line
%  if vector 2xn grid sampling points

%%%%%% LINE %%%%%%% (nx2) 
%%%%%%              (it can also be random sampling points) 
%%%%%%              (2 columns: X and Y positions)

% from a file:
line1 = load('data/Line1_positions');
line2 = load('data/Line2_positions');

% from here
X_A = 594.097 ; Y_A = 2418.75 ; X_B = 609.42 ; Y_B = 2419.22;
line1 = line_creation( X_A,Y_A,X_B,Y_B,n_sampling_positions );
X_A = 620.824 ; Y_A = 2394.28 ; X_B = 608.712 ; Y_B = 2439.159;
line2 =  line_creation( X_A,Y_A,X_B,Y_B,n_sampling_positions );

%%%%%% GRID %%%%%%% (2xn)
% from data:
X_liminf = min(prec_dat(:,1))-5*(max(prec_dat(:,1)) - min(prec_dat(:,1)))/100;
X_limsup = max(prec_dat(:,1))+5*(max(prec_dat(:,1)) - min(prec_dat(:,1)))/100;
Y_liminf = min(prec_dat(:,2))-5*(max(prec_dat(:,2)) - min(prec_dat(:,2)))/100;
Y_limsup = max(prec_dat(:,2))+5*(max(prec_dat(:,2)) - min(prec_dat(:,2)))/100;

grid = zeros(2,n_sampling_positions);
for i=1:1:n_sampling_positions
    grid(1,i) = X_liminf + i*(X_limsup - X_liminf)/(n_sampling_positions-1);
    grid(2,i) = Y_liminf + i*(Y_limsup - Y_liminf)/(n_sampling_positions-1);
end

% from a file:
% grid  whether 2 lines of n points respectively representing the X sampling and Y sampling of the grid or take the transpose of a 'line' file...
grid_IRSN = load('data/Grid_IRSN_50');

%% Variogram parameters definition
precise_vario_params = [0.02 0.08 11];
imprecise_vario_params = [0 0.04 ; 0.06 0.1 ; 10 12];

%% imprecise kriging with precise data
if(1)
    tmpImpData = [prec_dat prec_dat(:,3)];
    % intervallist kriging 
    t = cputime; 
	[pred_inf_grid,pred_sup_grid] = intervallist_kriging(tmpImpData,imprecise_vario_params,AUTO,grid_IRSN,HYBRID,n_sims);
	e = cputime-t ;
    fprintf('time for intervallist kriging: %f sec.\n',e);
    
	figure(nFig); nFig = nFig + 1 ; 
	surf(grid_IRSN(1,:),grid_IRSN(2,:),pred_inf_grid);
    hold on;
	surf(grid_IRSN(1,:),grid_IRSN(2,:),pred_sup_grid);
    title('intervallist kriging');
	figure(nFig); nFig = nFig + 1 ; 
	surf(grid_IRSN(1,:),grid_IRSN(2,:),pred_sup_grid-pred_inf_grid);
    title('intervallist kriging imprecision');
end


%% Kriging on a grid
if(0)
    
    % precise kriging 
    t = cputime; 
    [pred_grid,var_grid] = precise_kriging(prec_dat,precise_vario_params,AUTO,grid_IRSN);
    e = cputime-t ;
    fprintf('time for precise kriging (not optimized): %f sec.\n',e);
    
	figure(nFig); nFig = nFig + 1 ;
	surf(grid_IRSN(1,:),grid_IRSN(2,:),pred_grid);
    title('precise kriging');
    figure(nFig); nFig = nFig + 1 ;
	surf(grid_IRSN(1,:),grid_IRSN(2,:),var_grid);
    title('variance of the precise kriging');
    
    % intervallist kriging 
    t = cputime; 
	[pred_inf_grid,pred_sup_grid] = intervallist_kriging(imp_dat,imprecise_vario_params,AUTO,grid_IRSN,HYBRID,n_sims);
	e = cputime-t ;
    fprintf('time for intervallist kriging: %f sec.\n',e);
    
	figure(nFig); nFig = nFig + 1 ; 
	surf(grid_IRSN(1,:),grid_IRSN(2,:),pred_inf_grid);
    hold on;
	surf(grid_IRSN(1,:),grid_IRSN(2,:),pred_sup_grid);
    title('intervallist kriging');
	figure(nFig); nFig = nFig + 1 ; 
	surf(grid_IRSN(1,:),grid_IRSN(2,:),pred_sup_grid-pred_inf_grid);
    title('intervallist kriging imprecision');
end

%% Comparison on a line
if(0)
	[pred_line1,var_line1] = precise_kriging(prec_dat,precise_vario_params,AUTO,line1);
	[pred_inf_line1,pred_sup_line1] = intervallist_kriging(imp_dat,imprecise_vario_params,AUTO,line1,HYBRID,n_sims);

	figure(nFig); nFig = nFig + 1 ; 
	hold on;
	
	ax = plotyy( line1(:,1) , var_line1 , line1(:,1) , pred_sup_line1 - pred_inf_line1);
	xlabel ('X positions');
	ylabel (ax(1),'variance');
	ylabel (ax(2),'imprecision');
	if(1)
		for k=1:10

			[prec,imp] = random_change(prec_dat,imp_dat);
			[pred_line1,var_line1] = precise_kriging(prec,precise_vario_params,AUTO,line1);
			[pred_inf_line1,pred_sup_line1] = intervallist_kriging(imp,imprecise_vario_params,AUTO,line1,HYBRID,n_sims);
			ax = plotyy( line1(:,1) , var_line1 , line1(:,1) , pred_sup_line1 - pred_inf_line1);
			%plot( line1(:,1) , pred_sup_line1 - pred_inf_line1);
			%plot(line1(:,1),10*);
			%plot(,'r','-;change k;');

		end

	end

end


