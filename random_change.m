function [changed_prec_dat,changed_imp_dat] = random_change(p_dat,i_dat)

    changed_prec_dat = p_dat;
    changed_imp_dat = i_dat;
    s = size(changed_imp_dat) ;
    
    un = unifrnd (zeros(s(1),1),changed_imp_dat(:,3));

    changed_prec_dat(:,3) = changed_prec_dat(:,3) - un ;
    changed_imp_dat(:,3) = changed_imp_dat(:,3) - un ;
    changed_imp_dat(:,4) = changed_imp_dat(:,4) - un ;

end
