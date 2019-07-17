function [Uijs,Uis,Ujs,order] = HCCA_PPI(Xs, co, perc, Ls)
    
    Uijs = {};
    Uis = {};
    Ujs = {};
    order = [];

    imin = -1;
    jmin = -1;
    m = inf;
    for i=1:length(Xs)-1
        if isempty(Xs{i})
            continue
        end
        N = size(Xs{i},1);
        Xi = normc(Xs{i});
        for j=i+1:length(Xs)
            if isempty(Xs{j})
                continue
            end
            Xj = normc(Xs{j});
            c = cond(Xi'*Xj);
            if c < m
                m = c;
                imin = i;
                jmin = j;
            end
        end
    end
    if imin == -1
        return
    end
    
    Li = [];
    Lj = [];
    if imin <= length(Ls)
        Li = Ls{imin};
    end
    if jmin <= length(Ls)
        Lj = Ls{jmin};
    end
    
    order = [imin;jmin];
    [Ui,Uj,~,~,~] = MyCCAPPI(Xs{imin}, Xs{jmin}, co, perc, Li, Lj);
    Uij = [Ui Uj];
    Uijs = [Uijs; Uij];
    Uis = [Uis; Ui];
    Ujs = [Ujs; Uj];
    
    Xs{imin} = [];
    Xs{jmin} = [];
    Xs = [Xs; Uij];
    
    [Uijs2,Uis2,Ujs2,order2] = HCCA_PPI(Xs, co, perc, Ls);
    
    Uijs = [Uijs; Uijs2];
    Uis = [Uis; Uis2];
    Ujs = [Ujs; Ujs2];
    order = [order;order2];
    
end
