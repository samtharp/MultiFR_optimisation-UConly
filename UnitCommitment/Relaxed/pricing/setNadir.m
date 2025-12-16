function [DV_Total_FR_atTd,Inertia_term,PLoss_term,FR_term,Bounds,...
    Nadir_constraints] = ...
    setNadir(InputData,FR,H_total,DV_PLoss,Bounds)

% Author: Luis Badesa
% Pricing adaptation (Samuel Tharp):
% - Fix nadir interval i (no implies, no binary selection)
% - Keep the SAME vector SOC form labeled 'nadir' so study_duals.m works

    Td = InputData.Td;
    Pmax = InputData.Gen_limits(:,2);
    PLossMax = InputData.PLossMax;
    nadir_req = InputData.nadir_req;
    FR_capacity = InputData.FR_capacity;
    H_const = InputData.H_const;
    H_const_Wind = InputData.H_const_Wind;
    H_const_Nuclear = InputData.H_const_Nuclear;
    P_Wind = InputData.P_Wind;
    num_gen = InputData.num_gen;

    %% Total FR delivered by each time Td(i)
    DV_Total_FR_atTd = sdpvar(1,length(Td));
    for i = 1:length(Td)
        DV_Total_FR_atTd(i) = 0;
        for j = 1:length(Td)
            if Td(j) > Td(i)
                DV_Total_FR_atTd(i) = DV_Total_FR_atTd(i) + FR(j)*Td(i)/Td(j);
            else
                DV_Total_FR_atTd(i) = DV_Total_FR_atTd(i) + FR(j);
            end
        end
    end
    Bounds = [Bounds,...
              0 <= DV_Total_FR_atTd <= (num_gen'*FR_capacity)*ones(1,length(Td))];
    clear i j

    %% Define LHS/RHS terms for each possible interval i
    Inertia_term = sdpvar(1,length(Td));
    PLoss_term   = sdpvar(1,length(Td));
    FR_term      = sdpvar(1,length(Td));

    for i = 1:length(Td)
        Inertia_term(i) = H_total;
        PLoss_term(i)   = DV_PLoss;

        FR_term_defined = false;
        for j = 1:length(Td)
            if Td(j) < Td(i)
                Inertia_term(i) = Inertia_term(i) - FR(j)*Td(j)/(4*nadir_req);
                PLoss_term(i)   = PLoss_term(i) - FR(j);
            else
                if ~FR_term_defined
                    FR_term(i) = FR(j)/Td(j);
                    FR_term_defined = true;
                else
                    FR_term(i) = FR_term(i) + FR(j)/Td(j);
                end
            end
        end
    end
    clear i j FR_term_defined

    %% Bounds (same style as original)
    UpperBound_FR_term = num_gen'*(FR_capacity./Td);
    UpperBound_FR_term = UpperBound_FR_term*ones(1,length(Td));
    Bounds = [Bounds, 0 <= FR_term <= UpperBound_FR_term];

    UpperBound_Inertia_term = ((num_gen.*H_const)'*Pmax)/InputData.f_0 + ...
        (H_const_Wind*P_Wind)/InputData.f_0 + ...
        (H_const_Nuclear*PLossMax)/InputData.f_0;
    UpperBound_Inertia_term = UpperBound_Inertia_term*ones(1,length(Td));

    LowerBound_Inertia_term = - num_gen'*(FR_capacity.*Td/(4*nadir_req));
    LowerBound_Inertia_term = LowerBound_Inertia_term*ones(1,length(Td));

    LowerBound_PLoss_term = (PLossMax - num_gen'*FR_capacity);
    LowerBound_PLoss_term = LowerBound_PLoss_term*ones(1,length(Td));

    Bounds = [Bounds,...
        LowerBound_Inertia_term <= Inertia_term <= UpperBound_Inertia_term,...
        LowerBound_PLoss_term   <= PLoss_term   <= PLossMax*ones(1,length(Td))];

   %% Pricing nadir enforcement (fixed interval, no implies)
i_fix = 1;

t = Inertia_term(i_fix) + FR_term(i_fix);
u = [Inertia_term(i_fix) - FR_term(i_fix); ...
     2*sqrt(1/(4*nadir_req))*PLoss_term(i_fix)];

% Single SOC constraint in YALMIP cone primitive
c_nadir = cone(u, t);

% Label it so Constraints('nadir') returns THIS object
Nadir_constraints = [c_nadir : 'nadir'];

% (Optional but recommended for the rotated-SOC derivation assumptions)
Nadir_constraints = [Nadir_constraints, Inertia_term(i_fix) >= 0, FR_term(i_fix) >= 0];

end

