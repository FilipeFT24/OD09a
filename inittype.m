function [Fd] = inittype(type, f, xdc, xqc, Wbf_m)
switch type
    case 0 % INTERPOLATION (NODAL!):
        Fd = f(xdc);
    case 1 % PROJECTION:
        Fq = f(xqc);
        Fd = Fq*Wbf_m;
    otherwise
        return
end
end