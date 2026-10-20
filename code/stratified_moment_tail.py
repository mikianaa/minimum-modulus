"""Certified tail stratification by retained product bins and their complement."""
from fractional_moment_tail import dyadic_power,fractional_moment,fractional_tail_bound
from verify_minimum_modulus import I,S,ONE,ceildiv,fixed_text


def assemble(full,orders,coefficients,moments,mass_error=0,moment_error=0):
    rows=[];bin_bounds=[];outside_bounds=[]
    for (num,den),coefficient,U in zip(orders,coefficients,moments):
        retained_lower=0;bounds=[]
        for b,u,mass,moment in zip(full.bins,full.upper,full.mass,full.moment):
            if U is not None:
                lower=dyadic_power(I.rat(b),num-den,den).lo*moment//S
                if mass:
                    lower=max(lower,dyadic_power(I.rat(moment,mass),num-den,den).lo*moment//S)
                retained_lower+=lower
            upower=dyadic_power(I.rat(u),num,den).hi
            upper_moment=ceildiv(upower*mass,S)
            if u>b:
                bpower=dyadic_power(I.rat(b),num,den).hi
                aa=max(0,u*mass-moment)
                bb=max(0,moment+moment_error-b*max(0,mass-mass_error))
                upper_moment=min(upper_moment,ceildiv(bpower*aa+upower*bb,(u-b)*S))
            bounds.append(ceildiv(coefficient*upper_moment,S))
        outside=max(0,U.hi-retained_lower) if U is not None else None
        outside_bounds.append(ceildiv(coefficient*outside,S) if outside is not None else None);bin_bounds.append(bounds)
        rows.append(dict(order=f'{num}/{den}',coefficient_upper=fixed_text(coefficient),moment_upper=fixed_text(U.hi) if U is not None else None,retained_moment_lower=fixed_text(retained_lower) if U is not None else None,overflow_tail_upper=fixed_text(outside_bounds[-1]) if outside is not None else None))
    bins_total=sum(min(values) for values in zip(*bin_bounds))
    available=[j for j,v in enumerate(outside_bounds) if v is not None];assert available
    best=min(available,key=lambda j:outside_bounds[j]);outside=outside_bounds[best]
    choices=[min(range(len(orders)),key=lambda j:bin_bounds[j][i]) for i in range(len(full.bins))]
    # Keep a certified moment upper bound for each group of bins assigned
    # to one order. This allows a COMMON future distortion to be changed
    # later without reconstructing the finite distribution.
    grouped=[]
    for j,(order,coefficient) in enumerate(zip(orders,coefficients)):
        cost=sum(bin_bounds[j][i] for i,c in enumerate(choices) if c==j)
        assert coefficient>0
        grouped.append(dict(order=list(order),tail_upper=fixed_text(cost),
                            moment_upper=fixed_text(ceildiv(cost*S,coefficient))))
    state=dict(bins=list(full.bins),upper=list(full.upper),mass=list(map(str,full.mass)),moment=list(map(str,full.moment)),mass_error=str(mass_error),moment_error=str(moment_error))
    return bins_total+outside,dict(stratified_retained_tail_upper=fixed_text(bins_total),stratified_overflow_tail_upper=fixed_text(outside),stratified_overflow_order=rows[best]['order'],stratified_bin_choices=[choices.count(j) for j in range(len(orders))],stratified_orders=rows,stratified_selected_groups=grouped,stratified_bin_order_indices=choices,stratified_state=state,stratified_mass_rounding_error=fixed_text(mass_error),stratified_moment_rounding_error=fixed_text(moment_error))


def tail_bound(full,primes,start,ds,orders,D,mean_upper,overflow_orders=None):
    assert orders and all(isinstance(a,int) and isinstance(b,int) and b in (4,8,16,32) and b<a<=5*b for a,b in orders)
    n=start+len(ds);coefficients=[];moments=[]
    if overflow_orders is None:overflow_orders=orders
    assert overflow_orders and all(isinstance(a,int) and isinstance(b,int) and b>0 for a,b in overflow_orders)
    for num,den in orders:
        coefficients.append(fractional_tail_bound(primes[n-1],n,ONE,D,num,den)[0])
        use=any(num*b==a*den for a,b in overflow_orders)
        moments.append(fractional_moment(primes,start,ds,num,den) if use else None)
    # At most 109 additions per bin per prime, each rounded by < 1 ulp.
    # Probability transitions contract total error; moment transitions
    # amplify it by at most the full Euler first moment.
    mass_error=n*len(full.bins)*109
    moment_error=ceildiv(mass_error*mean_upper,S)
    return assemble(full,orders,coefficients,moments,mass_error,moment_error)
