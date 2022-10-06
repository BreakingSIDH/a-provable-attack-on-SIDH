def generate_torsion_points(E, la, ea, lb, eb):
    def get_l_torsion_basis(E, l):
        n = (p+1) // l
        return (n*G for G in E.gens())
    P2, Q2 = get_l_torsion_basis(E, la^ea)
    P3, Q3 = get_l_torsion_basis(E, lb^eb)
    return P2, Q2, P3, Q3

def check_torsion_points(E, la, ea, lb, eb, P2, Q2, P3, Q3):
    # Make sure Torsion points are
    # generated correctly
    return all([la^(ea-1)*P2 != E(0),
                lb^(eb-1)*P3 != E(0),
                P2.weil_pairing(Q2, la^ea)^(la^(ea-1)) != 1,
                P3.weil_pairing(Q3, lb^eb)^(lb^(eb-1)) != 1])

def gen_keypair(E, P, Q, P_other, Q_other, l, e):
    # generate challenge key
    key = randint(0,l^e)
    K = P + key*Q
    phi = E.isogeny(K, algorithm="factored")
    E_new = phi.codomain()
    E_new.set_order((p+1)^2, num_checks=0)
    P_new, Q_new = phi(P_other), phi(Q_other)
    return key, (E_new, P_new, Q_new)
