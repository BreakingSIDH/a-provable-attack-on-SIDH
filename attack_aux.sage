def find_smallest_extension(q):
    """
    Computes the smallest integer k such
    that q | # E / F_p^2k
    Returns the int k and the order n
    """
    if not q.is_prime():
        raise NotImplementedError
    # Fp2-Frobenius equals [-p], hence Fp(2k)-Frobenius equals [(-p)^k],
    # hence trace is 2(-p)^k, hence #E(Fp(2k)) is (p^k - (-1)^k)^2 which
    # is divisible by q iff (p^k) ≡ (-1)^k mod q, equivalently (-p)^k ≡ 1
    k = Mod(-p,q).multiplicative_order()
    # assert EA.cardinality(extension_degree=k) % q == 0
    # assert min(k for k,n in enumerate(EA.count_points(k+5),1) if n%q==0) == k
    n = (p^k - (-1)^k)^2
    # assert EA.cardinality(extension_degree=k) == n
    return k, n

def gen_point_order_q(E, q):
    """
    Given a prime q, we can find the minimum
    k such that there is a point of order q.
    Then, we can generate such a point by
    computing the extension Fpk and looking at
    E / Fpk
    To find a point of order q, we generate
    a random point on the curve and multiply by
    (n // q^2)
    If P != E(0) then it must be a point of
    order q.
    """
    if not q.is_prime():
        raise NotImplementedError

    k, n = find_smallest_extension(q)
    Fp2k.<T> = Fp2.extension(k)
    E_new = E.base_extend(Fp2k)
    E_new.set_order(n, num_checks=0)
    infty = E_new(0)

    assert n.is_square() and n%q^2 == 0
    scalar = (n // q^2).isqrt()
    while True:
        P = scalar * E_new.random_point()
        if P != infty:
            break
    P.set_order(q)
    return P

def kerpoly_from_point(P):
    """
    Given a point P, generates a kernel polynomial for Kohel's algorithm
    """
    q = P.order()
    xs = []
    Q = P
    for i in range(q//2):
        xs.append(Q.xy()[0])
        Q += P
    assert len(set(xs)) == len(xs)
    R.<X> = P.curve().base_ring()[]
    poly = prod(X-x for x in xs)
    return poly.change_ring(Fp2)
    
def divide_by(pt, d, all=False):
    """
    Finds all points P s.t. d*P == pt
    If all == False, it returns only one point
    """
    d = ZZ(d)
    if all:
        pt = [pt]
        fun = lambda l, t: sum((u.division_points(l) for u in t), [])
    else:
        fun = lambda l, t: choice(t.division_points(l))
    for l,m in d.factor():
        for _ in range(m):
            pt = fun(l, pt)
    return pt

def find_secret_key(T, P, Q, ea):
    """
    Finds k s.t. <P+k*Q> = <T>, where T,P,Q are 3^ea-torsion points
    This is needed to retrieve the secret key in SIDH once the dual 
    of the secret (partial)-isogeny is available 
    """
    a, b = 0, 0
    c = 3^(ea-1)
    m3iP, m3iQ = c*P, c*Q
    for i in reversed(range(ea)):     
        m3iT  = (3^i)*T
        m3iaP = (3^i*a)*P
        m3ibQ = (3^i*b)*Q
        RHS = m3iT - m3iaP - m3ibQ
        for s in range(3):
            for r in range(3):
                if s*m3iP+r*m3iQ == RHS:
                    a = a + s*3^(ea-1-i)
                    b = b + r*3^(ea-1-i)   
    assert a*P+b*Q == T
    return ( inverse_mod(a,3^ea) * b ) % 3^ea


def find_secret_key_pairing(T, P, Q, ea):
    """
    Similar to the above function but uses Weil-pairing
    """

    # e(P,Q)
    pair_PQ = P.weil_pairing(Q, 3^ea)

    pair_a = T.weil_pairing(Q, 3^ea)
    a = pair_a.log(pair_PQ)

    # e( alpha(aP + bQ), P ) = e(P, Q)^-bx
    pair_b = T.weil_pairing(P, 3^ea)
    b = - pair_b.log(pair_PQ) % 3^ea

    return ( inverse_mod(a,3^ea) * b ) % 3^ea
