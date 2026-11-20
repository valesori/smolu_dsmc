import numpy as np

from smolu_dsmc.population import Population, UniformStream


def check_invariants(pop: Population):
    masses = pop.masses()
    poly = sorted(pop._poly)
    assert poly == sorted(np.flatnonzero(masses > 1).tolist())
    for k, pos in enumerate(pop._poly):
        assert pop._poly_index[pos] == k
    assert all(pop._poly_index[p] == -1 for p in np.flatnonzero(masses == 1))


def test_random_operations_keep_invariants():
    rng = np.random.default_rng(0)
    uniform = UniformStream(rng, block_size=97)
    pop = Population.monomers(200)
    for _ in range(5000):
        if pop.n > 1 and (pop.n_poly == 0 or uniform() < 0.5):
            pop.merge(*pop.random_pair(uniform))
        else:
            pos = pop.random_polymer(uniform)
            mk = pop.mass(pos)
            m1 = 1 + uniform.index(mk - 1)
            pop.split(pos, m1, mk - m1)
        assert pop.total_mass == 200
    check_invariants(pop)


def test_random_pair_is_uniform_over_ordered_pairs():
    uniform = UniformStream(np.random.default_rng(1))
    pop = Population([1, 2, 3])
    counts = {}
    for _ in range(60_000):
        pair = pop.random_pair(uniform)
        counts[pair] = counts.get(pair, 0) + 1
    assert len(counts) == 6
    assert all(abs(c - 10_000) < 500 for c in counts.values())


def test_histogram_and_zero_masses():
    pop = Population([0, 1, 3, 1, 0, 3, 2])
    assert pop.n == 5
    assert pop.histogram().tolist() == [[1, 2], [2, 1], [3, 2]]
