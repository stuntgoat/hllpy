from hllpy import HLL


class TestUnion(object):
    def test_union(self):

        h2 = HLL(16)
        h3 = HLL(16)
        h4 = HLL(16)
        h5 = HLL(16)

        for i in xrange(50000):
            h2.add(str(i))

        for i in xrange(50000):
            h3.add(str(i) + 'hi')

        for i in xrange(50000):
            h4.add(str(i) + 'bye')

        for i in xrange(50000):
            h5.add(str(i) + 'happy')

        estimate = h2.union(h3, h4, h5)
        diff = estimate / 200000.0

        # Check that the union estimate is < 1%
        assert abs(1 - diff) < .009
