from hashlib import sha1
from math import log


def get_bits(number, bits):
    """
    Count the leading 1 bits in the value 'h'
    args:
    - number: the value to check.
    - bits: the bits to shift 'number' to the right
            before counting. These number of bits
            were used to create the HLL precision.
    """
    c = 1

    # NOTE: when converting to a binary string,
    # we chop off the first 3 characters since
    # we know that the first 2 are '0b' and the
    # third is a '1'.
    for bit in bin(number >> (bits - 1))[3:]:
        if bit == '1':
            c += 1
            continue
        break

    return c


def _get_alpha(bits):
    """
    'bits' is used for the exponent when calculating
    bucket size, where:

        number of buckets = 2 ^ bits

    """
    bits = int(bits)
    if bits < 4:
        raise Exception("value needs to be > 4" % bits)

    if bits == 4:
        return 0.673

    if bits == 5:
        return 0.697

    if bits == 6:
        return 0.709

    return 0.7213 / (1.0 + 1.079 / (1 << bits))


class HLL(object):
    def __init__(self, bits):
        bits = int(bits)
        assert bits > 4, 'HLL instance bucket count exponent > 4'

        # Bits of precision
        self._bits = bits

        # Number of HLL buckets.
        self.num_buckets = float(2 ** bits)

        self.buckets = [0] * (2 ** bits)

        # Select a constant, defined in the original paper.
        self.alpha = _get_alpha(bits)

    def get_bucket(self, _hash):
        return _hash & (1 << self._bits) - 1

    def get_bits(self, _hash):
        return get_bits(_hash, self._bits)

    def add(self, value):
        """
        Add a value to the hll set.

        """
        _hash = long(sha1(value).hexdigest(), 16)
        bucket = self.get_bucket(_hash)
        self.buckets[bucket] = max(self.get_bits(_hash), self.buckets[bucket])

    def _validate(self, hll):
        """
        Check the precision and bits and alpha of an HLL instance.
        This is to insure it's been set with the same parameters
        before performing union or intersection operations.

        """
        alpha = hll.alpha == self.alpha
        bits = hll._bits == self._bits
        bucket_size = hll.num_buckets == self.num_buckets

        return alpha and bits and bucket_size

    def _estimate(self, buckets):
        numerator = self.alpha * (self.num_buckets ** 2)
        return int(numerator / sum(pow(2.0, -x) for x in buckets))

    def _zeros_estimate(self, num_zeros, buckets):
        E = self.num_buckets * log(self.num_buckets / float(num_zeros))
        if E > 5 * self.num_buckets:
            # TODO: Correct bias
            print 'not correcting bias yet'

        return int(E)

    def estimate(self):
        """
        Estimate the cardinality of values added.

        """
        num_zeros = self.buckets.count(0)
        if not num_zeros:
            return self._estimate(self.buckets)

        return self._zeros_estimate(num_zeros, self.buckets)

    def union(self, *hlls):
        """
        Estimate the union of other HLL instances.

        """
        for hll in hlls:
            if not self._validate(hll):
                raise Exception('HLL instance parameters do not match')

        all_hlls = hlls + (self,)
        all_buckets = [getattr(_, 'buckets') for _ in all_hlls]
        maximum = map(max, zip(*all_buckets))
        num_zeros = maximum.count(0)
        if not num_zeros:
            return self._estimate(maximum)
        return self._zeros_estimate(num_zeros, maximum)


__all__ = ['HLL']
