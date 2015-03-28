`hllpy`
------

HyperLogLog(++) draft in Python

Installation:

    $ python setup.py install

Example:

    from hllpy import HLL

    UNIQUES = 50000
    hll = HLL(16)
    # Add 50K unique elements.
    for i in xrange(UNIQUES):
        hll.add(str(i))

    # Is the set size estimate less than smaller than 1% error?
    assert abs(1 - (hll.estimate() / float(UNIQUES))) < .01

References:

[https://github.com/armon/hlld](https://github.com/armon/hlld)

[https://github.com/svpcom/hyperloglog](https://github.com/svpcom/hyperloglog)

[http://tech.adroll.com/blog/data/2013/07/10/hll-minhash.html](http://tech.adroll.com/blog/data/2013/07/10/hll-minhash.html)

[http://research.neustar.biz/2012/10/25/sketch-of-the-day-hyperloglog-cornerstone-of-a-big-data-infrastructure/](http://research.neustar.biz/2012/10/25/sketch-of-the-day-hyperloglog-cornerstone-of-a-big-data-infrastructure/)

-------------------------------------------------------------------


    TODO:
    - better intersections?
