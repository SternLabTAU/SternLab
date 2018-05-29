#! /usr/local/python_anaconda/bin/python3.4


def frange(start, stop, step):
    x = start
    while x < stop:
        yield round(x, 1)
        x += step