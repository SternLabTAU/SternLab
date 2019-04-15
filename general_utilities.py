#! /usr/local/python_anaconda/bin/python3.4


def frange(start, stop, step):
    x = start
    while x < stop:
        yield round(x, 1)
        x += step


def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z