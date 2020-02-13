#!/usr/bin/env python
# Copyright (C) 2017, Carlo de Franchis <carlo.de-franchis@polytechnique.org>

import os
import sys
import traceback
import multiprocessing

from s2p import common
from s2p.config import cfg


def show_progress(a):
    """
    Print the number of tiles that have been processed.

    Args:
        a: useless argument, but since this function is used as a callback by
            apply_async, it has to take one argument.
    """
    show_progress.counter += 1
    status = "done {:{fill}{width}} / {} tiles".format(show_progress.counter,
                                                       show_progress.total,
                                                       fill='',
                                                       width=len(str(show_progress.total)))
    if show_progress.counter < show_progress.total:
        status += chr(8) * len(status)
    else:
        status += '\n'
    sys.stdout.write(status)
    sys.stdout.flush()


def tilewise_wrapper(fun, *args, **kwargs):
    """
    """
    if not cfg['debug']:  # redirect stdout and stderr to log file
        f = open(kwargs['stdout'], 'a')
        sys.stdout = f
        sys.stderr = f

    try:
        out = fun(*args)
    except Exception:
        print("Exception in %s" % fun.__name__)
        traceback.print_exc()
        raise

    common.garbage_cleanup()
    if not cfg['debug']:  # close logs
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        f.close()

    return out


def launch_calls(fun, list_of_args, nb_workers, *extra_args, tilewise=True,
                 timeout=600):
    """
    Run a function several times in parallel with different given inputs.

    Args:
        fun: function to be called several times in parallel.
        list_of_args: list of (first positional) arguments passed to fun, one
            per call
        nb_workers: number of calls run simultaneously
        extra_args (optional): tuple containing extra arguments to be passed to
            fun (same value for all calls)
        tilewise (bool): whether the calls are run tilewise or not
        timeout (int): timeout for each function call (in seconds)

    Return:
        list of outputs
    """
    results = []
    outputs = []
    show_progress.counter = 0
    show_progress.total = len(list_of_args)
    pool = multiprocessing.Pool(nb_workers)
    for x in list_of_args:
        args = tuple()
        if type(x) == tuple:
            args += x
        else:
            args += (x,)
        args += extra_args
        if tilewise:
            if type(x) == tuple:  # we expect x = (tile_dictionary, pair_id)
                log = os.path.join(x[0]['dir'], 'pair_%d' % x[1], 'stdout.log')
            else:  # we expect x = tile_dictionary
                log = os.path.join(x['dir'], 'stdout.log')
            args = (fun,) + args
            results.append(pool.apply_async(tilewise_wrapper, args=args,
                                            kwds={'stdout': log},
                                            callback=show_progress))
        else:
            results.append(pool.apply_async(fun, args=args, callback=show_progress))

    for r in results:
        try:
            outputs.append(r.get(timeout))
        except KeyboardInterrupt:
            pool.terminate()
            sys.exit(1)

    pool.close()
    pool.join()
    common.print_elapsed_time()
    return outputs
