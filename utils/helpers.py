import logging, os
import pickle
from functools import wraps
import numpy as np
import matplotlib
if not os.getenv('DISPLAY'): matplotlib.use('agg') # to have it work in remote machines
import matplotlib.pyplot as plt

logging.basicConfig(level=getattr(logging, os.getenv('DEBUG', 'INFO')))
logger = logging.getLogger('rocket')

def filter(pred, lst):
  # Filters a list and returns the matching indices
  matches = []
  for index in range(len(lst)-1, -1, -1):
    if pred(lst[index]):
      matches.append(index)
  return matches

# mimic matlab's isreal()
def isreal(inp):
  # if input is not a list
  res = np.isreal(inp)
  return res
  # if inp is a list:
  # return all(it == True for it in np.isreal(2))

# memoize helper
def memoize(f):
    cache = {}
    @wraps(f)
    def decorated(*args):
        key = (f, str(args))
        result = cache.get(key, None)
        if result is None:
            result = f(*args)
            cache[key] = result
        return result
    return decorated

def plot_attempts(values, fname=None, xlabel='Iteration', ylabel=None, title=None, line_label=None):
  iters = np.linspace(1,len(values), len(values), dtype=int)
  plt.scatter(iters, values, c='green')
  plt.plot(iters, values, c='b', linewidth=1, label=line_label)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  if fname:
    plt.savefig('figs/' + fname)
    plt.close()


# finds the longest sub list
def _longest(l):
    if(not isinstance(l, list)): return(0)
    return(max([len(l),] + [len(subl) for subl in l if isinstance(subl, list)] +
        [_longest(subl) for subl in l]))

def sort_two_lists(list1, list2):
  list1, list2 = zip(*sorted(zip(list1, list2)))
  list1, list2 = (list(t) for t in zip(*sorted(zip(list1, list2))))
  return list1, list2

# plots all of angle attemps or landing altitude attempts
def plot_batch_attempts(values_arr, ratios, fname=None, xlabel='Iteration', ylabel=None, title=None):
  for attempt, values in enumerate(values_arr):
    iters = np.linspace(1,len(values), len(values), dtype=int)
    plt.scatter(iters, values)
    plt.plot(iters, values, linewidth=1, label=f'#{attempt}: {ratios[attempt]:.5}')

  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  plt.legend()
  if fname:
    plt.savefig('figs/' + fname)
    plt.close()


def dump(fname, data):
  with open(fname, 'wb') as f:
    pickle.dump(data, f)

def load(fname):
  with open(fname, 'rb') as f:
    return pickle.load(f)
