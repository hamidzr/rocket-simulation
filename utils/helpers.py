import json
import logging, os
from functools import wraps
import numpy as np
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

def plot_attempts(values, fname, xlabel='Try#', ylabel=None, title=None):
  tries = np.linspace(1,len(values), len(values), dtype=int)
  plt.scatter(tries, values, c='g')
  plt.plot(tries, values)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  plt.savefig('figs/' + fname)
  plt.close()


def dump(fname, data):
  with open(fname, 'w', encoding='UTF-8') as f:
    json.dump(data, f)

def load(fname):
  with open(fname, 'r', encoding='UTF-8') as f:
    return json.load(f)
