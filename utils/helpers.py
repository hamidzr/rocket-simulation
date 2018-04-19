import logging, os
import numpy as np

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
def isreal(arr):
  # if input is not a list
  return np.isreal(arr)
  # if inp is a list:
  # return all(it == True for it in np.isreal(2))
