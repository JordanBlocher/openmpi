#!bin/python
from random import seed, randint, uniform, random

seed(0)

print "20"
for i in range(0,20):
    print uniform(1, 200), " ", uniform(1, 200), " ", uniform(1, 200), " ", random(), " ", random(), " ", random(), " ", randint(1,5)
