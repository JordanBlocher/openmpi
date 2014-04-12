#!bin/python
from random import seed, randint, uniform

seed(0)

print "20"
for i in range(0,20):
    print uniform(1, 200), " ", uniform(1, 200), " ", uniform(1, 200), randint(1,5)
