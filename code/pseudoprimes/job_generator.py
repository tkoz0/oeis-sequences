# print out job arguments
# job size of 2**30 up to 2**42

print('lo,hi')
for lo in range(0,2**42,2**30):
    hi = lo+2**30-1
    print('%d,%d'%(max(2,lo),hi))
