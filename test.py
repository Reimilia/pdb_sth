import matplotlib.pyplot as plt
import numpy as np
import vector_gen



#This shows how to use

'''
basically, when u run like that, the file will all generated.
when using get_one method, both vector and filename will return back to u
'''
generator = vector_gen.vector_generator('examplePDB.pdb')

list, filename = generator.get_one()

print list
print filename

count = sum(list)
print count


# this is for debugging, forget about that
'''
box = np.array(list).reshape([20,20,20])
scatter = [[i,j,k] for i in range(0,20) for j in range(0,20) for k in range(0,20) if box[i,j,k]!=0]
boxes = np.array(scatter)
boxes = boxes + generator.middles[data_set]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
xs = boxes[:,0]
ys = boxes[:,1]
zs = boxes[:,2]
ax.scatter(xs, ys, zs )
plt.show()
'''