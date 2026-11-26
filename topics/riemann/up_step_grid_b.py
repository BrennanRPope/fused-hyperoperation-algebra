import numpy as np, time, pickle
from multiprocessing import Pool
from up_step_grid import cell
if __name__=='__main__':
    jobs=[(x,y,1) for x in np.arange(0,60.01,1.0) for y in np.arange(21,46.01,1.0)]
    jobs+=[(x,y,2) for x in np.arange(0,60.01,1.0) for y in np.arange(0,46.01,1.0)]
    t=time.time()
    with Pool(4) as p: res=p.map(cell,jobs,chunksize=8)
    pickle.dump(res,open('up_step_grid_b.pkl','wb')); print('done',time.time()-t)
