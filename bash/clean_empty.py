result_PREFIX = '/n/scratch2/yw174/result/'
import os

if __name__=='__main__':
    real_dir = os.path.join(result_PREFIX,'fast')
    filedir = os.listdir(real_dir)
    for files in filedir:
        dirs = os.path.join(real_dir,files)
        if os.path.getsize(dirs) < 50:
            os.system('rm '+ dirs)
            print files + ' has been wiped up'