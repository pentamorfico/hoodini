from flagcsnap.utils import console,flat,desaturate
import shutil,os

class Cleaner:
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self):       
        if self.keep is None:
            if os.path.exists(self.output+"/assembly_folder"):
                shutil.rmtree(self.output+"/assembly_folder")
