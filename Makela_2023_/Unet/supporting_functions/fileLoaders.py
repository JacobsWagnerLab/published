import os
import re

def getFileList(path_to_file, chan_name = '', first_only = False, file_ending = False):
    """supporting function to walk in the experiment directory and return a list 
    of the all files existing in the master directory. 
    
    first_only - choose first image only
    file_ending - only files that ends with a specific string
    """
            
    file_list = []
    for (root,_,files) in os.walk(path_to_file):
        # load only the first image
        if first_only:
            if root.endswith(chan_name):
                name = sorted(files)[0]
                file_list.append(os.path.join(root,name))
        
        # if user has specific file ending
        if file_ending != False:
            for name in files:
                if re.search(chan_name,root) is not None:
                    if name.endswith(file_ending):
                        file_list.append(os.path.join(root,name))
            
        # or all images
        else:
            for name in files:
                if re.search(chan_name,root) is not None:
                    if name.endswith(('.tiff','.tif')):
                        file_list.append(os.path.join(root,name))
    
    print("{} files were loaded.".format(len(file_list)))

    return file_list 